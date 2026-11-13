# hifiasm 断点续跑 + 性能优化 — Handoff

## 项目背景

- 上游：https://github.com/chhylp123/hifiasm，基线 commit `ec9a8b2`
- 本地仓：`~/Documents/projects/tcuni/hifiasm`（分支 master，工作区干净）
- 动机：长任务（真实数据数天~两周）中途断电/被杀/OOM 得从零开始 → 加了阶段级 checkpoint
- 后续延伸：顺带做了可测的性能优化与全链路性能剖析
- 版本：0.25.0-r726 + 本分支 7 个 commit

## 已交付

### 1. 断点续跑（commit 4ad0d09 → 08b9525）

**机制**：数据库 WAL 思路的移植
- `*.r<N>.ecstr`：纠错每轮的 per-read WAL（worker 每完成一条 read 追加一条 `{magic,rid,plen,payload,crc32}`）
- `*.r<N>.snap`：轮级快照（tmp 写→fsync→rename 原子提交；全文件 crc + FIN 哨兵；fmt v2 含运行时变动的 `max_n_chain`）
- `*.ff.ovstr`：最终 overlap 阶段同款 WAL
- 参数指纹 `params_hash`：-k/-w/输入文件等改动 → 旧 checkpoint 拒绝混用
- 幂等边界：只固化 `cal_ec_multiple`/`cal_ov_r` 已提交结果；`sl_ec_r`/`post_rev` 等非幂等 pass 一律重放
- `--ckpt-int N`（默认 60s fsync 间隔；0=关闭）；`-i` 旁路一切
- `ha_ckpt_purge()`：成功后清理 snap/ecstr + 半成品 `.tmp`

**已验证场景**（全部 GFA bit-exact）：round 内 SIGKILL、resume 中再 SIGKILL（嵌套）、ecstr 尾部垃圾、snap 中间翻字节、-k 变更拒绝、-i 旁路、ff 阶段 kill、snap 写入中 kill。

**实测开销（9036, 1Gb reads, -t16）**：214.7s vs 212.3s ≈ **+1.1%**；内存增量 ≈0。

### 2. 性能优化（commit 344599a）

`htab.cpp` k-mer/minimizer 计数流水线：step0 的串行 ingest（`recover_UC_Read`/`ha_compress_base`/N位点/名字拷贝）下沉到 step1 的 `kt_for` worker 并行；w==1 枚举用每线程私有桶行。

**实测（9036, 1Gb reads）**：
- -t16: 214.7→200.8s（**-6.5%**），-t64: 117.2→107.4s（**-8.3%**）
- 5 个 GFA 与基线逐字节一致；峰值内存不变（17.2GB）
- 性能改动后 kill/resume 复验通过（426Mb reads：56s 处 SIGKILL → resume 36s，WAL 回放 33793/38743，GFA bit-exact）

## 测试机与实测数据

| 机器 | 规格 | 状态 |
|---|---|---|
| 9036 | 2×EPYC 7402, 96t, 128GB | 测试机，仓在 `~/hifiasm-ckpt/`，数据 `t/`（1Gb）、`q/`（各次跑批） |
| 9009 | 同规格 | **生产机，共享使用**（minio/jenkins/docker 常驻，load ~16），只做只读观察 |

线程伸缩实测（9036）：-t16=212s → -t32=143s → -t64=117s → -t96=110s。

## 性能剖析结论（全部实测，非推测）

阶段分布（-t16, ~200s）：纠错 117.5s(58%) | 读入+计数 21s | 建图+输出 24s | 索引 13s | 其余 ~25s

### 试过且无效的（别再试）

- `-march=native`（AVX2/BMI2）：201.9s ≈ 持平
- jemalloc LD_PRELOAD：214.3s，**-7% 更慢**
- `-flto`：203.5s ≈ 持平
- NUMA 钉 socket0 -t48：124.3s，不如散跑 -t64

### 不能安全并行的（查实原因）

- `normalize_ma_hit_t_single_side_advance`×2(1.6s)：双向写 `del` 标志，并行→写序不确定→输出非确定。上游 `_mult` 并行版被注释就是这个原因
- `output_trio_unitig_graph` FATHER/MOTHER(4.8s)：共享写 `ruIndex->index`/`sg->seq[].c`
- 读入+计数：内存带宽瓶颈非 CPU(-t64 与 -t16 同速）

### 关键硬件发现（30% 级杠杆，唯一被数据支持的大头）

**9036 和 9009 都有同一个硬件缺陷：socket 1 没插内存**

```
numactl --hardware:
  node0: CPU 0-23,48-71 + 128GB
  node1: CPU 24-47,72-95 + 0 GB   ← socket1 全远端访问（距离32 vs 10）
```

含义：-t96 时一半线程每次内存访问跨 UPI 走 socket0。这是 -t64→-t96 收益趋零、计数阶段"带宽饱和"的根因。**修复=给 socket1 配 DIMM：带宽翻倍+全本地化，预期 20-40%**。

## 生产跑法建议（9009）

共享机（load ~16 常驻 daemon），建议：

```bash
numactl --cpunodebind=0 ./hifiasm -o <prefix> -t48 <reads>   # 钉 socket0，全本地访问，socket1 留给 daemon
```

9009 上**不要**用 -t96 裸跑（会与 daemon 抢 socket1 远端带宽+超额订阅）。错峰跑更好。

## 复测/交接操作

```bash
# 构建
cd hifiasm && make -j8

# 断点续跑验证（kill/resume harness）
cd test_ckpt && ./kill_resume_test.sh ../hifiasm <reads.fa> <outdir> 16 [ecstr阈值]

# 9036 上的既有产物
~/hifiasm-ckpt/t/reads_1g.fa   # 1Gb 测试集
~/hifiasm-ckpt/q/q1..q12.*     # 各轮跑批日志+输出
```

## 遗留事项

- [ ] **硬件工单**：9009/9036 socket1 配内存（唯一 30% 级机会）
- [ ] graph 阶段（建图+输出，~24s）无阶段内 checkpoint——崩溃整段重放；`all.debug.*.bin` 机制可升级复用
- [ ] resume 每轮重建 `ha_pt_gen` 索引（~4-9s/轮）；可选把索引也落盘（常态 IO 税换小概率收益，倾向不做）
- [ ] PGO 未测完（插桩慢 10 倍，预期 <10%，性价比低）
- [ ] 纠错 58% 的剩余加速空间在算法层（chain/align 向量化、跳已收敛 read），不在调度层

## 注意

- 所有结论均为实测；不要再假设未验证的加速比
- `ec.bin` 在多轮正常跑之间字节也不同（非确定性），正确性判据=GFA bit-exact
- 9009 是生产机：只观察，不要拷文件、不要跑测试

---

## 多机分片可行性验证（2026 最新，已验证）

**结论：纠错阶段按 read 分片跨机可行，集合级无损；但 read 内条目顺序是调度依赖的，合并需要规范化排序。**

实验（9036，426Mb reads，-t16 -r1，`HIFIASM_SHARD=k/N` + `HIFIASM_PAF_DUMP` 探针）:

| 检查 | 结果 |
|---|---|
| shard0 + shard1 PAF 合并 vs 完整跑 | **817,984 行，集合逐字节相等** |
| 两次同配置完整跑互比 | 顺序 0 差异（同配置确定） |
| shard vs 完整：同一 read 内条目顺序 | 98% 的 read 有差异（worker 线程缓冲调度依赖，原二进制固有，跨 -t 也存在；GFA 不受影响） |

含义/工程要求：
- 分片产出**正确的结果集**，但原始 PAF 不能 byte-diff；合并工具需**按 read 分组 + 组内规范排序**(sort by tn/coords)，或直接用最终 GFA 作验收
- index 各机本地重建已验证等价
- 预期收益：纠错占 58% 可 50/50 切分，2 机省 ~25-30%，3 机 ~40%（待端到端实测）
- 未验证项：多轮同步（每轮后需交换修正序列）、index 跨机构建逐位一致性、真实双机端到端耗时

## ckpt 代码评审问题处理（本轮已修）

| 级别 | 问题 | 处理 |
|---|---|---|
| 中 | `ckpt_read_one` 失败路径泄漏/半成品状态（tmp 不释放、hom/het_cov 提前覆盖、max_n_chain 先写后验） | **已修**：全部加载入局部变量，CRC 全量验证后才提交全局；`free_partial_reads`/`free_cc_v`/`free_paf_arr` 清理；指针数组改 calloc 保证部分释放安全 |
| 中 | fwrite/fflush/fsync 返回值不检查，磁盘满会把截断文件 rename 成正文件 | **已修**:`ckpt_rename_sync`/`write_All_reads`/`write_ma_hit_ts` 加 ferror+fflush+fsync+fclose 检查，失败不 rename、不打印成功 |
| 中 | flt_tab 正文不在 CRC 内，中段位翻转静默出错 | **已修**:open_memstream 缓冲+段级 CRC+长度前缀；且整体改为**两遍读**(pass1 先验 footer+全体 CRC,pass2 才解析）——损坏文件到不了任何解析器，**FMT_VER 2→3**（旧快照自动拒绝） |
| 低 | `(uInt)(n)` 截断，>4GB section CRC 只覆盖低 32 位 | **已修**:`ckpt_crc32` 分块 |
| 低 | rd_cc/rd_paf_arr/stream plen 先按文件长度字段分配后校验 | **已修**：全部按 fstat 文件大小做 sanity bound |
| 低 | `load_All_reads_fp` 指针数组 malloc 未初始化，部分加载无法安全释放 | **已修**:N_site/read_sperate/rsc 改 calloc;hom/het_cov 末尾才提交 |
| 低 | stream_append 持锁 fsync | **已修**:fsync 移出 mutex（覆盖语义不变） |
| 低 | stream_rescan realloc 丢旧指针 | **已修** |
| 低 | ha_ckpt_purge 清不掉不同 -r 的遗留 | **已修**：固定清扫 r0..r63 |
| 低 | "stale ignored" 死代码、HA_CKPT_NO_FLUSH 未用 | **已删** |
| 低 | params_hash 缺 rl_cut/sc_cut | **已补** |
| — | ckpt 默认开启的开销（报 173% 是载机上的误测） | 9036 安静机实测 **~1-3%**，维持默认开；--ckpt-int 0 可关 |
| — | ha_assemble_pair 回放陈旧 ecstr | 已确认安全：`ha_ckpt_on` 排除 sec_in 路径 |
| — | params_hash 字段运行中是否被改 | max_n_chain 单独快照处理；flag 的 HA_F_PARTITION 变更在建图阶段（ckpt 写入之后），无影响 |

验收：harness 重跑通过（SIGKILL→resume→4×GFA bit-exact,WAL 回放 323/2924)。
