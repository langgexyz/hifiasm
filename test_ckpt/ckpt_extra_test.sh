#!/bin/bash
# hifiasm ckpt 补充覆盖测试 —— 主脚本(kill_resume_test.sh)之外的缺口:
#   [F] final-overlap(OV)阶段 kill + resume     -> 必须 bit-exact
#   [G] 快照字节损坏                             -> 拒绝加载 + 干净重跑 + 仍 bit-exact
#   [H] 双崩溃(kill->resume->kill->resume)幂等   -> 必须 bit-exact
#   [I] 参数变更 resume(--rl-cut)               -> 拒绝旧 ckpt, 干净完成
#   [J] resume 换线程数                          -> 正常 resume + bit-exact
#   [K] ONT 数据 kill + resume                  -> 必须 bit-exact
#   [L] resume 加 -i                            -> 忽略 ckpt, 干净完成 + bit-exact
#   [M] resume 加 --ckpt-int 0                  -> 不 resume, 干净完成 + bit-exact
#
# 用法: ./ckpt_extra_test.sh <hifiasm> <reads.fa> <ont_reads.fq> <工作目录> [线程] [ecstr阈值] [ovstr阈值]
set -u
HIFIASM=${1:-../hifiasm}
READS=${2:-reads.fa}
ONT_READS=${3:-reads_ont.fq}
DIR=${4:-.}
T=${5:-4}
THRESH=${6:-1000000}
OVTH=${7:-100000}
cd "$DIR" || exit 1

pick_gfa() { ls ${1}.*.gfa 2>/dev/null; }
fsz() { stat -f%z "$1" 2>/dev/null || stat -c%s "$1" 2>/dev/null; }
md() { md5 -q "$1" 2>/dev/null || md5sum "$1" | cut -d' ' -f1; }
now() { date +%s; }
fail=0

# 等文件长到阈值; 进程提前退出视为失败(数据集太小/跑太快)
wait_file() { # $1=file $2=minsize $3=pid $4=timeout_s
    local t=0
    while [ $t -lt ${4:-3600} ]; do
        kill -0 $3 2>/dev/null || return 2
        [ -f "$1" ] && [ "$(fsz "$1")" -ge "$2" ] && return 0
        sleep 0.5; t=$((t+1))
    done
    return 1
}
wait_grep() { # $1=log $2=pattern $3=timeout_s
    local t=0
    while [ $t -lt ${3:-600} ]; do
        grep -q "$2" "$1" 2>/dev/null && return 0
        sleep 0.5; t=$((t+1))
    done
    return 1
}
cmp_gfa() { # $1=ref前缀 $2=目标前缀 $3=标签
    local n=0 ok=1 g f a b
    for g in $(pick_gfa "$1"); do
        f=${g#$1.}
        [ -f "$2.$f" ] || { echo "DIFF $3: $2.$f 缺失"; ok=0; continue; }
        a=$(md "$1.$f"); b=$(md "$2.$f")
        if [ "$a" = "$b" ] && [ -n "$a" ]; then echo "OK   $3 $f"; else echo "DIFF $3 $f  $1=$a $2=$b"; ok=0; fi
        n=$((n+1))
    done
    [ $n -eq 0 ] && { echo "DIFF $3: $1 无任何 GFA"; ok=0; }
    [ $ok -eq 1 ] || fail=1
}
kill_run() { # $1=prefix $2=触发文件 $3=阈值 $4=完整命令(字符串)
    rm -f $1.log
    $4 -o $1 2>$1.log & local pid=$!
    wait_file "$2" "$3" $pid || { echo "!! $1: 进程提前退出或超时(数据集太小?), 跳过kill"; kill -9 $pid 2>/dev/null; wait $pid 2>/dev/null; return 1; }
    kill -9 $pid; wait $pid 2>/dev/null; sleep 1; return 0
}
kill_run_log() { # $1=prefix $2=log正则(命中即杀) $3=完整命令 —— 用于一轮<1s的快速阶段
    rm -f $1.log
    $3 -o $1 2>$1.log & local pid=$!
    local t=0
    while [ $t -lt 3600 ]; do
        kill -0 $pid 2>/dev/null || { echo "!! $1: 触发前进程已退出"; return 1; }
        grep -qE "$2" $1.log 2>/dev/null && break
        sleep 0.1; t=$((t+1))
    done
    kill -9 $pid; wait $pid 2>/dev/null; sleep 1
    echo "[$1] 日志驱动 kill @ '$(grep -E "$2" $1.log | tail -1 | cut -c1-80)'"
    return 0
}
SECS=${SECS:-ALL}
want() { [ "$SECS" = "ALL" ] || echo " $SECS " | grep -q " $1 "; }

echo "===== [0] 公共基线 bb (-t$T, 带 ckpt) ====="
rm -f bb.*
T0=$(now); $HIFIASM -o bb -t$T "$READS" 2>bb.log || { echo "基线跑失败"; exit 1; }
echo "[0] 基线耗时: $(( $(now) - T0 ))s, GFA: $(pick_gfa bb | wc -l | tr -d ' ') 个"

if want F; then
echo "===== [F] final-overlap 阶段 kill + resume (日志驱动, ~50%处杀) ====="
rm -f ov.*
if kill_run_log ov "final overlap.*[3-9]0%" "$HIFIASM -t$T $READS"; then
    echo "[F] kill 现场: ovstr=$(fsz ov.ff.ovstr 2>/dev/null)B, snap: $(ls ov.*.snap 2>/dev/null | tr '\n' ' ')"
    $HIFIASM -o ov -t$T "$READS" 2>ov_r.log
    grep -E "resumed" ov_r.log | head -3
    cmp_gfa bb ov F
else echo "SKIP F"; fi
fi

echo "===== [G] 快照字节损坏 -> 拒绝 + 干净重跑 ====="
rm -f cs.*
if kill_run cs "cs.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    SNAP=$(ls cs.*.snap 2>/dev/null | head -1)
    if [ -n "$SNAP" ]; then
        OFF=$(( $(fsz "$SNAP") / 2 ))
        printf '\xFF' | dd of="$SNAP" bs=1 seek=$OFF conv=notrunc 2>/dev/null
        echo "[G] 已在 $SNAP 偏移 $OFF 处翻 1 字节"
    else echo "[G] !! 无 snap 可损(数据集太小?)"; fi
    $HIFIASM -o cs -t$T "$READS" 2>cs_r.log
    if grep -q "resumed checkpoint" cs_r.log; then echo "FAIL G: 损坏快照仍被 resume"; fail=1
    elif grep -qi "crc mismatch\|failed" cs_r.log; then echo "OK   G: 损坏快照被拒, 干净重跑"
    else echo "??   G: 未 resume 也未见拒绝信息, 查 cs_r.log"; fi
    cmp_gfa bb cs G
else echo "SKIP G"; fi

echo "===== [H] 双崩溃幂等 kill->resume->kill->resume ====="
rm -f dc.*
if kill_run dc "dc.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    SZ1=$(fsz dc.r1.ecstr)
    $HIFIASM -o dc -t$T "$READS" 2>dc_r1.log & pid=$!
    # resume 回放 WAL 后 ecstr 重新长过原长 => 确实续算了新工作, 此时再杀
    t=0; grown=0
    while [ $t -lt 3600 ]; do
        kill -0 $pid 2>/dev/null || break
        [ -f dc.r1.ecstr ] && [ "$(fsz dc.r1.ecstr)" -gt "$SZ1" ] && { grown=1; break; }
        sleep 0.5; t=$((t+1))
    done
    if kill -0 $pid 2>/dev/null && [ $grown -eq 1 ]; then
        kill -9 $pid; wait $pid 2>/dev/null; sleep 1
        echo "[H] 第二次 kill 于 resume 续算中 (ecstr $(fsz dc.r1.ecstr)B > ${SZ1}B)"
        $HIFIASM -o dc -t$T "$READS" 2>dc_r2.log
        grep -E "resumed" dc_r2.log | head -3
    else
        wait $pid 2>/dev/null; echo "[H] resume 一次跑完(数据太小), 只验证单次恢复"
    fi
    cmp_gfa bb dc H
else echo "SKIP H"; fi

echo "===== [I] 参数变更 resume (--hg-size 1m->2m) ====="
# 注: 不能用 --rl-cut —— 非 ONT 输入下被强制改写为 -1, 参数实际不变.
# --hg-size 已进入 params_hash, 变更必须导致快照+WAL 被拒.
rm -f pm.*
if kill_run pm "pm.r1.ecstr" "$THRESH" "$HIFIASM -t$T --hg-size 1m $READS"; then
    $HIFIASM -o pm -t$T --hg-size 2m "$READS" 2>pm_r.log && RC=0 || RC=1
    N=$(grep -c "resumed" pm_r.log)
    if [ $RC -eq 0 ] && [ "$N" -eq 0 ]; then echo "OK   I: 参数失配->拒 resume,干净跑完"
    else echo "FAIL I: rc=$RC resumed行数=$N"; fail=1; fi
else echo "SKIP I"; fi

echo "===== [J] resume 换线程数 ($T -> 2) ====="
rm -f tc.*
if kill_run tc "tc.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    $HIFIASM -o tc -t2 "$READS" 2>tc_r.log
    grep -q "resumed" tc_r.log && echo "OK   J: -t2 正常 resume" || { echo "FAIL J: 未 resume"; fail=1; }
    cmp_gfa bb tc J
else echo "SKIP J"; fi

echo "===== [L] resume 加 -i (忽略全部 ckpt) ====="
rm -f li.*
if kill_run li "li.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    $HIFIASM -o li -t$T -i "$READS" 2>li_r.log
    grep -q "resumed" li_r.log && { echo "FAIL L: -i 仍 resume"; fail=1; } || echo "OK   L: -i 忽略 ckpt"
    cmp_gfa bb li L
else echo "SKIP L"; fi

echo "===== [M] resume 加 --ckpt-int 0 (关闭 ckpt) ====="
rm -f nc.*
if kill_run nc "nc.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    $HIFIASM -o nc -t$T --ckpt-int 0 "$READS" 2>nc_r.log
    grep -q "resumed" nc_r.log && { echo "FAIL M: ckpt-int=0 仍 resume"; fail=1; } || echo "OK   M: ckpt 关闭"
    cmp_gfa bb nc M
else echo "SKIP M"; fi

if want K; then
echo "===== [K] ONT 模式 kill + resume (日志驱动, round2 ~50%处杀) ====="
if [ -f "$ONT_READS" ]; then
    rm -f ont.* obb.*
    $HIFIASM -o obb -t$T --ont "$ONT_READS" 2>obb.log || { echo "ONT 基线失败"; fail=1; }
    if kill_run_log ont "round 2/3.*[3-9]0%" "$HIFIASM -t$T --ont $ONT_READS"; then
        $HIFIASM -o ont -t$T --ont "$ONT_READS" 2>ont_r.log
        grep -E "resumed" ont_r.log | head -3
        cmp_gfa obb ont K
    else echo "SKIP K(kill)"; fi
else echo "SKIP K: 无 $ONT_READS"; fi
fi

if want N; then
echo "===== [N] graph 阶段 kill + resume (r2.snap+ovstr 接力) ====="
rm -f gr.*
if kill_run_log gr "graph: (build initial graph|ug extend|write contig)" "$HIFIASM -t$T $READS"; then
    $HIFIASM -o gr -t$T "$READS" 2>gr_r.log
    grep -q "resumed checkpoint of round 2" gr_r.log && echo "OK   N: 从 r2.snap 恢复" || { echo "FAIL N: 未按预期 resume"; fail=1; }
    cmp_gfa bb gr N
else echo "SKIP N"; fi
fi

if want O; then
echo "===== [O] WAL 中段单字节损坏 -> 截断回放 ====="
rm -f wc.*
if kill_run wc "wc.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    OFF=$(( $(fsz wc.r1.ecstr) / 2 ))
    printf '\xFF' | dd of=wc.r1.ecstr bs=1 seek=$OFF conv=notrunc 2>/dev/null
    echo "[O] 已在 r1.ecstr 偏移 $OFF 处翻 1 字节"
    $HIFIASM -o wc -t$T "$READS" 2>wc_r.log
    N2=$(grep -o "resumed [0-9]*/[0-9]* reads" wc_r.log | tail -1)
    echo "[O] WAL 回放: ${N2:-无}(截断点之前有效)"
    cmp_gfa bb wc O
else echo "SKIP O"; fi
fi

if want P; then
echo "===== [P] 残留 .tmp 垃圾文件 -> 忽略 ====="
rm -f st.*
echo "GARBAGE_NOT_A_SNAPSHOT" > st.r0.snap.tmp
head -c 4096 /dev/urandom > st.r1.ecstr.tmp 2>/dev/null || dd if=/dev/urandom of=st.r1.ecstr.tmp bs=4096 count=1 2>/dev/null
if kill_run st "st.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    $HIFIASM -o st -t$T "$READS" 2>st_r.log
    grep -q "resumed" st_r.log && echo "OK   P: .tmp 不影响 resume" || { echo "FAIL P"; fail=1; }
    cmp_gfa bb st P
else echo "SKIP P"; fi
fi

if want Q; then
echo "===== [Q] fmt_ver 字节损坏 -> 拒绝 ====="
rm -f fv.*
if kill_run fv "fv.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    SNAP=$(ls fv.*.snap | head -1)
    printf '\xEE' | dd of="$SNAP" bs=1 seek=8 conv=notrunc 2>/dev/null   # fmt_ver @offset8
    $HIFIASM -o fv -t$T "$READS" 2>fv_r.log
    grep -q "resumed checkpoint" fv_r.log && { echo "FAIL Q: fmt_ver 损坏仍 resume"; fail=1; } || echo "OK   Q: fmt_ver 不符被拒"
    cmp_gfa bb fv Q
else echo "SKIP Q"; fi
fi

if want R; then
echo "===== [R] 同名输入文件被改 -> 拒绝 ====="
rm -f df.*
cp "$READS" reads_bak.fa
if kill_run df "df.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    printf 'A' | dd of="$READS" bs=1 seek=100 conv=notrunc 2>/dev/null
    $HIFIASM -o df -t$T "$READS" 2>df_r.log && RC=0 || RC=1
    mv reads_bak.fa "$READS"   # 立刻恢复, 后续场景还要用
    N=$(grep -c "resumed" df_r.log)
    [ $RC -eq 0 ] && [ "$N" -eq 0 ] && echo "OK   R: 输入变化->拒 resume" || { echo "FAIL R: rc=$RC resumed=$N"; fail=1; }
else mv reads_bak.fa "$READS"; echo "SKIP R"; fi
fi

if want S; then
echo "===== [S] 快照回退链: 最新坏+旧的好 -> 用旧的 ====="
rm -f fb.*
if kill_run fb "fb.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    cp fb.r0.snap fb_snap.bak
    if kill_run fb "fb.r2.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
        OFF=$(( $(fsz fb.r1.snap) / 2 ))
        printf '\xFF' | dd of=fb.r1.snap bs=1 seek=$OFF conv=notrunc 2>/dev/null
        cp fb_snap.bak fb.r0.snap
        $HIFIASM -o fb -t$T "$READS" 2>fb_r.log
        grep -q "resumed checkpoint of round 0" fb_r.log && echo "OK   S: 回退到 r0.snap" || { echo "FAIL S: 未回退 r0"; fail=1; }
        cmp_gfa bb fb S
    else echo "SKIP S(二段kill)"; fi
else echo "SKIP S"; fi
fi

if want U; then
echo "===== [U] reads+count 早期 kill -> 无 ckpt 干净重跑 ====="
rm -f er.*
if kill_run_log er "reads\+count" "$HIFIASM -t$T $READS"; then
    $HIFIASM -o er -t$T "$READS" 2>er_r.log
    cmp_gfa bb er U
else echo "SKIP U"; fi
fi

if want V; then
echo "===== [V] 成功完成后同 prefix 重跑 -> 幂等 ====="
# bb 已完成且 ckpt 已 purge; 同 prefix 再跑必须全新计算且不恢复任何残留
for g in $(pick_gfa bb); do echo "$(md "$g") $g"; done > bb.md5
$HIFIASM -o bb -t$T "$READS" 2>bb2.log
grep -q "resumed" bb2.log && { echo "FAIL V: 完成后同prefix仍 resume"; fail=1; } || echo "OK   V: 无残留可恢复"
vok=1
for g in $(pick_gfa bb); do
    f=$(basename "$g")
    a=$(grep " $f" bb.md5 | cut -d' ' -f1); b=$(md "$g")
    [ "$a" = "$b" ] && echo "OK   V $f" || { echo "DIFF V $f"; vok=0; }
done
[ $vok -eq 1 ] || fail=1
fi

if want W; then
echo "===== [W] 仅 WAL 无 snap 的 resume (round0 杀) ====="
rm -f wo.*
if kill_run wo "wo.r0.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    ls wo.*.snap 2>/dev/null && echo "[W] 注意: kill 时已有 snap"
    $HIFIASM -o wo -t$T "$READS" 2>wo_r.log
    N=$(grep -o "resumed [0-9]*/[0-9]* reads" wo_r.log | tail -1)
    echo "[W] 无 snap, 仅 WAL 回放: ${N:-无}"
    grep -q "resumed checkpoint" wo_r.log && echo "[W] 注: kill 偏晚, snap 已写 -> 走了 snap+WAL 路径"
    cmp_gfa bb wo W
else echo "SKIP W"; fi
fi

if want X; then
echo "===== [X] 收尾阶段 kill (write contig output) ====="
rm -f xl.*
if kill_run_log xl "graph: write contig output" "$HIFIASM -t$T $READS"; then
    $HIFIASM -o xl -t$T "$READS" 2>xl_r.log
    cmp_gfa bb xl X
else echo "SKIP X"; fi
fi

if want Y; then
echo "===== [Y] ovstr 中段单字节损坏 -> ov_apply 截断回放 ====="
rm -f yc.*
if kill_run yc "yc.ff.ovstr" "$OVTH" "$HIFIASM -t$T $READS"; then
    OFF=$(( $(fsz yc.ff.ovstr) / 2 ))
    printf '\xFF' | dd of=yc.ff.ovstr bs=1 seek=$OFF conv=notrunc 2>/dev/null
    $HIFIASM -o yc -t$T "$READS" 2>yc_r.log
    N=$(grep -o "resumed [0-9]*/[0-9]* reads" yc_r.log | tail -1)
    echo "[Y] ovstr 回放: ${N:-无}"
    cmp_gfa bb yc Y
else echo "SKIP Y"; fi
fi

if want Z1; then
echo "===== [Z1] --dbg-ovec 门控 -> 不写任何 ckpt ====="
rm -f dg.*
$HIFIASM -o dg -t$T --dbg-ovec "$READS" 2>dg.log && RC=0 || RC=1
NC=$(ls dg.*.snap dg.*.ecstr dg.*.ovstr 2>/dev/null | wc -l | tr -d ' ')
[ $RC -eq 0 ] && [ "$NC" -eq 0 ] && echo "OK   Z1: 门控生效, 0 个 ckpt 文件" || { echo "FAIL Z1: rc=$RC ckpt文件=$NC"; fail=1; }
fi

if want Z2; then
echo "===== [Z2] --write-ec + kill/resume 共存 ====="
rm -f we.*
if kill_run we "we.r1.ecstr" "$THRESH" "$HIFIASM -t$T --write-ec $READS"; then
    $HIFIASM -o we -t$T --write-ec "$READS" 2>we_r.log
    grep -q "resumed" we_r.log && echo "OK   Z2: --write-ec 下正常 resume" || { echo "FAIL Z2"; fail=1; }
    cmp_gfa bb we Z2
    ls we.ec.bin >/dev/null 2>&1 && echo "OK   Z2: .ec.bin 已产出" || echo "[Z2] 注: 无 we.ec.bin"
else echo "SKIP Z2"; fi
fi

if want Z3; then
echo "===== [Z3] --ckpt-int 1 (激进 fsync) kill/resume ====="
rm -f ci.*
if kill_run ci "ci.r1.ecstr" "$THRESH" "$HIFIASM -t$T --ckpt-int 1 $READS"; then
    $HIFIASM -o ci -t$T --ckpt-int 1 "$READS" 2>ci_r.log
    cmp_gfa bb ci Z3
else echo "SKIP Z3"; fi
fi

if want Z4; then
echo "===== [Z4] --dbg-gfa 门控 -> 不写任何 ckpt ====="
rm -f dg4.*
$HIFIASM -o dg4 -t$T --dbg-gfa "$READS" 2>dg4.log && RC=0 || RC=1
NC=$(ls dg4.*.snap dg4.*.ecstr dg4.*.ovstr 2>/dev/null | wc -l | tr -d ' ')
[ $RC -eq 0 ] && [ "$NC" -eq 0 ] && echo "OK   Z4: 门控生效, 0 个 ckpt 文件" || { echo "FAIL Z4: rc=$RC ckpt文件=$NC"; fail=1; }
fi

if want AA; then
echo "===== [AA] snap 在但 WAL 被删 -> snap resume + 当轮重算 ====="
rm -f aa.*
if kill_run aa "aa.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    rm -f aa.r1.ecstr
    $HIFIASM -o aa -t$T "$READS" 2>aa_r.log
    grep -q "resumed checkpoint of round 0" aa_r.log && echo "OK   AA: snap 恢复" || { echo "FAIL AA"; fail=1; }
    grep -o "resumed [0-9]*/[0-9]* reads" aa_r.log | tail -1 || echo "[AA] 当轮无 WAL, 全量重算(预期)"
    cmp_gfa bb aa AA
else echo "SKIP AA"; fi
fi

if want AB; then
echo "===== [AB] -o 带子目录路径 -> ckpt 落在子目录 ====="
mkdir -p subd; rm -f subd/* 2>/dev/null
if kill_run subd/sd "subd/sd.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    ls subd/*.snap subd/*.ecstr >/dev/null 2>&1 && echo "OK   AB: ckpt 随 prefix 进子目录" || { echo "FAIL AB"; fail=1; }
    $HIFIASM -o subd/sd -t$T "$READS" 2>sd_r.log
    for g in $(pick_gfa bb); do
        f=${g#bb.}
        [ -f "subd/sd.$f" ] && [ "$(md subd/sd.$f)" = "$(md $g)" ] && echo "OK   AB $f" || { echo "DIFF AB $f"; fail=1; }
    done
else echo "SKIP AB"; fi
fi

if want AC; then
echo "===== [AC] snap 写入瞬间 kill -> 半成品 .tmp 不得被当快照 ====="
# 大数据集下 snap 序列化要秒级; 等 round prog 到 100% 立刻杀, 落在写盘窗口内
rm -f sw.*
$HIFIASM -o sw -t$T "$READS" 2>sw.log & pid=$!
t=0; hit=0
while [ $t -lt 7200 ]; do
    kill -0 $pid 2>/dev/null || break
    grep -qE "ec round 1/[0-9]+.*100%" sw.log 2>/dev/null && { hit=1; break; }
    sleep 0.05; t=$((t+1))
done
if [ $hit -eq 1 ] && kill -0 $pid 2>/dev/null; then
    kill -9 $pid; wait $pid 2>/dev/null; sleep 1
    echo "[AC] kill 现场: $(ls sw.*.snap* sw.*.ecstr 2>/dev/null | tr '\n' ' ')"
    $HIFIASM -o sw -t$T "$READS" 2>sw_r.log
    grep -E "resumed" sw_r.log | head -3
    cmp_gfa bb sw AC
else
    kill -9 $pid 2>/dev/null; wait $pid 2>/dev/null
    echo "SKIP AC(没逮到写盘窗口, 需要大数据集)"
fi
fi

if want AD; then
echo "===== [AD] 深变异输入(>4KB, size/mtime 相同) -> 必须拒绝 ckpt ====="
# 旧指纹只哈希前4KB; 在 1MB 处改一个碱基并保持 size/mtime 不变, 验证全文件指纹
cp "$READS" "$WORK/reads_ad.orig"
rm -f dp.*
$HIFIASM -o dp -t$T "$READS" 2>dp_a.log & P=$!
wait_file "$WORK/dp.r0.ecstr" 1024 $P && kill -9 $P
sleep 1
command -v python3 >/dev/null || { echo "SKIP AD(无 python3)"; cp "$WORK/reads_ad.orig" "$READS"; }
if command -v python3 >/dev/null; then
python3 - "$READS" <<'PY'
import sys
d = bytearray(open(sys.argv[1],"rb").read())
i = 1000000
while d[i:i+1] not in (b"A",b"C",b"G",b"T"): i += 1
d[i:i+1] = b"C" if d[i:i+1] != b"C" else b"A"
open(sys.argv[1],"wb").write(d)
PY
touch -r "$WORK/reads_ad.orig" "$READS" 2>/dev/null
$HIFIASM -o dp -t$T "$READS" 2>dp_b.log
if grep -qE "resumed" dp_b.log; then
    echo "FAIL AD(深变异仍被 resume)"; fail=1
else
    echo "OK(深变异被拒)"
    cmp_gfa bb dp AD
fi
fi
cp "$WORK/reads_ad.orig" "$READS"
fi

echo "===== 汇总 ====="
[ $fail -eq 0 ] && echo "  全部通过" || echo "  !! 存在失败项, 见上方 FAIL/DIFF"
exit $fail
