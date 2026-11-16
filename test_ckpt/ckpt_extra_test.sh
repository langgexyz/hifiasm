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

if want G; then
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

fi
if want H; then
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

fi
if want I; then
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

fi
if want J; then
echo "===== [J] resume 换线程数 ($T -> 2) ====="
rm -f tc.*
if kill_run tc "tc.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    $HIFIASM -o tc -t2 "$READS" 2>tc_r.log
    grep -q "resumed" tc_r.log && echo "OK   J: -t2 正常 resume" || { echo "FAIL J: 未 resume"; fail=1; }
    cmp_gfa bb tc J
else echo "SKIP J"; fi

fi
if want L; then
echo "===== [L] resume 加 -i (忽略全部 ckpt) ====="
rm -f li.*
if kill_run li "li.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    $HIFIASM -o li -t$T -i "$READS" 2>li_r.log
    grep -q "resumed" li_r.log && { echo "FAIL L: -i 仍 resume"; fail=1; } || echo "OK   L: -i 忽略 ckpt"
    cmp_gfa bb li L
else echo "SKIP L"; fi

fi
if want M; then
echo "===== [M] resume 加 --ckpt-int 0 (关闭 ckpt) ====="
rm -f nc.*
if kill_run nc "nc.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    $HIFIASM -o nc -t$T --ckpt-int 0 "$READS" 2>nc_r.log
    grep -q "resumed" nc_r.log && { echo "FAIL M: ckpt-int=0 仍 resume"; fail=1; } || echo "OK   M: ckpt 关闭"
    cmp_gfa bb nc M
else echo "SKIP M"; fi

fi
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
    # 两条合法恢复路径: ckpt(r2.snap+ovstr) 或 .ec/.ovlp 缓存(kill 前已写完)
    if grep -q "resumed checkpoint of round 2" gr_r.log; then echo "OK   N: 从 r2.snap 恢复"
    elif grep -q "loaded corrected reads" gr_r.log; then echo "OK   N: 从 ovlp 缓存恢复"
    else echo "FAIL N: 未按预期 resume"; fail=1; fi
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
cp "$READS" "reads_ad.orig"
rm -f dp.*
$HIFIASM -o dp -t$T "$READS" 2>dp_a.log & P=$!
wait_file "dp.r0.ecstr" 1024 $P && kill -9 $P
sleep 1
command -v python3 >/dev/null || { echo "SKIP AD(无 python3)"; cp "reads_ad.orig" "$READS"; }
if command -v python3 >/dev/null; then
python3 - "$READS" <<'PY'
import sys
d = bytearray(open(sys.argv[1],"rb").read())
i = 1000000
while d[i:i+1] not in (b"A",b"C",b"G",b"T"): i += 1
d[i:i+1] = b"C" if d[i:i+1] != b"C" else b"A"
open(sys.argv[1],"wb").write(d)
PY
touch -r "reads_ad.orig" "$READS" 2>/dev/null
$HIFIASM -o dp -t$T "$READS" 2>dp_b.log
if grep -qE "resumed" dp_b.log; then
    echo "FAIL AD(深变异仍被 resume)"; fail=1
else
    echo "OK(深变异被拒)"
    cmp_gfa bb dp AD
fi
fi
cp "reads_ad.orig" "$READS"
fi

if want AE; then
echo "===== [AE] WAL rid字段损坏 -> CRC拒绝(不得错槽位应用) ====="
# rec = MAGIC(2)|stage(4)|rid(8)|plen(8)|payload|crc(4); 翻第二条记录的 rid 低位
rm -f rf.*
if kill_run rf "rf.r1.ecstr" "$THRESH" "$HIFIASM -t$T $READS"; then
    python3 - rf.r1.ecstr <<'PY' 2>/dev/null
import sys, struct
f = sys.argv[1]; d = bytearray(open(f,'rb').read())
off = 40                       # 跳过文件头
plen = struct.unpack('<Q', d[off+14:off+22])[0]
off += 26 + plen               # 第2条记录起点
d[off+6] ^= 0x01               # rid 最低位字节翻转 (rid±1, 仍在界内)
open(f,'wb').write(d)
print("flipped rid byte at", off+6)
PY
    $HIFIASM -o rf -t$T $READS 2>rf_r.log
    grep -E "resumed" rf_r.log | head -2
    cmp_gfa bb rf AE
else echo "SKIP AE"; fi
fi

if want AF; then
echo "===== [AF] .ec/.ovlp 缓存身份: 匹配加载/参数变拒绝/stamp缺拒绝 ====="
# a) 同参重跑 -> 应加载缓存(快路径)
rm -f af.*
$HIFIASM -o af -t$T $READS 2>af_a.log
if grep -q "loaded corrected reads" af_a.log; then echo "FAIL AF-a(首跑就加载?)"; fail=1; fi
$HIFIASM -o af -t$T $READS 2>af_b.log
if grep -q "loaded corrected reads" af_b.log; then
    echo "OK AF-a(同参缓存命中)"
    cmp_gfa bb af AF-a
else echo "FAIL AF-a(同参缓存未命中)"; fail=1; fi
# b) 改 --hg-size 重跑 -> 必须拒绝缓存, 重新纠错
$HIFIASM -o af -t$T --hg-size 100m $READS 2>af_c.log
if grep -q "loaded corrected reads" af_c.log; then
    echo "FAIL AF-b(参数变了还加载旧缓存)"; fail=1
elif grep -qE "ec round" af_c.log; then
    echo "OK AF-b(参数变->拒缓存重算)"
else echo "FAIL AF-b(既没加载也没重算?)"; fail=1; fi
# c) 删 stamp 重跑 -> 必须拒绝(不可验证来源)
rm -f af.ovlp.params
$HIFIASM -o af -t$T $READS 2>af_d.log
if grep -q "loaded corrected reads" af_d.log; then
    echo "FAIL AF-c(无stamp仍加载)"; fail=1
else echo "OK AF-c(无stamp->拒)"; fi
fi

if want AG; then
echo "===== [AG] --dbg-gfa 重载 debug graph 不被 ovlp.params stamp 误伤 ====="
rm -f dg2.*
$HIFIASM -o dg2 -t$T --dbg-gfa $READS 2>dg2_a.log
$HIFIASM -o dg2 -t$T --dbg-gfa $READS 2>dg2_b.log
if grep -q "debug gfa has been loaded" dg2_b.log; then
    echo "OK AG(dbg-gfa 重载正常)"
else
    echo "FAIL AG(dbg-gfa 重载被挡)"; fail=1
fi
fi

if want AH; then
echo "===== [AH] 多输入: resume 加第二个文件 -> 必须拒绝 ====="
head -2000 "$READS" > mi_a.fa; tail -1200 "$READS" > mi_b.fa
rm -f mi.*
$HIFIASM -o mi -t$T mi_a.fa 2>mi_a.log & P=$!
wait_file "mi.r0.ecstr" 1024 $P && kill -9 $P; sleep 1
$HIFIASM -o mi -t$T mi_a.fa mi_b.fa 2>mi_b.log
if grep -qE "resumed" mi_b.log; then
    echo "FAIL AH(输入集变了还 resume)"; fail=1
else echo "OK AH(输入集变更被拒)"; fi
rm -f mi_a.fa mi_b.fa
fi

if want AI; then
echo "===== [AI] .ovlp.source.bin 中段损坏 -> CRC拒+重算 ====="
rm -f ai.*
$HIFIASM -o ai -t$T "$READS" 2>ai_a.log
SZ=$(fsz ai.ovlp.source.bin); OFF=$((SZ/2))
printf '\xFF' | dd of=ai.ovlp.source.bin bs=1 seek=$OFF conv=notrunc 2>/dev/null
$HIFIASM -o ai -t$T "$READS" 2>ai_b.log
if grep -qE 'loaded corrected reads' ai_b.log; then
    echo "FAIL AI(损坏的 ovlp.source 被接受)"; fail=1
elif grep -qE 'corrupt|bad crc' ai_b.log; then
    echo "OK AI(损坏被拒并重算)"
else echo "FAIL AI(未见预期拒绝)"; fail=1; fi
cmp_gfa bb ai AI
fi

if want AJ; then
echo "===== [AJ] .ec.bin 中段损坏 -> footer CRC拒+重算 ====="
rm -f aj.*
$HIFIASM -o aj -t$T "$READS" 2>aj_a.log
SZ=$(fsz aj.ec.bin); OFF=$((SZ*2/3))
printf '\xAA' | dd of=aj.ec.bin bs=1 seek=$OFF conv=notrunc 2>/dev/null
$HIFIASM -o aj -t$T "$READS" 2>aj_b.log
if grep -qE 'loaded corrected reads' aj_b.log; then
    echo "FAIL AJ(损坏的 .ec.bin 被接受)"; fail=1
elif grep -qE 'bad crc|corrupt' aj_b.log; then
    echo "OK AJ(.ec.bin 损坏被拒并重算)"
else echo "FAIL AJ(未见预期拒绝)"; fail=1; fi
cmp_gfa bb aj AJ
fi

if want AK; then
echo "===== [AK] 部分缓存: 只有 .ec.bin -> R_INF 归零重算 ====="
rm -f ak.*
$HIFIASM -o ak -t$T "$READS" 2>ak_a.log
rm -f ak.ovlp.source.bin ak.ovlp.reverse.bin   # 只留 .ec + params
$HIFIASM -o ak -t$T "$READS" 2>ak_b.log
if grep -qE 'loaded corrected reads' ak_b.log; then
    echo "FAIL AK(半个缓存被当成完整缓存)"; fail=1
elif grep -qE 'Reads has been loaded' ak_b.log; then
    echo "OK AK(.ec 加载后 ovlp 缺失, 归零重算)"
else echo "FAIL AK(未见预期路径)"; fail=1; fi
cmp_gfa bb ak AK
fi

if want AL; then
echo "===== [AL] .ec.bin 截尾 -> 提前 fread 失败 -> 归零重算 ====="
rm -f al.*
$HIFIASM -o al -t$T "$READS" 2>al_a.log
SZ=$(fsz al.ec.bin); TR=$((SZ*3/4))
if ! truncate -s $TR al.ec.bin 2>/dev/null; then
    dd if=al.ec.bin of=al.ec.bin.tr bs=$TR count=1 2>/dev/null && mv al.ec.bin.tr al.ec.bin
fi
$HIFIASM -o al -t$T "$READS" 2>al_b.log
if grep -qE 'loaded corrected reads' al_b.log; then
    echo "FAIL AL(截断的 .ec.bin 被接受)"; fail=1
else echo "OK AL(截断被拒并重算)"; fi
cmp_gfa bb al AL
fi

if want AM; then
echo "===== [AM] 部分缓存: 只有 .ovlp.* -> 无 .ec 重算 ====="
rm -f am.*
$HIFIASM -o am -t$T "$READS" 2>am_a.log
rm -f am.ec.bin
$HIFIASM -o am -t$T "$READS" 2>am_b.log
if grep -qE 'loaded corrected reads' am_b.log; then
    echo "FAIL AM(无 .ec 却走了缓存)"; fail=1
else echo "OK AM(缺 .ec 正常重算)"; fi
cmp_gfa bb am AM
fi

if want AN; then
echo "===== [AN] 并发写缓存: 同 prefix 双进程 -> 第三跑必须 bit-exact ====="
rm -f an.*
# 两个进程同 prefix 同输入并发, .tmp 交错写, 后 rename 者赢
$HIFIASM -o an -t2 "$READS" 2>an_a.log & PA=$!
$HIFIASM -o an -t2 "$READS" 2>an_b.log & PB=$!
wait $PA; wait $PB
# 第三跑: 缓存要么有效要么被 CRC 拒, 输出都必须对
$HIFIASM -o an -t$T "$READS" 2>an_c.log
if grep -qE 'loaded corrected reads' an_c.log; then
    echo "OK AN(并发后缓存完好命中)"
elif grep -qE 'corrupt|bad crc|WARN' an_c.log; then
    echo "OK AN(交错缓存被 CRC 拒并重算)"
else
    echo "FAIL AN(第三跑既未命中也未拒绝)"; fail=1
fi
cmp_gfa bb an AN
fi

if want AO; then
echo "===== [AO] 缓存写入中 kill -> .tmp 忽略 + params 缺席 -> 重算 ====="
rm -f ao.*
# AO-a: 在 "Writing reads to disk" 出现瞬间杀 (.ec.bin.tmp 写一半)
$HIFIASM -o ao -t$T "$READS" 2>ao_a.log & P=$!
t=0
while [ $t -lt 3600 ]; do
    kill -0 $P 2>/dev/null || break
    grep -q 'Writing reads to disk' ao_a.log 2>/dev/null && break
    sleep 0.1; t=$((t+1))
done
kill -9 $P 2>/dev/null; wait $P 2>/dev/null; sleep 1
TMPF=$(ls ao.*.tmp 2>/dev/null | head -1)
echo "[AO-a] kill @ $(ls ao.*.tmp 2>/dev/null | tr '\n' ' ')"
$HIFIASM -o ao -t$T "$READS" 2>ao_b.log
if [ -n "$TMPF" ] && grep -qE 'loaded corrected reads' ao_b.log; then
    echo "FAIL AO-a(半写缓存被接受)"; fail=1
elif [ -z "$TMPF" ]; then
    echo "OK AO-a(kill 未落在写窗口, 跑完即合法缓存)"
else echo "OK AO-a(无 params 拒缓存重算)"; fi
cmp_gfa bb ao AO-a
# AO-b: 在 "Writing ma_hit_ts" 出现瞬间杀 (.ec 已提交, ovlp.tmp 写一半)
rm -f ao.*
$HIFIASM -o ao -t$T "$READS" 2>ao_c.log & P=$!
t=0
while [ $t -lt 3600 ]; do
    kill -0 $P 2>/dev/null || break
    grep -q 'ma_hit_ts to disk' ao_c.log 2>/dev/null && break
    sleep 0.1; t=$((t+1))
done
kill -9 $P 2>/dev/null; wait $P 2>/dev/null; sleep 1
TMPF=$(ls ao.*.tmp 2>/dev/null | head -1)
echo "[AO-b] kill @ $(ls ao.*.tmp ao.ec.bin 2>/dev/null | tr '\n' ' ')"
$HIFIASM -o ao -t$T "$READS" 2>ao_d.log
if [ -n "$TMPF" ] && grep -qE 'loaded corrected reads' ao_d.log; then
    echo "FAIL AO-b(半写状态被当成完整缓存)"; fail=1
elif [ -z "$TMPF" ]; then
    echo "OK AO-b(kill 未落在写窗口)"
else echo "OK AO-b(.ec 已提交但无 params -> 归零重算)"; fi
cmp_gfa bb ao AO-b
fi

if want AP; then
echo "===== [AP] header-only WAL: 续写不重复 header、二次崩溃后记录存活 ====="
rm -f ap.*
$HIFIASM -o ap -t$T "$READS" >ap_a.log 2>&1 & P=$!
wait_file "ap.r1.ecstr" 1024 $P && kill -9 $P; wait $P 2>/dev/null; sleep 1
if [ -f ap.r1.ecstr ]; then
    # 手术式截断成"只有 40B header、0 条记录"的 WAL(模拟写完 header 即断电)
    truncate -s 40 ap.r1.ecstr
    # run2: 在 header-only WAL 上续写, 不得追加第二个 header
    $HIFIASM -o ap -t$T "$READS" >ap_b.log 2>&1 & P=$!
    wait_file "ap.r1.ecstr" 100000 $P; sleep 2; kill -9 $P 2>/dev/null; wait $P 2>/dev/null; sleep 1
    HDRS=$(LC_ALL=C grep -ao 'HICKPT01' ap.r1.ecstr 2>/dev/null | wc -l | tr -d ' ')
    # run3: 再 resume —— run2 续写的记录应被回放(旧 bug 会因其前多出的第二
    # 个 header 在 rescan 时把全部记录判为坏尾而丢弃)
    $HIFIASM -o ap -t$T "$READS" >ap_c.log 2>&1
    R3=$(sed -n 's/.*resumed \([0-9]*\)\/[0-9]* reads.*/\1/p' ap_c.log | head -1)
    HIT=$(grep -c 'loaded corrected reads and overlaps' ap_c.log)
    cmp_gfa bb ap AP
    if [ "$HDRS" = "1" ] && { [ "${R3:-0}" -gt 0 ] || [ "$HIT" -gt 0 ]; }; then
        echo "OK AP(header-only WAL 续写无重复 header, 二次崩溃记录存活/缓存命中)"
    else echo "FAIL AP hdrs=$HDRS r3_resumed=${R3:-NA} hit=$HIT"; fail=1; fi
else echo "SKIP AP(未产生 r1.ecstr)"; kill -9 $P 2>/dev/null; wait $P 2>/dev/null; fi
fi

if want AQ; then
echo "===== [AQ] --hom-cov 变更 -> 快照/WAL/缓存一律拒绝重算 ====="
rm -f aq.* aq40.*
$HIFIASM -o aq -t$T --hom-cov 20 "$READS" >aq_a.log 2>&1 & P=$!
wait_file "aq.r1.ecstr" 1024 $P && kill -9 $P; wait $P 2>/dev/null; sleep 1
$HIFIASM -o aq -t$T --hom-cov 40 "$READS" >aq_b.log 2>&1
if grep -q 'resumed checkpoint\|loaded corrected reads and overlaps' aq_b.log; then
    echo "FAIL AQ(--hom-cov 变更后旧状态被复用)"; fail=1
else echo "OK AQ(--hom-cov 变更 -> 拒绝旧状态, 干净重算)"; fi
# 输出正确性:与 --hom-cov 40 全量基线 bit-exact
$HIFIASM -o aq40 -t$T --hom-cov 40 "$READS" >aq40.log 2>&1
cmp_gfa aq40 aq AQ
fi

# ---- 反例组:构造"校验通过但语义恶意/陈旧混合"的脏状态,验证整体拒绝 ----
# 共用:取真实 WAL 的前 40B header(自带匹配 params_hash),追加伪造记录
if want AR; then
echo "===== [AR] 伪造记录: CRC合法但 scc.n=2^63 (2*n 回绕绕过边界) ====="
if command -v python3 >/dev/null; then
rm -f fr.*
$HIFIASM -o fr -t$T "$READS" >fr_a.log 2>&1 & P=$!
wait_file "fr.r1.ecstr" 1024 $P && kill -9 $P; wait $P 2>/dev/null; sleep 1
if [ -f fr.r1.ecstr ]; then
python3 - <<'PYEOF'
import zlib, struct
d = open('fr.r1.ecstr','rb').read()
hdr = d[:40]                                   # 保留真实 header(参数指纹匹配)
# payload: trio(1)|het(4)|scc.n=2^63(8)|scc.f(1)  -> 2*n 回绕为 0
payload = struct.pack('<B',0) + struct.pack('<I',0) + struct.pack('<Q',1<<63) + struct.pack('<B',0)
rh = struct.pack('<H',0x5243) + struct.pack('<I',0) + struct.pack('<Q',0) + struct.pack('<Q',len(payload))
crc = zlib.crc32(rh[2:]); crc = zlib.crc32(payload, crc)
open('fr.r1.ecstr','wb').write(hdr + rh + payload + struct.pack('<I', crc & 0xffffffff))
PYEOF
$HIFIASM -o fr -t$T "$READS" >fr_b.log 2>&1
RC=$?
cmp_gfa bb fr AR
if [ $RC -eq 0 ]; then echo "OK AR(伪造溢出记录被拒, 无内存投毒, 干净重算)"; else echo "FAIL AR rc=$RC"; fail=1; fi
else echo "SKIP AR(未产生 r1.ecstr)"; fi
else echo "SKIP AR(无 python3)"; fi
fi

if want AT; then
echo "===== [AT] 伪造记录: rid 越界(合法 CRC) ====="
if command -v python3 >/dev/null; then
rm -f fo.*
$HIFIASM -o fo -t$T "$READS" >fo_a.log 2>&1 & P=$!
wait_file "fo.r1.ecstr" 1024 $P && kill -9 $P; wait $P 2>/dev/null; sleep 1
if [ -f fo.r1.ecstr ]; then
python3 - <<'PYEOF'
import zlib, struct
d = open('fo.r1.ecstr','rb').read()
hdr = d[:40]
payload = struct.pack('<B',0) + struct.pack('<I',0) + struct.pack('<Q',0) + struct.pack('<B',0)
rh = struct.pack('<H',0x5243) + struct.pack('<I',0) + struct.pack('<Q',99999999) + struct.pack('<Q',len(payload))
crc = zlib.crc32(rh[2:]); crc = zlib.crc32(payload, crc)
open('fo.r1.ecstr','wb').write(hdr + rh + payload + struct.pack('<I', crc & 0xffffffff))
PYEOF
$HIFIASM -o fo -t$T "$READS" >fo_b.log 2>&1
RC=$?
cmp_gfa bb fo AT
if [ $RC -eq 0 ]; then echo "OK AT(越界 rid 记录在 apply 前被拒)"; else echo "FAIL AT rc=$RC"; fail=1; fi
else echo "SKIP AT(未产生 r1.ecstr)"; fi
else echo "SKIP AT(无 python3)"; fi
fi

if want AS; then
echo "===== [AS] 跨run陈旧WAL移植: 死run的r1记录注入新run -> 回放仍一致 ====="
rm -f mg.* wal_A.bin
# run A: r1 中途杀 -> snap r0 + wal r1(A代), 保存 WAL
$HIFIASM -o mg -t$T "$READS" >mg_a.log 2>&1 & P=$!
wait_file "mg.r1.ecstr" 2048 $P && kill -9 $P; wait $P 2>/dev/null; sleep 1
if [ -f mg.r1.ecstr ]; then
cp mg.r1.ecstr wal_A.bin
rm -f mg.*
# run B: 全新跑, r1 中途杀 -> snap r0 + wal r1(B代,部分)
$HIFIASM -o mg -t$T "$READS" >mg_b.log 2>&1 & P=$!
wait_file "mg.r1.ecstr" 512 $P && kill -9 $P; wait $P 2>/dev/null; sleep 1
# 移植: 用 A 代 WAL 覆盖 B 代 WAL(同输入同参数 -> header 有效)
cp wal_A.bin mg.r1.ecstr
$HIFIASM -o mg -t$T "$READS" >mg_c.log 2>&1
RC=$?
R=$(sed -n 's/.*resumed \([0-9]*\)\/[0-9]* reads.*/\1/p' mg_c.log | head -1)
cmp_gfa bb mg AS
if [ $RC -eq 0 ]; then echo "OK AS(跨run陈旧WAL回放 resumed=${R:-NA}, bit-exact)"; else echo "FAIL AS rc=$RC"; fail=1; fi
else echo "SKIP AS(未产生 r1.ecstr)"; kill -9 $P 2>/dev/null; wait $P 2>/dev/null; fi
fi

echo "===== 汇总 ====="
[ $fail -eq 0 ] && echo "  全部通过" || echo "  !! 存在失败项, 见上方 FAIL/DIFF"
exit $fail
