#!/bin/bash
# hifiasm WAL-checkpoint 断电续跑验证脚本
# 用法: ./kill_resume_test.sh <hifiasm二进制> <reads.fa> <工作目录> [线程数] [ecstr阈值字节]
#
# 对比效果:
#   [A] 基线完整跑耗时          vs  [C] kill 后 resume 耗时
#   日志 "resumed N/M reads" 显示捡回的工作量; GFA md5 必须全部一致
set -u
HIFIASM=${1:-../hifiasm}
READS=${2:-reads.fa}
DIR=${3:-.}
T=${4:-4}
THRESH=${5:-1000000}   # r1.ecstr 长到该字节数才 kill; 大数据集调大(如 512M)
cd "$DIR" || exit 1

pick_gfa() { ls ${1}.*.gfa 2>/dev/null; }
wait_file() { # $1=file $2=minsize $3=timeout_s
    local t=0
    while [ $t -lt ${3:-300} ]; do
        [ -f "$1" ] && [ "$(fsz "$1")" -ge "${2:-1}" ] && return 0
        sleep 0.5; t=$((t+1))
    done
    return 1
}
fsz() { stat -f%z "$1" 2>/dev/null || stat -c%s "$1" 2>/dev/null; }
md() { md5 -q "$1" 2>/dev/null || md5sum "$1" | cut -d' ' -f1; }
now() { date +%s; }

echo "===== [A] 基线完整跑 (-t$T) ====="
rm -f base.*
T0=$(now); $HIFIASM -o base -t$T "$READS" 2>base.log; T1=$(now)
BASE_S=$((T1-T0)); echo "[A] 基线耗时: ${BASE_S}s"

echo "===== [B] kill 于纠错第2轮中途 (SIGKILL, 等价断电) ====="
rm -f kr.*
$HIFIASM -o kr -t$T "$READS" 2>kr.log & PID=$!
T0=$(now)
wait_file kr.r1.ecstr "$THRESH" 7200 || { echo "没等到 kr.r1.ecstr, 数据集太小?"; ls kr.*; kill -9 $PID 2>/dev/null; exit 1; }
kill -9 $PID; T1=$(now); LOST_S=$((T1-T0))
sleep 1
echo "[B] 崩溃前已耗时: ${LOST_S}s"
echo "[B] 崩溃现场: $(ls kr.*.snap kr.*.ecstr kr.*.ovstr 2>/dev/null | tr '\n' ' ')  [r1.ecstr=$(fsz kr.r1.ecstr)B]"

echo "===== [C] 原命令重跑 (自动 resume) ====="
T0=$(now); $HIFIASM -o kr -t$T "$READS" 2>kr_r.log; T1=$(now)
RESUME_S=$((T1-T0)); echo "[C] resume 耗时: ${RESUME_S}s"
grep -E "resumed|checkpoint" kr_r.log | head -8
R_N=$(grep -o "resumed [0-9]*/[0-9]* reads" kr_r.log | tail -1)
echo "[C] WAL 回放: ${R_N:-无}"

echo "===== [D] ckpt 开销: --ckpt-int 0 跑一遍 ====="
rm -f nock.*
T0=$(now); $HIFIASM -o nock -t$T --ckpt-int 0 "$READS" 2>nock.log; T1=$(now)
NOCK_S=$((T1-T0)); echo "[D] 无 checkpoint 耗时: ${NOCK_S}s (基线${BASE_S}s, 开销=$((BASE_S-NOCK_S))s ≈$(( (BASE_S-NOCK_S)*100/(NOCK_S+1) ))%)"

echo "===== [E] GFA 一致性对比 ====="
fail=0
for g in $(pick_gfa base); do
    f=${g#base.}
    a=$(md base.$f); b=$(md kr.$f 2>/dev/null); c=$(md nock.$f 2>/dev/null)
    if [ "$a" = "$b" ] && [ "$a" = "$c" ] && [ -n "$a" ]; then echo "OK   $f"; else echo "DIFF $f  base=$a kr=$b nock=$c"; fail=1; fi
done

echo "===== 汇总 ====="
echo "  基线(带ckpt): ${BASE_S}s | 崩溃前耗时: ${LOST_S}s | resume: ${RESUME_S}s | 无ckpt: ${NOCK_S}s"
echo "  断电损失上界 = resume 中重算的部分; WAL 捡回: ${R_N:-无}"
[ $fail -eq 0 ] && echo "  全部一致: 断点续跑结果与完整跑 bit-exact" || echo "  !! 存在差异"
exit $fail
