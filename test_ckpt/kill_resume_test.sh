#!/bin/bash
# hifiasm WAL-checkpoint 断电续跑验证脚本
# 用法: ./kill_resume_test.sh <hifiasm二进制> <reads.fa> <工作目录> [线程数]
# 产出: 基线 GFA 与 kill+resume GFA 的 md5 对比；resume 的省时从日志可读
set -u
HIFIASM=${1:-../hifiasm}
READS=${2:-reads.fa}
DIR=${3:-.}
T=${4:-4}
cd "$DIR" || exit 1

pick_gfa() { ls ${1}.bp.p_ctg.gfa ${1}.bp.r_utg.gfa ${1}.bp.hap1.p_ctg.gfa ${1}.bp.hap2.p_ctg.gfa 2>/dev/null; }

wait_file() { # $1=file $2=minsize $3=timeout_s
    local t=0
    while [ $t -lt ${3:-300} ]; do
        [ -f "$1" ] && [ "$(stat -f%z "$1" 2>/dev/null || stat -c%s "$1" 2>/dev/null)" -ge "${2:-1}" ] && return 0
        sleep 0.5; t=$((t+1))
    done
    return 1
}
fsz() { stat -f%z "$1" 2>/dev/null || stat -c%s "$1" 2>/dev/null; }
md() { md5 -q "$1" 2>/dev/null || md5sum "$1" | cut -d' ' -f1; }

echo "===== [1/4] 基线完整跑 (-t$T) ====="
rm -f base.*
{ time $HIFIASM -o base -t$T "$READS" ; } 2>base.log
tail -2 base.log

echo "===== [2/4] kill 于纠错第2轮中途 ====="
rm -f kr.*
$HIFIASM -o kr -t$T "$READS" 2>kr.log & PID=$!
# r1.ecstr = 第2轮纠错(r=1)的 per-read WAL; 等到它长到 ~20% 再杀
wait_file kr.r1.ecstr 1000000 3600 || { echo "没等到 r1.ecstr, 数据集太小?"; ls kr.*; kill -9 $PID 2>/dev/null; exit 1; }
kill -9 $PID
sleep 1
echo "崩溃现场: $(ls kr.*.snap kr.*.ecstr kr.*.ovstr 2>/dev/null | tr '\n' ' ')  [r1.ecstr=$(fsz kr.r1.ecstr)B]"

echo "===== [3/4] 原命令重跑 (自动 resume) ====="
{ time $HIFIASM -o kr -t$T "$READS" ; } 2>kr_r.log
grep -E "resumed|checkpoint" kr_r.log | head -8
tail -2 kr_r.log

echo "===== [4/4] GFA 一致性对比 ====="
fail=0
for g in $(pick_gfa base); do
    f=${g#base.}
    a=$(md base.$f); b=$(md kr.$f 2>/dev/null)
    if [ "$a" = "$b" ] && [ -n "$a" ]; then echo "OK   $f"; else echo "DIFF $f  base=$a kr=$b"; fail=1; fi
done
[ $fail -eq 0 ] && echo "全部一致: 断点续跑结果与完整跑 bit-exact" || echo "!! 存在差异"
exit $fail
