#!/bin/bash
# 多进程干扰测试 —— 同 prefix 并发/交叉运行, 验证不串数据 + 可恢复
#   [MP1] 同 prefix 同输入, 双进程并发          -> 共享 WAL/快照, 看行为与输出
#   [MP2] A 被 kill -> B 同 prefix 跑完(purge) -> A resume: 进度被清但不得串数据
#   [MP3] 同 prefix 不同输入, 双进程并发        -> params 互斥, 输出不得混合
#   [MP4] 不同 prefix 并发 (对照)              -> 必须互不影响
#
# 用法: ./ckpt_mp_test.sh <hifiasm> <reads.fa> <reads_sub.fa> <工作目录> [线程] [阈值]
set -u
HIFIASM=${1:-../hifiasm}
READS=${2:-reads.fa}
READS2=${3:-reads_sub.fa}
DIR=${4:-.}
T=${5:-4}
THRESH=${6:-200000}
cd "$DIR" || exit 1

fsz() { stat -f%z "$1" 2>/dev/null || stat -c%s "$1" 2>/dev/null; }
md() { md5 -q "$1" 2>/dev/null || md5sum "$1" | cut -d' ' -f1; }
pick_gfa() { ls ${1}.*.gfa 2>/dev/null; }
wait_file() { local t=0; while [ $t -lt 3600 ]; do kill -0 $3 2>/dev/null || return 2; [ -f "$1" ] && [ "$(fsz "$1")" -ge "$2" ] && return 0; sleep 0.5; t=$((t+1)); done; return 1; }
fail=0
note() { echo "     $*"; }
# 只清 hifiasm 产物, 不能 rm -f P.* —— 会误删脚本自身输出 P.out
wipe() { rm -f $1.*.gfa $1.*.bed $1.*.bin $1.*.snap $1.*.ecstr $1.*.ovstr $1.*.tmp $1.*.lowQ.bed; }

cmp_set() { # $1=gfa前缀 $2=基准前缀 $3=标签 -> 全同=0 否则=1
    local ok=1 n=0 g f
    for g in $(pick_gfa "$1"); do f=${g#$1.}; [ "$(md $g)" = "$(md $2.$f 2>/dev/null)" ] || ok=0; n=$((n+1)); done
    [ $n -eq 0 ] && ok=0
    echo "  $3: $([ $ok -eq 1 ] && echo '同基准' || echo '不同/缺失') ($n 文件)"
    return $ok
}

echo "===== [0] 基线 bb / bs(子集输入) ====="
wipe bb; wipe bs
$HIFIASM -o bb -t$T "$READS" 2>bb.log || { echo "基线失败"; exit 1; }
$HIFIASM -o bs -t$T "$READS2" 2>bs.log || { echo "子集基线失败"; exit 1; }
NBB=$(pick_gfa bb | wc -l | tr -d ' '); echo "[0] bb GFA ${NBB}个, bs GFA $(pick_gfa bs | wc -l | tr -d ' ')个"

echo "===== [MP1] 同 prefix 同输入双进程并发 ====="
wipe mp
$HIFIASM -o mp -t$T "$READS" 2>mp_a.log & PA=$!
$HIFIASM -o mp -t$T "$READS" 2>mp_b.log & PB=$!
wait $PA; RCA=$?; wait $PB; RCB=$?
echo "  进程A rc=$RCA, 进程B rc=$RCB"
note "A日志尾部: $(tail -2 mp_a.log | tr '\n' '|')"
note "B日志尾部: $(tail -2 mp_b.log | tr '\n' '|')"
if [ $RCA -eq 0 ] && [ $RCB -eq 0 ]; then
    if cmp_set mp bb "MP1 输出 vs bb"; then echo "  MP1: 双进程输出与基线全等"; else echo "FAIL MP1: 输出与基线不等"; fail=1; fi
else
    echo "  MP1: 有进程非零退出(记录行为, 并发同写本就未定义)"
fi
ls mp.*.snap mp.*.ecstr mp.*.ovstr mp.*.tmp 2>/dev/null | while read f; do note "残留: $f ($(fsz $f)B)"; done

echo "===== [MP2] A kill -> B 同prefix跑完 -> A resume ====="
wipe m2
$HIFIASM -o m2 -t$T "$READS" 2>m2_a.log & PA=$!
wait_file m2.r1.ecstr "$THRESH" $PA || { echo "m2 没等到 ecstr"; kill -9 $PA 2>/dev/null; }
kill -9 $PA 2>/dev/null; wait $PA 2>/dev/null; sleep 1
note "A 被杀时现场: $(ls m2.* 2>/dev/null | tr '\n' ' ')"
$HIFIASM -o m2 -t$T "$READS" 2>m2_b.log   # B 跑完(会 purge 掉 A 的 ckpt)
note "B 完成后残留 ckpt: $(ls m2.*.snap m2.*.ecstr m2.*.ovstr 2>/dev/null | tr '\n' ' ')"
$HIFIASM -o m2 -t$T "$READS" 2>m2_c.log && RC=0 || RC=1
N=$(grep -c "resumed" m2_c.log)
echo "  A resume: rc=$RC, resumed行数=$N ($(grep -o 'resumed [0-9]*/[0-9]*' m2_c.log | tail -1))"
cmp_set m2 bb "MP2 输出 vs bb"
for g in $(pick_gfa m2); do f=${g#m2.}; [ "$(md $g)" = "$(md bb.$f)" ] || fail=1; done

echo "===== [MP3] 同 prefix 不同输入双进程并发 ====="
wipe m3
$HIFIASM -o m3 -t$T "$READS" 2>m3_a.log & PA=$!
$HIFIASM -o m3 -t$T "$READS2" 2>m3_b.log & PB=$!
wait $PA; RCA=$?; wait $PB; RCB=$?
echo "  A(全集) rc=$RCA, B(子集) rc=$RCB"
note "A resumed/reject: $(grep -cE 'resumed|incompatible' m3_a.log)"
note "B resumed/reject: $(grep -cE 'resumed|incompatible' m3_b.log)"
if [ $RCA -eq 0 ] && [ $RCB -eq 0 ]; then
    if cmp_set m3 bb >/dev/null; then W=bb;
    elif cmp_set m3 bs >/dev/null; then W=bs; else W=NONE; fi
    echo "  最终 GFA 归属: $W (必须等于某个完整基准, 不得是混合态)"
    [ "$W" = "NONE" ] && { echo "FAIL MP3: 混合态输出"; fail=1; } || echo "  MP3: 输出未串数据(归属 $W)"
else
    echo "  MP3: 有进程非零退出(记录行为)"; fi
ls m3.*.snap m3.*.ecstr 2>/dev/null | while read f; do note "残留: $f ($(fsz $f)B)"; done

echo "===== [MP4] 不同 prefix 并发(对照) ====="
wipe d1; wipe d2
$HIFIASM -o d1 -t$T "$READS" 2>d1.log & PA=$!
$HIFIASM -o d2 -t$T "$READS2" 2>d2.log & PB=$!
wait $PA; RCA=$?; wait $PB; RCB=$?
ok=1
for g in $(pick_gfa d1); do f=${g#d1.}; [ "$(md $g)" = "$(md bb.$f)" ] || ok=0; done
for g in $(pick_gfa d2); do f=${g#d2.}; [ "$(md $g)" = "$(md bs.$f)" ] || ok=0; done
[ $RCA -eq 0 ] && [ $RCB -eq 0 ] && [ $ok -eq 1 ] && echo "  MP4: 完全隔离, 双边 bit-exact" || { echo "FAIL MP4"; fail=1; }

echo "===== 汇总 ====="
[ $fail -eq 0 ] && echo "  全部通过" || echo "  !! 存在失败项"
exit $fail
