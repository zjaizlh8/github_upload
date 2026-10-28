#!/usr/bin/env bash
# tooth_x8.sh —— X8（tol 有数值 ⟺ 对象面）的齿：反例必须被咬住（2026-10-29 · §179.7）
#
# 三个篡改（都只改一行，其余 197 行原样）：
#   c1 过程面带 tol：hl:7（FLOW）tol 改成 3e+00 —— 这个值【正是输出侧真印出来的】⟹ X2 对不上不算错
#                    （旧审计 FAIL=0 实测；靠 X8 才有判别力）
#   c2 对象面缺 tol：bsd:5（FORCED）tol 改成 NA —— 旧审计已由 X6 误撞，但仍须被 X8b/X8a 直接咬住
#   c3 未归面桶   ：首行 class 改成 WEIRD —— 旧审计 FAIL=0（面桶漏格无人管）
# 通过判据：c1/c2/c3 各 FAIL>=1 且 RESULT 非 0；真台账 FAIL=0。
set -u
cd "$(dirname "$0")/.." || exit 1
LED=matrix/forcedness_ledger.csv
TMP=$(mktemp -d)
n_pass=0; n_fail=0
ok(){ n_pass=$((n_pass+1)); printf '  [PASS] %s\n' "$1"; }
ng(){ n_fail=$((n_fail+1)); printf '  [FAIL] %s\n' "$1"; }

awk -F, 'BEGIN{OFS=","} $1=="hl"&&$2=="7"{$5="3e+00"} {print}' "$LED" > "$TMP/c1.csv"
awk -F, 'BEGIN{OFS=","} $1=="bsd"&&$2=="5"{$5="NA"} {print}' "$LED" > "$TMP/c2.csv"
awk -F, 'BEGIN{OFS=","} NR==2{$3="WEIRD"} {print}' "$LED" > "$TMP/c3.csv"

echo "[齿] X8 · tol 有值 ⟺ 对象面（产物 audit/x8_face.awk）"
for c in c1 c2 c3; do
    x8out=$(awk -f audit/x8_face.awk "$TMP/$c.csv" 2>/dev/null)
    awk_x=$?
    lines=$(bash audit/verify_ledger.sh "$TMP/$c.csv" 2>/dev/null | grep -c '^  \[FAIL\]')
    if [ "$awk_x" -ne 0 ] && [ "$lines" -ge 1 ]
    then ok "$c 被咬住（x8_face exit=$awk_x · verify_ledger FAIL 行=$lines · $x8out）"
    else ng "$c 逃过（x8_face exit=$awk_x · verify_ledger FAIL 行=$lines · $x8out）"; fi
done

real_out=$(awk -f audit/x8_face.awk "$LED" 2>/dev/null)
real_x=$?
real_fail=$(bash audit/verify_ledger.sh 2>/dev/null | grep -E '^RESULT FAIL=')
if [ "$real_x" -eq 0 ] && [ "$real_fail" = "RESULT FAIL=0" ]
then ok "真台账放行（$real_out · $real_fail）"
else ng "真台账被误判（exit=$real_x · $real_fail）"; fi

rm -rf "$TMP"
echo
echo "RESULT PASS=$n_pass FAIL=$n_fail"
if [ "$n_fail" -eq 0 ]; then echo "EXIT=0"; exit 0; else echo "EXIT=1"; exit 1; fi
