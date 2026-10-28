#!/usr/bin/env bash
# tooth_x9.sh —— X9（台账 tol = 输出 tol，对象面不许背离）的齿（2026-10-30 · §184.6）
#
# 三个篡改（都只改一行，其余 197 行原样）：
#   c1 对象面 tol 改成【输出里没有的数】：bsd:5 tol=1.0e-03 ⟹ X9 咬住（X2 也一起咬）
#   c2 对象面 tol 改成【本篇输出里真有的另一个 tol】（bsd 有 1.0e-08/1.0e-12 等）
#      ⟹ 这类篡改 X2 只比 (dev,tol) 对是否在本篇出现过，容易被「同一对恰好在别处也有」放过；
#        X9 是逐行对齐比，必然咬住 —— 这一条专门测 X9 是否【超出 X2 的判别力】。
#   c3 过程面带数值 tol（hl:7 FLOW → 3e+00）：X9 只比两侧都是数值的行，合法放行（不误咬）。
# 通过判据：c1/c2 各 FAIL>=1；c3 不得因 X9 翻红；真台账 FAIL=0。
set -u
cd "$(dirname "$0")/.." || exit 1
LED=matrix/forcedness_ledger.csv
TMP=$(mktemp -d)
n_pass=0; n_fail=0
ok(){ n_pass=$((n_pass+1)); printf '  [PASS] %s\n' "$1"; }
ng(){ n_fail=$((n_fail+1)); printf '  [FAIL] %s\n' "$1"; }

awk -F, 'BEGIN{OFS=","} $1=="bsd"&&$2=="5"{$5="1.0e-03"} {print}' "$LED" > "$TMP/c1.csv"
awk -F, 'BEGIN{OFS=","} $1=="bsd"&&$2=="5"{$5="1.0e-08"} {print}' "$LED" > "$TMP/c2.csv"
awk -F, 'BEGIN{OFS=","} $1=="hl"&&$2=="7"{$5="3e+00"} {print}' "$LED" > "$TMP/c3.csv"

echo "[齿] X9 · 台账 tol = 输出 tol（产物 audit/x9_tol.awk）"
for c in c1 c2; do
    x9out=$(awk -f audit/x9_tol.awk audit/census2.csv "$TMP/$c.csv" matrix/profile_v2.csv 2>/dev/null)
    bad=$(printf '%s\n' "$x9out" | awk '/^BAD/{print $2}')
    fails=$(bash audit/verify_ledger.sh "$TMP/$c.csv" 2>/dev/null | grep -c '^  \[FAIL\]')
    if [ "${bad:-0}" -ge 1 ]
    then ok "$c 被 X9 咬住（$x9out · verify_ledger FAIL 行=$fails）"
    else ng "$c 逃过 X9（$x9out · verify_ledger FAIL 行=$fails）"; fi
done

x9c3=$(awk -f audit/x9_tol.awk audit/census2.csv "$TMP/c3.csv" matrix/profile_v2.csv 2>/dev/null)
c3_bad=$(printf '%s\n' "$x9c3" | awk '/^BAD/{print $2}')
if [ "${c3_bad:-1}" = "0" ]
then ok "c3 不误咬（过程面带数值 tol 合法：$x9c3）"
else ng "c3 被误咬（$x9c3）"; fi

real_out=$(awk -f audit/x9_tol.awk audit/census2.csv "$LED" matrix/profile_v2.csv 2>/dev/null)
real_fail=$(bash audit/verify_ledger.sh 2>/dev/null | grep -E '^RESULT FAIL=')
if [ "$real_out" = "BAD 0 TOT 161 PROC 0" ] && [ "$real_fail" = "RESULT FAIL=0" ]
then ok "真台账放行（$real_out · $real_fail）"
else ng "真台账被误判（$real_out · $real_fail）"; fi

rm -rf "$TMP"
echo
if [ "$n_fail" -eq 0 ]
then echo "RESULT PASS=$n_pass FAIL=0"; echo "EXIT=0"; exit 0
else echo "RESULT PASS=$n_pass FAIL=$n_fail"; echo "EXIT=1"; exit 1; fi
