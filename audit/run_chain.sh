#!/usr/bin/env bash
# 七支验证器 + 审计器串跑，落 audit/outs/*.txt。判 FAIL 用结构匹配，不用词匹配。
set -u
cd "$(dirname "$0")/.." || exit 1
ok=1
run_one() {
  p="$1"; d="$2"
  gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o "bin/verify_$p" "$d/code/verify_$p.c" -lm 2>/dev/null
  "./bin/verify_$p.exe" > "audit/outs/$p.txt" 2>&1
  r=$(grep '^RESULT' "audit/outs/$p.txt" | tail -1)
  f=$(printf '%s' "$r" | sed 's/.*FAIL=\([0-9]*\).*/\1/')
  printf '%-7s %s\n' "$p" "$r"
  [ "$f" = "0" ] || ok=0
}
run_one bsd   bsd-conjecture
run_one hl    hl
run_one hodge hodge-conjecture
run_one ns    navier-stokes
run_one pvsnp p-vs-np
run_one rh    riemann-hypothesis
run_one ym    yang-mills
gcc -Wall -Wextra -O2 -o bin/verify_audit audit/verify_audit.c -lm 2>/dev/null
./bin/verify_audit.exe matrix/audit_layer_gap.csv > audit/outs/audit.txt 2>&1
a=$(grep '^RESULT' audit/outs/audit.txt | tail -1)
af=$(printf '%s' "$a" | sed 's/.*FAIL=\([0-9]*\).*/\1/')
printf 'audit   %s\n' "$a"
[ "$af" = "0" ] || ok=0
# 齿：台账 tol 的机制界（bsd:46/74 的 8.0e-10）必须咬得住 —— 拿掉它，tol 就退回圆整手给
# 作者确认: 2026-10-29 · YG（新增一支齿，不改任何既有判据与判定方向；permit=abe9d92858e7b26e）
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o bin/tooth_tol_mech audit/tooth_tol_mech.c -lm 2>/dev/null
./bin/tooth_tol_mech.exe > audit/outs/tooth_tol_mech.txt 2>&1
t=$(grep '^RESULT' audit/outs/tooth_tol_mech.txt | tail -1)
tf=$(printf '%s' "$t" | sed 's/.*FAIL=\([0-9]*\).*/\1/')
printf 'tooth   %s\n' "$t"
[ "$tf" = "0" ] || ok=0
# 齿：X8（tol 有值 ⟺ 对象面）必须咬得住三种篡改 —— 拿掉它，过程面就能带 tol、面桶就能漏格
# 作者确认: 2026-10-29 · YG（新增一支齿，不改任何既有判据与判定方向；permit=dbd25dea34d96e67）
bash audit/tooth_x8.sh > audit/outs/tooth_x8.txt 2>&1
x=$(grep '^RESULT' audit/outs/tooth_x8.txt | tail -1)
xf=$(printf '%s' "$x" | sed 's/.*FAIL=\([0-9]*\).*/\1/')
printf 'tooth8  %s\n' "$x"
[ "$xf" = "0" ] || ok=0

# 齿：X9（台账 tol = 输出 tol）必须咬得住对象面的任何背离 —— 拿掉它，声明与代码可以各判各的界
# 作者确认: 2026-10-30 · YG（新增一支齿，不改任何既有判据与判定方向；AUDIT §184）
bash audit/tooth_x9.sh > audit/outs/tooth_x9.txt 2>&1
y=$(grep '^RESULT' audit/outs/tooth_x9.txt | tail -1)
yf=$(printf '%s' "$y" | sed 's/.*FAIL=\([0-9]*\).*/\1/')
printf 'tooth9  %s\n' "$y"
[ "$yf" = "0" ] || ok=0
echo "ALL_FAIL_ZERO=$ok"
