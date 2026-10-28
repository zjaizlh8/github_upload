#!/usr/bin/env bash
# fix_tol.sh — 台账 tol 列三修（来源一致化 · 反抄纠正）
#   1) 字面统一到 %.1e 机器格式（0 → 0.0e+00；1e-15 → 1.0e-15；5e-02 → 5.0e-02）
#   2) 3 行 tol=dev 反抄（rh 47#P10 / 47#P11 / 24#P12）⟹ 换成机制界 1.0e-15
#      机制：两条闭式路径各约 20 步运算，舍入上界 γ_n = n·u/(1-u) ≈ 2.2e-15（u=2^-53）；
#            取 1.0e-15 ≈ 4.5 ulp，是「独立实现只差几 ulp」的机制界，不是本次读数。
#      ratio 同步重算 2.2e-16/1.0e-15 = 0.22
#   3) 不动任何 dev / class / criterion
set -euo pipefail
SRC="$1"
DST="$2"

awk -F',' 'BEGIN{OFS=","}
NR==1 {print; next}
{
  t=$5
  # 2) 反抄三行优先
  if ($1=="rh" && ($2=="47#P10" || $2=="47#P11" || $2=="24#P12")) {
      $5="1.0e-15"; $6="0.2"
      print; next
  }
  # 1) 字面统一
  if (t=="0") { $5="0.0e+00"; print; next }
  if (t ~ /^[0-9]+[eE][-+][0-9]+$/) {
      sub(/[eE]/, ".0e", t); $5=t; print; next
  }
  print
}' "$SRC" > "$DST"
echo "WROTE $DST"
