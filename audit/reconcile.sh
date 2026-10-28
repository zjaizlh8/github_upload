#!/usr/bin/env bash
# §173b 台账 vs 实跑普查：按 dev 多重集对账（贪心消耗）
# 判据：同篇内每个台账 dev 必须在普查里找到一条还没被用掉的同值行；找不到=台账无实跑支撑。
set -u
cd "$(dirname "$0")/.." || exit 1
for p in bsd hl hodge ns pvsnp rh ym; do
  awk -F, -v P="$p" 'NR>1&&$1==P{v=$4+0; print sprintf("%.6e",v)}' matrix/forcedness_ledger.csv | sort > /tmp/L_$p.txt
  awk -F, -v P="$p" '$1==P && $3!="NA"{print sprintf("%.6e",$3+0)}' audit/census.csv | sort > /tmp/C_$p.txt
  comm -23 /tmp/L_$p.txt /tmp/C_$p.txt > /tmp/miss_$p.txt   # 台账有、普查无（未计重数）
  comm -13 /tmp/L_$p.txt /tmp/C_$p.txt > /tmp/extra_$p.txt  # 普查有、台账无
  printf '%-7s 台账 %3d · 普查(带dev) %3d · 台账无支撑 %2d · 普查多出 %2d\n' \
    "$p" "$(wc -l < /tmp/L_$p.txt)" "$(wc -l < /tmp/C_$p.txt)" "$(wc -l < /tmp/miss_$p.txt)" "$(wc -l < /tmp/extra_$p.txt)"
  [ -s /tmp/miss_$p.txt ] && sed 's/^/        台账独有 /' /tmp/miss_$p.txt
done
