#!/usr/bin/env bash
# §172b 台账 line 列的地址有效性：全 198 行逐行问「那一行存在且非空吗」。
set -u
cd "$(dirname "$0")/.." || exit 1
for pair in "bsd:bsd-conjecture" "hl:hl" "hodge:hodge-conjecture" "ns:navier-stokes" "pvsnp:p-vs-np" "rh:riemann-hypothesis" "ym:yang-mills"; do
  p=${pair%%:*}; d=${pair##*:}
  f="$d/paper.zh.md"; [ -f "$f" ] || { echo "$p 缺 $f"; continue; }
  tot=$(awk -F, -v P="$p" 'NR>1&&$1==P' matrix/forcedness_ledger.csv | wc -l)
  awk -F, -v P="$p" 'NR>1&&$1==P{print $2}' matrix/forcedness_ledger.csv \
  | sed 's/#.*//' | sort -u > /tmp/ln_$p.txt
  ok=0; bad=0
  while read -r ln; do
    t=$(sed -n "${ln}p" "$f" 2>/dev/null)
    if [ -n "${t// /}" ]; then ok=$((ok+1)); else bad=$((bad+1)); fi
  done < /tmp/ln_$p.txt
  printf '%-7s 行 %3d · 不同地址 %3d · 有效 %3d · 落空 %3d\n' "$p" "$tot" "$(wc -l < /tmp/ln_$p.txt)" "$ok" "$bad"
done
