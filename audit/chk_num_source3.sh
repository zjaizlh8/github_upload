#!/usr/bin/env bash
# §170c 修正版：先把输出里的数按 token 抽出（不再逐字段清洗，避免 e 被当成指数），再比 %.6e。
set -u
cd "$(dirname "$0")/.." || exit 1
TOT=0; HIT=0; NOWHERE=0
for pair in "bsd:bsd-conjecture" "hl:hl" "hodge:hodge-conjecture" "ns:navier-stokes" "pvsnp:p-vs-np" "rh:riemann-hypothesis" "ym:yang-mills"; do
  p=${pair%%:*}; d=${pair##*:}
  grep -oE '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?' "audit/outs/$p.txt" \
    | awk '{print sprintf("%.6e",$1+0)}' | sort -u > /tmp/onum_$p.txt
  n=0; h=0
  for v in $(grep -oE '[0-9]+\.[0-9]{4,}(e[-+][0-9]+)?|[0-9]\.[0-9]{3,}e[-+][0-9]{2,}' "$d/paper.zh.md" | sort -u); do
    n=$((n+1))
    t=$(awk -v V="$v" 'BEGIN{print sprintf("%.6e",V+0)}')
    if grep -qxF "$t" /tmp/onum_$p.txt; then h=$((h+1)); continue; fi
    hits=$(grep -rlF "$v" --include=*.c --include=*.md --include=*.txt --include=*.csv . 2>/dev/null \
           | grep -v "^./$d/" | grep -v "forcedness_ledger.csv" | grep -v "AUDIT_2026" | grep -v "README.md" | head -2)
    if [ -n "$hits" ]; then
      printf '   %-7s %-16s 间接: %s\n' "$p" "$v" "$(echo "$hits" | tr '\n' ' ')"
    else
      printf '   %-7s %-16s ★全仓无来源\n' "$p" "$v"; NOWHERE=$((NOWHERE+1))
    fi
  done
  printf '%-7s 读数样 %3d · 输出可复现 %3d · 缺口 %3d\n' "$p" "$n" "$h" "$((n-h))"
  TOT=$((TOT+n)); HIT=$((HIT+h))
done
echo "== 合计 $TOT · 可复现 $HIT · 缺口 $((TOT-HIT)) · 全仓无来源 $NOWHERE =="
