#!/usr/bin/env bash
# §170b 对 §170 的 MISS 数做「仓内是否有任何来源」的二分：论文自身不算来源。
set -u
cd "$(dirname "$0")/.." || exit 1
for pair in "bsd:bsd-conjecture" "hl:hl" "hodge:hodge-conjecture" "ns:navier-stokes" "pvsnp:p-vs-np" "rh:riemann-hypothesis" "ym:yang-mills"; do
  p=${pair%%:*}; d=${pair##*:}
  [ -f "/tmp/nums_$p.txt" ] || continue
  n=0; h=0
  while read -r v; do
    if awk -v V="$v" 'BEGIN{target=sprintf("%.6e",V+0)} {for(i=1;i<=NF;i++){gsub(/[^0-9eE.+-]/,"",$i); if($i!="" && sprintf("%.6e",$i+0)==target) found=1}} END{exit(found?0:1)}' "audit/outs/$p.txt"; then continue; fi
    n=$((n+1))
    # 在仓内找「原样字面」出现（排除该篇论文自身与台账）
    hits=$(grep -rlF "$v" --include=*.c --include=*.md --include=*.txt --include=*.csv . 2>/dev/null \
           | grep -v "^./$d/" | grep -v "forcedness_ledger.csv" | grep -v "AUDIT_2026" | head -3)
    if [ -n "$hits" ]; then
      h=$((h+1)); printf '   %-7s %-16s 仓内另有来源: %s\n' "$p" "$v" "$(echo "$hits" | tr '\n' ' ')"
    else
      printf '   %-7s %-16s ★全仓无来源\n' "$p" "$v"
    fi
  done < <(grep -oE '[0-9]+\.[0-9]{4,}(e[-+][0-9]+)?|[0-9]\.[0-9]{3,}e[-+][0-9]{2,}' "$d/paper.zh.md" | sort -u)
  printf '%-7s 无输出来源 %3d · 其中仓内另有来源 %3d · 全仓无来源 %3d\n' "$p" "$n" "$h" "$((n-h))"
done
