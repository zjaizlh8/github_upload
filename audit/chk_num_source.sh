#!/usr/bin/env bash
# §170 论文里的「像读数的数」是否出现在自己的验证器输出里
# 判据：≥5 位有效数字的数（含科学记数法）视为读数；落在 outs/<paper>.txt 里算可复现，否则=无来源。
set -u
cd "$(dirname "$0")/.." || exit 1
tot=0; hit=0
for pair in "bsd:bsd-conjecture" "hl:hl" "hodge:hodge-conjecture" "ns:navier-stokes" "pvsnp:p-vs-np" "rh:riemann-hypothesis" "ym:yang-mills"; do
  p=${pair%%:*}; d=${pair##*:}
  paper="$d/paper.zh.md"; out="audit/outs/$p.txt"
  [ -f "$paper" ] || continue
  # 抓 ≥5 位有效数字的小数（跳过纯整数、跳过 1e-16 之类的容差常数）
  grep -oE '[0-9]+\.[0-9]{4,}(e[-+][0-9]+)?|[0-9]\.[0-9]{3,}e[-+][0-9]{2,}' "$paper" \
    | sort -u > /tmp/nums_$p.txt
  n=0; h=0
  while read -r v; do
    n=$((n+1))
    # 数值比对：把输出里的每个数都规格化成 %.6e 再比
    if awk -v V="$v" 'BEGIN{target=sprintf("%.6e",V+0)} {for(i=1;i<=NF;i++){gsub(/[^0-9eE.+-]/,"",$i); if($i!="" && sprintf("%.6e",$i+0)==target) found=1}} END{exit(found?0:1)}' "$out"; then
      h=$((h+1))
    else
      printf '   %-7s MISS %s\n' "$p" "$v"
    fi
  done < /tmp/nums_$p.txt
  printf '%-7s 读数样 %3d · 可复现 %3d · 无来源 %3d\n' "$p" "$n" "$h" "$((n-h))"
  tot=$((tot+n)); hit=$((hit+h))
done
echo "== 合计 $tot · 可复现 $hit · 无来源 $((tot-hit)) =="
