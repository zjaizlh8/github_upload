#!/usr/bin/env bash
# §172 从「台账 line 指向的那一行论文原文」里机械抽出判据名，只看 unambiguous 的。
# 判据：该行内恰好命中一个 P<数字> 或 P<数字><字母> 记号 ⟹ 可填；零个或多个 ⟹ 留 NA。
set -u
cd "$(dirname "$0")/.." || exit 1
fill=0; amb=0; none=0; bad=0
awk -F, 'NR>1 && $9 ~ /^NA/{print $1"\t"$2"\t"$4}' matrix/forcedness_ledger.csv | while IFS=$'\t' read -r p ln dev; do
  d=$(case $p in bsd) echo bsd-conjecture;; hl) echo hl;; hodge) echo hodge-conjecture;; ns) echo navier-stokes;; pvsnp) echo p-vs-np;; rh) echo riemann-hypothesis;; ym) echo yang-mills;; esac)
  f="$d/paper.zh.md"
  txt=$(sed -n "${ln}p" "$f" 2>/dev/null)
  [ -z "$txt" ] && { printf 'BADADDR\t%s\t%s\n' "$p" "$ln"; continue; }
  toks=$(printf '%s' "$txt" | grep -oE '\bP[0-9]+[a-c]?\b' | sort -u | tr '\n' ' ')
  cnt=$(printf '%s' "$toks" | wc -w)
  if [ "$cnt" = "1" ]; then printf 'FILL\t%s\t%s\t%s\n' "$p" "$ln" "$toks"
  elif [ "$cnt" = "0" ]; then printf 'NONE\t%s\t%s\n' "$p" "$ln"
  else printf 'AMB\t%s\t%s\t%s\n' "$p" "$ln" "$toks"; fi
done > /tmp/nam_fill.txt
printf '%-10s %s\n' "可填" "$(grep -c '^FILL' /tmp/nam_fill.txt)"
printf '%-10s %s\n' "多义" "$(grep -c '^AMB' /tmp/nam_fill.txt)"
printf '%-10s %s\n' "行内无 P 记号" "$(grep -c '^NONE' /tmp/nam_fill.txt)"
printf '%-10s %s\n' "地址落空" "$(grep -c '^BADADDR' /tmp/nam_fill.txt)"
echo "--- 可填的前 12 条 ---"; grep '^FILL' /tmp/nam_fill.txt | head -12
