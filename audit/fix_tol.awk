#!/usr/bin/env awk -f
# §174 把台账 tol 列换成【实跑输出里的原始串】
# 规则（逐条可核，不确定就不动）：
#   R1 dev 为 NA(implicit) ⟹ tol 也记 NA(implicit)（实跑那一行既无 dev 也无 tol）
#   R2 dev 为非零数值 ⟹ 同篇内按 dev 精确（%.6e）找候选；唯一候选则抄其 tol；
#      多候选但 tol 全同则抄；否则 FLAG 不动
#   R3 dev 为 0（实跑印出 0.000e+00）⟹ 候选多，只在【该篇零行的 tol 全同】时抄，否则不动
BEGIN{ FS=","; OFS=",";
  while ((getline ln < "audit/census.csv") > 0) {
    n = split(ln, a, ","); if (n < 4) continue;
    p = a[1]; key = p SUBSEP (a[3] == "NA" ? "IMP" : sprintf("%.6e", a[3] + 0));
    tolstr[key] = tolstr[key] "|" ((a[4] == "NA") ? "IMP" : a[4]);
    cnt[key]++;
  }
}
NR == 1 { print; next }
{
  p = $1; d = $4;
  if (d == "NA(implicit)") { if ($5 != "NA(implicit)") { $5 = "NA(implicit)"; r1++ } print; next }
  key = p SUBSEP (($4 + 0) == 0 ? sprintf("%.6e", 0) : sprintf("%.6e", $4 + 0));
  if (!(key in cnt)) { print; next }
  s = tolstr[key]; sub(/^\|/, "", s);
  k = split(s, ta, "|"); uniq = 1;
  for (i = 2; i <= k; i++) if (ta[i] != ta[1]) uniq = 0;
  if (uniq) {
    if ($5 != ta[1]) { printf("REP %-6s dev=%-14s tol %-10s -> %-10s\n", p, $4, $5, ta[1]) > "/dev/stderr"; r2++ }
    else r3++;
    $5 = ta[1];
  } else {
    printf("FLAG %-6s dev=%-14s tol=%-10s 候选 %d 条且 tol 不同（%s）\n", p, $4, $5, cnt[key], s) > "/dev/stderr"; r4++;
  }
  print;
}
END{ printf("R1 隐含标记 %d · R2 替换 %d · R3 已一致 %d · R4 不动 %d\n", r1, r2, r3, r4) > "/dev/stderr" }
