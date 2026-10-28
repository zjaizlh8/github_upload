#!/usr/bin/env awk -f
# §173d 把台账 dev 列换成【实跑输出里的原始字符串】（按计算实际为准）
# 规则：同篇内，取数值最接近的那条 census 原始串；相对差 ≤ 0.5% 才替换；dev=0 的行不动。
# 输出：原 9 列不变，只改第 4 列；同时把「替换/未替换/差多少」报到 stderr。
BEGIN{
  FS=","; OFS=",";
  while ((getline line < "audit/census.csv") > 0) {
    n=split(line, a, ",");
    if (n < 4) continue;
    if (a[3] == "NA") continue;
    k  = a[1] SUBSEP (++cnt[a[1]]);
    cp[k] = a[1]; cd[k] = a[3] + 0; cs[k] = a[3];
  }
  nrep=0; nkeep=0; nflag=0;
}
FNR==1{ next }
{
  p=$1; v=$4+0;
  if (v != 0) {
    best=-1; bestrel=1e9;
    for (i=1; i<=cnt[p]; i++) {
      k = p SUBSEP i;
      dv = cd[k]; if (dv == 0) continue;
      rel = (v - dv) / (dv < 0 ? -dv : dv); if (rel < 0) rel = -rel;
      if (rel < bestrel) { bestrel = rel; best = i; }
    }
    if (best > 0 && bestrel <= 5e-3) {
      k = p SUBSEP best;
      if (cs[k] != $4) { printf("REP %-6s %-14s -> %-14s (rel %.2e)\n", p, $4, cs[k], bestrel) > "/dev/stderr"; nrep++ }
      else nkeep++;
      $4 = cs[k];
    } else {
      printf("FLAG %-6s %-14s 无 0.5%% 内的对应（最近 rel %.2e）\n", p, $4, bestrel) > "/dev/stderr"; nflag++;
    }
  }
  print;
}
END{ printf("替换 %d · 已一致 %d · 无对应 %d\n", nrep, nkeep, nflag) > "/dev/stderr" }
