#!/usr/bin/env awk -f
# §173e 把 24 行「隐含零」的 dev 由 0.000000e+00 / 0 标成 NA(implicit)
# 依据：这 24 行对应的实跑判决行里【没有 dev= token】（census.csv 的 $3=="NA"），
#       dev=0 是「PASS ⟹ 无残差」推出来的，不是读出来的（G25：声明值必须真在输出里）。
# 识别（两路独立且一致）：① 普查侧 hl 6 + pvsnp 18 行无 dev；② 台账侧按 line 精确命中。
BEGIN{
  FS=","; OFS=",";
  split("375#pi2_1e5 376#pi2_1e6 377#pi2_1e7 378#pi2_1e8 379#pi2_1e9 375#center_6m", a, " ");
  for (i in a) mark["hl" SUBSEP a[i]] = 1;
  for (n = 2; n <= 8; n++) { mark["pvsnp" SUBSEP ("197#n" n)] = 1; mark["pvsnp" SUBSEP ("211#n" n)] = 1 }
  mark["pvsnp" SUBSEP "276#P8a"] = 1; mark["pvsnp" SUBSEP "276#P8b"] = 1;
  mark["pvsnp" SUBSEP "276#P8c"] = 1; mark["pvsnp" SUBSEP "209#P6b"] = 1;
}
NR == 1 { print; next }
{
  k = $1 SUBSEP $2;
  if (k in mark) {
    if (($4 + 0) != 0) printf("REFUSE %s,%s dev=%s 非零，不动\n", $1, $2, $4) > "/dev/stderr";
    else { $4 = "NA(implicit)"; hit++; delete mark[k] }
  }
  print;
}
END{ printf("标记 %d 行 · 未命中 %d 条\n", hit, length(mark)) > "/dev/stderr" }
