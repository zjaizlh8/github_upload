# fix_ratio.awk —— 把台账第 6 列 ratio 改成【由 dev/tol 重算】的机器串
#
# 为什么必须重算（§174.8）：
#   1. 陈旧：§174 把 rh 5 行的 tol 从 1e-12/1e-15 改成 1.000000e-08，ratio 没跟着走
#      ⟹ 表内自相矛盾（ratio 说 0.0020，dev/tol 说 1.998e-07，差 10^4）。
#   2. 截位：另有 7 行把 1.578e-03 写成 0.0016（只留 2 位有效 ⟹ 差 1.4%）。
#   3. 格式两义：零记成 "0.0000"（104 行）与 "0"（34 行）两种写法。
#
# 口径（与 dev/tol 两列同一口径）：
#   dev = NA(implicit)  ⟹ ratio = NA(implicit)（零是 PASS 推的，不是读数）
#   tol == 0            ⟹ ratio = 0.000000e+00（精确声明，分母为 0 不定义比值）
#   其余                ⟹ ratio = %.6e(dev/tol)
#
# 用法：awk -F, -f fix_ratio.awk 台账.csv > 新.csv（由外层脚本补表头并替换）

BEGIN { FS = ","; OFS = "," }
NR == 1 { print; next }
{
    if ($4 == "NA(implicit)") {
        $6 = "NA(implicit)"; nimpl++
    } else {
        d = $4 + 0; t = $5 + 0
        if (t == 0) { $6 = "0.000000e+00"; nzero++ }
        else { $6 = sprintf("%.6e", d / t); nreal++ }
    }
    print
}
END {
    printf "RATIO-REBUILD implicit=%d zero=%d real=%d\n", nimpl, nzero, nreal > "/dev/stderr"
}
