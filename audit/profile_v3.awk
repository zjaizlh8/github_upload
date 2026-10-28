# profile_v3.awk —— 结构数据画像 v3：把六条轴并进【一张表】，并保证每条划分都能加回 198
#
# 为什么要 v3：v2 只有 join/criterion 两轴。本轮又建了三条独立结构（§184）：
#   · 面（face）      ：tol 有数值 ⟺ 对象面（X8）
#   · 零的种类（zk）  ：几何零 / 计数零 / 隐含零 / 非零残差（§184.7）
#   · 残差型（ra）    ：R 型 dev = k·2^-52（k=m/2^j, m≤64, j≤8）/ A 型（dev 恰 0）/ X 型（都不是）
# 三条都是【对 198 行的划分】，但彼此正交（面 vs 零的种类 vs 二进制可表性）。
# 画像的用处不是好看，是让「还剩什么没被解释」能被机器点出来：每一个 X 型行、
# 每一个 count 零、每一个 NA(addr_dead) 都是一个【尚未被解释的位置】。
#
# 用法：awk -F, -f audit/profile_v3.awk matrix/profile_v2.csv matrix/forcedness_ledger.csv
# 输出列：paper,seq,name,dev,tol,join,class,face,zk,ra,addr,criterion_final
BEGIN { FS = ","; OFS = "," }
FNR == 1 { next }

# 台账只用来补 class（profile_v2 第 10 列已有 ledger_class，这里仍以台账为真源）
FILENAME ~ /forcedness_ledger/ { CLS[$1 "|" $2] = $3; next }

{
    cls = $10; if (cls == "") cls = CLS[$1 "|" $9]

    # ① 面：对象面 / 过程面 / 重复标记
    if (cls == "FLOW" || cls == "TREND") face = "process"
    else if (cls == "DUP")             face = "dup"
    else                               face = "object"

    # ② 零的种类
    dev = $5
    if (dev == "NA(implicit)" || dev == "NA") zk = "implicit"
    else if ((dev + 0) == 0) {
        # 计数零：dev 是「违例个数」不是残差 —— 判据文本自带这个说法
        crit = $13; for (i = 14; i <= NF; i++) crit = crit "," $i
        if (crit ~ /违例|violation|无残差可报|no residual/) zk = "count"
        else zk = "geo"
    } else zk = "nonzero"

    # ③ 残差型 R/A/X：R = dev 能在二进制下写成 k·2^-52（k 为 m/2^j，m≤64，j≤8）
    # 容差必须由【印位数】给：台账 dev 只有 7 位有效数字 ⟹ 相对误差 ≤5e-7；
    # 换算到 m 的单位上，允许 |m − round(m)| ≤ 5e-7·m（第一版写成固定 1e-3 是错的：
    # m 可到 64，1e-3 的绝对窗等于 1e-3 的相对窗 —— 比允许的宽 2000 倍，会把 23 行误判成 R）。
    if (zk != "nonzero") ra = "A"
    else {
        x = (dev + 0) / 2.220446049250313e-16      # 以 2^-52 为单位
        ra = "X"
        for (j = 0; j <= 8; j++) {
            mm = x * (2 ^ j)
            mr = int(mm + 0.5)
            if (mr >= 1 && mr <= 64) {
                dm = mm - mr; if (dm < 0) dm = -dm
                if (dm <= 5e-7 * mm) { ra = "R"; break }
            }
        }
    }

    crit2 = $13; for (i = 14; i <= NF; i++) crit2 = crit2 "," $i
    print $1, $2, $4, dev, $6, $7, cls, face, zk, ra, $12, crit2
}
