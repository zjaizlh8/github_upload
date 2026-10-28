# x9_tol.awk —— X9：画像对齐行里「台账 tol」与「输出 tol」的不一致必须全部落在【过程面】
#
# 为什么需要它：census2/outs 里的 tol 是【这台机器真判的界】，台账 tol 是【声明的界】。
# 两侧不一致本身可以合法：过程面（FLOW/TREND）按 §179/X8 主动把 tol 撤回 NA。
# 但若【对象面】也冒出一处不一致，说明声明和代码已经不是同一个界 —— 这正是
# §181 之后 bsd 三条藏了整整一轮的那类事故（census2 陈旧，X2 只自证 dev 不证 tol）。
#
# 判据：对象面（class ∉ {FLOW,TREND,DUP}）且两侧 tol 都是数值 且 相对差 > 1e-12 ⟹ 计一处违规。
# 用法：awk -F, -f audit/x9_tol.awk audit/census2.csv matrix/forcedness_ledger.csv matrix/profile_v2.csv
#       （三个文件按此序给；识别靠【文件名】——census2 / profile_v2 / 其余=台账，
#         这样齿测算在临时路径上的台账副本也能被正确识别，不会被文件名坑掉。）
# 输出：BAD <违规数> TOT <参与比较行数> PROC <落在过程面的不一致数>
BEGIN { FS = "," }
FNR == 1 { next }

FILENAME ~ /census2/ { CT[$1 "|" $2] = $6; next }

FILENAME ~ /profile_v2/ {
    if ($7 == "none") next
    ln = $9
    if (ln == "") next
    key = $1 "|" ln
    if (!(key in LT)) next
    t_out = CT[$1 "|" $2]
    t_led = LT[key]
    if (t_out == "NA" || t_led == "" || t_led == "NA" || t_led == "NA(implicit)") next
    a = t_out + 0; b = t_led + 0
    d = a - b; if (d < 0) d = -d
    m = (a > b) ? a : b
    tot++
    if (m > 0 && d / m > 1e-12) {
        if (LC[key] == "FLOW" || LC[key] == "TREND") proc++
        else { bad++; if (bad <= 5) printf "    [X9] %s[%s] 输出 tol=%s 台账 tol=%s class=%s\n", $1, $2, t_out, t_led, LC[key] }
    }
    next
}

{ LT[$1 "|" $2] = $5; LC[$1 "|" $2] = $3 }

END { printf "BAD %d TOT %d PROC %d\n", bad + 0, tot + 0, proc + 0 }
