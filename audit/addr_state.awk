# addr_state.awk —— 台账 line 列的地址状态（valid / blank / oob）
# 用法：awk -F, -f audit/addr_state.awk paper.zh.md matrix/forcedness_ledger.csv
# 输出：paper,line,state,text（text 截 60 字，只给 valid 行；供 profile_v2 用）
# 说明：判据 = 「该行号上有没有非空正文」。空行/越界都算 dead（§172：68 个落空）。
BEGIN { FS = "," }
FNR == NR {
    # 第一个文件 = 论文正文（逐行读入）
    pl[FNR] = $0
    pn = FNR
    next
}
{
    if ($1 == "paper") next
    # line 列可能是 `209#P6b` / `276#P8a` 这类「行号#标签」写法：取前导数字，不是整列转数
    if (match($2, /^[0-9]+/)) ln = substr($2, 1, RLENGTH) + 0
    else ln = 0
    if (ln < 1) { state = "oob" }
    else if (ln > pn) { state = "oob" }
    else {
        t = pl[ln]
        gsub(/^[ \t]+/, "", t); gsub(/[ \t]+$/, "", t)
        if (t == "") state = "blank"
        else if (t ~ /^#/) state = "blank"
        else state = "valid"
    }
    if (state == "valid") {
        txt = t
        gsub(/,/, ";", txt)
        if (length(txt) > 60) txt = substr(txt, 1, 60)
        printf "%s,%s,%s,%s\n", $1, $2, state, txt
    } else {
        printf "%s,%s,%s,%s\n", $1, $2, state, "-"
    }
}
