# fill_sheet.awk —— 把「地址落空」的台账行整理成一张填名清单（供作者 5 分钟收口）
# 用法：awk -F, -f audit/fill_sheet.awk matrix/profile_v2.csv audit/census2.csv
# 输出：paper,seq,name,dev,tol,ledger_row,class,slice,n_cand,candidates
#   candidates = 同篇内 (dev,tol) 相同的 census 行，最多列 6 个 ⟹ 作者只需二选一/多选一
BEGIN { FS = "," }
FILENAME ~ /census2/ {
    if ($0 ~ /paper,/) next
    p = $1; k = p SUBSEP $5 SUBSEP $6
    CL[k] = CL[k] (CL[k] == "" ? "" : "|") $2 ":" substr($4, 1, 34)
    CC[k]++
    next
}
{
    if ($0 ~ /^paper,/) {
        print "paper,seq,name,dev,tol,ledger_row,class,slice,n_cand,candidates"
        next
    }
    if ($13 != "NA(addr_dead)") next
    # profile 列：1 paper 2 seq 3 src 4 name 5 dev 6 tol 7 join 8 ledger_row 9 ledger_line 10 class 11 inv 12 addr 13 crit
    k = $1 SUBSEP $5 SUBSEP $6
    printf "%s,%s,%s,%s,%s,%s,%s,%s,%d,%s\n", $1, $2, $4, $5, $6, $8, $10, $11, CC[k] + 0, CL[k]
}
