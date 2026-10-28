# profile_v2.awk —— 结构数据画像 v2：把「输出侧」（census2）与「声明侧」（台账）按
# 保序最大匹配（LCS）对齐，逐行给出身份、对象名、地址状态与最终判据文本。
#
# 用法（由 audit/build_profile.sh 调用）：
#   awk -F, -f audit/profile_v2.awk audit/census2.csv matrix/addr_state.csv matrix/forcedness_ledger.csv
#
# 三档 join：
#   pair = (dev,tol) 在本篇内两侧都只出现一次 ⟹ 配对被迫
#   lcs  = 落在保序最大匹配的路径上 ⟹ 行序约定成立（但 dev 重复时配对不唯一）
#   none = 不在路径上 ⟹ 本篇存在插入/删除（rh 的 part1/part2 混排、bsd 的三条新增）
#
# 名字传递规则（不许猜）：
#   name_ok = pair  或  (lcs 且 census 侧 dev≠0 且该 dev 在本篇只出现一次)
#   其余行只允许用论文自己的话（paper:<line> <text>，地址必须 valid）或标 NA(addr_dead)。
BEGIN { FS = ","; NPAPER = 7; split("bsd hl hodge ns pvsnp rh ym", PL, " ") }

function norm(x) {
    # 两种「零」不是同一个对象：Z0 = 输出里真印出来的零（读出的），
    # ZI = 判决行里没有 dev token、零是由 PASS 推出来的（推断的，§173）。
    # 旧版把两者压成同一个 "Z" ⟹ LCS 会把「读出的零」配到「推断的零」上
    # （实测 hl 21–26 与 29–35 互换、pvsnp 1–9 与 19–27 互换）⟹ 已分家。
    if (x == "NA") return "ZI"
    if (x == "NA(implicit)") return "ZI"
    if (x == "") return "ZI"
    if (x ~ /^0\.0+e\+00$/) return "Z0"
    if (x == "0") return "Z0"
    return x
}
function iszero(x) { return (norm(x) == "Z0" || norm(x) == "ZI") }

FILENAME ~ /census2/ {
    if ($0 ~ /paper,/) next
    p = $1
    CNT[p]++
    i = CNT[p]
    CD[p, i] = norm($5)
    CT[p, i] = $6
    CN[p, i] = $4
    CS[p, i] = $3
    CRAW[p, i] = $5
    next
}
FILENAME ~ /addr_state/ {
    p = $1; ln = $2
    AD[p, ln] = $3
    AT[p, ln] = $4
    next
}
{
    if ($0 ~ /^paper,/) next
    p = $1
    LNT[p]++
    i = LNT[p]
    LL[p, i] = $2
    LC[p, i] = $3
    LD[p, i] = norm($4)
    LT[p, i] = $5
    LI[p, i] = $8
    LX[p, i] = $9
    LRAW[p, i] = $4
}

END {
    printf "paper,seq,src,name,dev,tol,join,ledger_row,ledger_line,ledger_class,ledger_inv,addr,criterion_final\n"

    for (k = 1; k <= NPAPER; k++) {
        p = PL[k]
        na = LNT[p]; nb = CNT[p]

        # ---- dev 值在本篇内的出现次数（两侧分别数）----
        for (i = 1; i <= na; i++) LCNT[p, LD[p, i]]++
        for (j = 1; j <= nb; j++) CCNT[p, CD[p, j]]++

        # ---- 保序最大匹配（LCS）· 允许整篇相位偏移 ----
        # 理由（实测）：rh 的台账行序是输出序的**循环移位**（后 10 行被写在前面），
        # 纯 LCS 对循环移位只能保住较长连续段（12/22）。相位 r 是自由参数，
        # 由「取使 LCS 最大的 r」定出 ⟹ 不是拟合，是补一个本来就在的自由度。
        best = -1; bestr = 0
        for (r = 0; r < na; r++) {
            for (i = 0; i <= na; i++) for (j = 0; j <= nb; j++) D[p, i, j] = 0
            for (i = 1; i <= na; i++) {
                ii = ((i - 1 + r) % na) + 1
                for (j = 1; j <= nb; j++) {
                    if (LD[p, ii] == CD[p, j]) D[p, i, j] = D[p, i - 1, j - 1] + 1
                    else D[p, i, j] = (D[p, i - 1, j] > D[p, i, j - 1]) ? D[p, i - 1, j] : D[p, i, j - 1]
                }
            }
            if (D[p, na, nb] > best) { best = D[p, na, nb]; bestr = r }
        }
        PHASE[p] = bestr; lcs = best; LCSLEN[p] = lcs
        for (i = 0; i <= na; i++) for (j = 0; j <= nb; j++) D[p, i, j] = 0
        for (i = 1; i <= na; i++) {
            ii = ((i - 1 + bestr) % na) + 1
            for (j = 1; j <= nb; j++) {
                if (LD[p, ii] == CD[p, j]) D[p, i, j] = D[p, i - 1, j - 1] + 1
                else D[p, i, j] = (D[p, i - 1, j] > D[p, i, j - 1]) ? D[p, i - 1, j] : D[p, i, j - 1]
            }
        }
        i = na; j = nb
        while (i > 0 && j > 0) {
            ii = ((i - 1 + bestr) % na) + 1
            if (LD[p, ii] == CD[p, j]) { MATCH[p, ii] = j; USED[p, j] = ii; i--; j-- }
            else if (D[p, i - 1, j] >= D[p, i, j - 1]) i--
            else j--
        }

        # ---- 逐台账行定 join 档与名字可否传递 ----
        for (i = 1; i <= na; i++) {
            j = MATCH[p, i]
            lv = LD[p, i]
            if (j == 0) { JK[p, i] = "none"; NOK[p, i] = 0 }
            else {
            pair_ok = (LCNT[p, lv] == 1 && CCNT[p, lv] == 1)
            if (pair_ok) JK[p, i] = "pair"
            else JK[p, i] = "lcs"
            if (pair_ok) NOK[p, i] = 1
            else if (!iszero(CD[p, j]) && CCNT[p, CD[p, j]] == 1) NOK[p, i] = 1
            else NOK[p, i] = 0
            }
        }

        # ---- 按 census 行序输出画像行 ----
        for (j = 1; j <= nb; j++) {
            i = USED[p, j]
            if (i == 0) {
                printf "%s,%d,%s,%s,%s,%s,%s,,%s,,,,\n", p, j, CS[p, j], CN[p, j], CRAW[p, j], CT[p, j], "none", "NA"
                ST[p, "none"]++
                continue
            }
            ln = LL[p, i]
            ad = (ln in AD[p]) ? AD[p, ln] : (AD[p, ln ""] != "" ? AD[p, ln] : "oob")
            cf = LX[p, i]
            if (cf ~ /^NA\(zero-dev\)/ || cf ~ /^line[0-9]+ points/) {
                if (NOK[p, i]) cf = "sigma(" p ":" j ") " CN[p, j]
                else if (ad == "valid") cf = "paper:" ln " " AT[p, ln]
                else cf = "NA(addr_dead)"
            }
            printf "%s,%d,%s,%s,%s,%s,%s,%d,%s,%s,%s,%s,%s\n",
                   p, j, CS[p, j], CN[p, j], CRAW[p, j], CT[p, j], JK[p, i], i, ln, LC[p, i], LI[p, i], ad, cf
            ST[p, JK[p, i]]++
            if (NOK[p, i]) NP_OK++
            if (JK[p, i] == "lcs" && !NOK[p, i]) NP_ORDER++
        }
    }

    # ---- 统计（stderr，人看）----
    printf "\n[profile v2] 逐篇 join 档（pair/lcs/none）与 LCS 长度\n" > "/dev/stderr"
    tp = 0; tl = 0; tn = 0
    for (k = 1; k <= NPAPER; k++) {
        p = PL[k]
        printf "  %-6s LCS=%2d/%2d · pair=%2d lcs=%2d none=%2d · 相位=%d\n", p, LCSLEN[p], CNT[p], ST[p, "pair"], ST[p, "lcs"], ST[p, "none"], PHASE[p] > "/dev/stderr"
        tp += ST[p, "pair"]; tl += ST[p, "lcs"]; tn += ST[p, "none"]
    }
    printf "  %-6s pair=%d lcs=%d none=%d · 名字可传=%d · 只能靠行序（不可传名）=%d\n",
           "合计", tp, tl, tn, NP_OK, NP_ORDER > "/dev/stderr"

    # ---- 相位表（X7f 的机器可读来源）----
    if (PHASEFILE != "") {
        printf "paper,phase,lcs,n\n" > PHASEFILE
        for (k = 1; k <= NPAPER; k++) {
            p = PL[k]
            printf "%s,%d,%d,%d\n", p, PHASE[p], LCSLEN[p], CNT[p] >> PHASEFILE
        }
        close(PHASEFILE)
    }
}
