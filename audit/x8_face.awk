# x8_face.awk —— X8：`tol` 有数值 ⟺ 该行在【对象面】（2026-10-29 · §179.7）
#
# 口径（用户定案）：
#   对象面 = FORCED / FORCED_eps / FORCED_trunc  ⟹ tol 必须是数值（机制或手给都得给，否则这行没法判）
#   过程面 = FLOW / TREND                        ⟹ tol 必须是 NA（过程只报走到哪儿，没有界可给）
#   DUP                                          ⟹ 重复标记，不属两面，不参与本判据（只计数）
#
# 这一个文件同时被两处调用：
#   ① audit/verify_ledger.sh  ⟹ 出 X8a / X8b / X8c 三行判定
#   ② audit/tooth_x8.sh       ⟹ 反例对照（改一行 class 或 tol ⟹ 必须 exit 1）
#
# 输出（stdout 一行，便于 shell 解析）：
#   X8 <obj> <obj_num> <proc> <proc_na> <dup> <unk> <bad>
# 退出码：bad + unk == 0 ⟹ 0；否则 1
BEGIN { FS="," }
FNR == 1 { next }
{
    tot++
    cls = $3; tol = $5
    isnum = (tol != "NA" && tol != "")
    if (cls == "FORCED" || cls == "FORCED_eps" || cls == "FORCED_trunc") {
        obj++
        if (isnum) obj_num++
        else { bad++; printf "X8-bad 对象面缺 tol: %s:%s class=%s tol=%s\n", $1, $2, cls, tol > "/dev/stderr" }
    } else if (cls == "FLOW" || cls == "TREND") {
        proc++
        if (!isnum) proc_na++
        else { bad++; printf "X8-bad 过程面带 tol: %s:%s class=%s tol=%s\n", $1, $2, cls, tol > "/dev/stderr" }
    } else if (cls == "DUP") {
        dup++
    } else {
        unk++; printf "X8-bad 未归面桶: %s:%s class=%s\n", $1, $2, cls > "/dev/stderr"
    }
}
END {
    printf "X8 %d %d %d %d %d %d %d\n", obj, obj_num, proc, proc_na, dup, unk, bad
    printf "X8 覆盖：对象面 %d + 过程面 %d + DUP %d = %d（缺口 %d 行）\n", obj, proc, dup, obj + proc + dup, tot - obj - proc - dup > "/dev/stderr"
    exit ((bad + unk) > 0) ? 1 : 0
}
