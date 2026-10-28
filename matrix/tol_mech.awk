# tol_mech.awk — 台账 tol 列的机制体检（2026-10-29）
# 用法：
#   awk -F',' -v MODE=report -f tol_mech.awk forcedness_ledger.csv        # 只出计数报告
#   awk -F',' -v MODE=list   -f tol_mech.awk forcedness_ledger.csv        # 逐行明细
#   awk -F',' -v MODE=apply  -f tol_mech.awk forcedness_ledger.csv > out   # 只把 R1 反抄行降 TREND
#
# 三规则（只认这三件，不拍脑袋）：
#   R1 反抄自比：tol ≈ 2·dev（±1e-6 相对）⟹ tol 抄自本次读数 ⟹ 永远不可能 FAIL
#                （与 G53「总和=行数」同型：自比）⟹ 必降 TREND
#   R2 机制界  ：tol/u ∈ [0.5,100]（u = 2^-52 = 2.221e-16，即「几 ulp」机器界）⟹ OK_MACH
#                或 criterion 明写 1/N、1/h、h^p、收敛 ⟹ OK_ORDER
#   R3 其余    ：NO（圆整界，未导出机制）⟹ 只报不判
BEGIN { OFS=","; U=2.220446049250313e-16; nR1=0; nOK=0; nORD=0; nNO=0; nZ=0 }
NR==1 { if (MODE=="apply") print; next }
{
  dev=$4+0; tol=$5+0; crit=$9
  st="NO"
  # 豁免：dev 严格为 0 ⟹ 整数域/恒等，无 tol 可绑（G54：整数域严格 0 本身就是律）
  if (dev==0) st="OK_ZERO"
  if (st=="NO" && dev>0 && tol>0) { r=tol/dev; if (r>1.9 && r<2.1) st="R1" }
  if (st=="NO") {
    if (tol>0) { k=tol/U; if (k>=0.5 && k<=100) st="OK_MACH" }
    if (st=="NO" && crit ~ /(1\/N|1\/h|h\^|收敛|convergence)/) st="OK_ORDER"
  }
  if (st=="R1") nR1++; else if (st=="OK_MACH") nOK++; else if (st=="OK_ORDER") nORD++; else if (st=="OK_ZERO") nZ++; else nNO++
  if (MODE=="list" && st!="OK_MACH" && st!="OK_ZERO") printf "#%-8s %-12s dev=%-14s tol=%-12s ratio=%-12s | %s\n", st, $1":"$2, $4, $5, $6, substr(crit,1,58)
  if (MODE=="apply") {
    if (st=="R1") { $3="TREND"; $5="NA"; $6="NA"; $9=crit " 〔tol 取在 dev 的 2 倍附近（无机制手给界）⟹ 降 TREND（R1 机制体检）〕" }
    print
  }
}
END {
  if (MODE!="apply") printf "TOTAL=%d · R1_反抄=%d · OK_MACH=%d · OK_ORDER=%d · OK_ZERO=%d · NO_MECH=%d\n", NR-1, nR1, nOK, nORD, nZ, nNO
}
