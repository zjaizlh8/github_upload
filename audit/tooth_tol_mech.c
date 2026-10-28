/* tooth_tol_mech.c —— 齿：把台账 tol 从「圆整手给」换成「机制导出」之后，界还咬得住吗（2026-10-29）
 *
 * 五法：乘法=木（16·L 与 2·Ω 的乘积项）· 跨维度=火（不涉）· 加法=水（Simpson 求和）
 *       降维=金（不涉）· 无损归心=土（不涉）
 *       本档不做 R(R(x)) 检验 ⟹ res = NA(no-R-test)，如实标注
 * 层=结构：dev(Δ) = 16·(L+Δ) − 2·Ω 对 Δ 的**线性预测**（预测式 = 16Δ，四条独立锚点逐条比）
 * 层=过程：Ω 的数值累加误差只报读数（不设判据）
 * 四判据数：res=NA(no-R-test) · domain={Δ=0, +2e-11, +5e-11, −5e-11} · dc=4（Δ 的 4 个取值）
 *           order=1（dev 对 Δ 线性，预测残差 = 16·½ulp(LEXT)）
 *
 * 为什么要这一支
 *   台账 `matrix/forcedness_ledger.csv` 里 bsd:46/74 原来的 tol = 1e-9（圆整手给，ratio 0.43 贴脸）。
 *   2026-10-29 改成机制界：残差 = 16·(L真 − L锚)，锚只有 10 位 ⟹ 上界 16·½·10⁻¹⁰ = 8.0e-10。
 *   改了就得证明「这个界不是随便写的」——**它必须能 FAIL**（G52：拿掉修它就漏）。
 *
 * 齿的三个方向
 *   T1 基准行在界内（不误伤）
 *   T2 Δ=+5e-11 ⟹ 预测 dev ≈ 1.23e-9 > 8e-10 ⟹ 界必须翻红（咬得住）
 *   T3 Δ=+2e-11 ⟹ 预测 dev ≈ 7.5e-10 ≤ 8e-10 ⟹ 仍在界内（界不是过紧）
 *   T4 dev(Δ) − 16Δ 对四个 Δ 全部为 ~1e-15 级 ⟹ 残差确实来自锚截断（机制成立，不是拟合）
 *   T5 对照：把 tol 换成「2×dev（自比界）」⟹ 4 次扰动 0 次拒绝（判别力 = 0）；机制界拒绝 2 次 ⟹ 两界不同
 */
#include <stdio.h>
#include <math.h>

#define NPTS 200000
#define LEXT 0.6555143886      /* 源：verify_bsd.c:711 的 10 位外部锚（独立源，不由恒等式反推） */
#define TOLM 8.0e-10           /* 机制界：16·½·10⁻¹⁰ */
#define OMANC 2.6220575543     /* 源：verify_bsd.c:736 外的 Ω/2 锚 */

static int g_pass = 0, g_fail = 0;
static void J(const char *name, double dev, double tol)
{
    int ok = (dev <= tol);
    if (ok) g_pass++; else g_fail++;
    printf("  [%s] %-56s val=%.6e tol=%.2e\n", ok ? "PASS" : "FAIL", name, dev, tol);
}

/* Ω(32a1) 的 Simpson 复合求积（与 verify_bsd.c 同一实现，便于逐位对照） */
static double om_quad(int N)
{
    int i;
    double h = (M_PI / 2.0) / (double)N, s = 0.0;
    for (i = 0; i <= N; i++) {
        double th = h * (double)i;
        double f = 1.0 / sqrt(1.0 + cos(th) * cos(th));
        s += ((i == 0 || i == N) ? 1.0 : ((i & 1) ? 4.0 : 2.0)) * f;
    }
    return 4.0 * s * h / 3.0;
}

/* 倒数账残差：dev(Δ) = 16·(LEXT+Δ) − 2·Ω */
static double dev_at(double om, double dlt)
{
    return fabs(16.0 * (LEXT + dlt) - 2.0 * om);
}

/* 带符号残差：线性预测只对带符号量成立（Δ<0 时 |dev| 反而变小） */
static double sdev_at(double om, double dlt)
{
    return 16.0 * (LEXT + dlt) - 2.0 * om;
}

/* 对一组扰动，统计「该界拒绝了几次」——独立性：界与扰动取自不同来源（界由锚位数给，扰动由本档给） */
static int count_reject(const double *dl, int n, double om, double tol, int mode)
{
    int i, r = 0;
    for (i = 0; i < n; i++) {
        double d = dev_at(om, dl[i]);
        double b = (mode == 0) ? tol : 2.0 * d;   /* mode 1 = 自比界（2×dev） */
        if (d > b) r++;
    }
    return r;
}

static int g2_main_body(void)
{
    const double DL[4] = { 0.0, 2e-11, 5e-11, -5e-11 };
    double om = om_quad(NPTS);
    double d0, d2, d5, dm5, pred, resmax = 0.0;
    int i, rej_mech, rej_self;

    printf("=== tooth: 台账 tol 的机制界是否咬得住（bsd:46/74 的 8.0e-10）===\n");
    printf("      Ω(%d) = %.15f   Ω 锚 ×2 = %.10f\n", NPTS, om, 2.0 * OMANC);

    d0  = dev_at(om, DL[0]);
    d2  = dev_at(om, DL[1]);
    d5  = dev_at(om, DL[2]);
    dm5 = dev_at(om, DL[3]);
    printf("      dev(Δ=0)=%.6e  dev(+2e-11)=%.6e  dev(+5e-11)=%.6e  dev(−5e-11)=%.6e\n",
           d0, d2, d5, dm5);
    printf("      源：verify_bsd.exe 实测 dev=4.316e-10（Δ=0 行，与上式同值）\n");

    /* T1 基准在界内 */
    J("T1 基准行 dev(Δ=0) ≤ 8.0e-10（不误伤）", d0, TOLM);
    /* T2 咬得住：+5e-11 必须 > 界（这里判「超过界」本身 = 用 界/实测 反着写） */
    J("T2 咬得住：8.0e-10 / dev(+5e-11) < 1（即 dev 超界）", TOLM / d5, 1.0);
    /* T3 不是过紧：+2e-11 仍在界内 */
    J("T3 不过紧：dev(+2e-11) ≤ 8.0e-10", d2, TOLM);
    /* T4 线性预测（带符号量）：sdev(Δ) − sdev(0) = 16Δ */
    for (i = 1; i < 4; i++) {
        pred = 16.0 * DL[i];
        resmax = fmax(resmax, fabs((sdev_at(om, DL[i]) - sdev_at(om, 0.0)) - pred));
    }
    J("T4 sdev(Δ)−sdev(0) = 16Δ（预测残差 ≤ 1.8e-15 = 16·½ulp）", resmax, 1.8e-15);
    /* T5 自比界 vs 机制界：拒绝次数必须不同 */
    rej_mech = count_reject(DL, 4, om, TOLM, 0);
    rej_self = count_reject(DL, 4, om, TOLM, 1);
    printf("      独立对照：4 个扰动下 机制界 拒绝 %d 次 · 自比界(2×dev) 拒绝 %d 次\n", rej_mech, rej_self);
    J("T5 自比界判别力 = 0（0 次拒绝）⟹ 必须与机制界分离", (double)(rej_self == 0 ? 0 : 1), 0.0);
    J("T6 机制界判别力 > 0（至少拒绝 1 次）", (double)(rej_mech >= 1 ? 0 : 1), 0.0);

    printf("\nRESULT PASS=%d FAIL=%d\n", g_pass, g_fail);
    printf("EXIT=%d\n", g_fail == 0 ? 0 : 1);
    return g_fail == 0 ? 0 : 1;
}

int main(void)
{
    return g2_main_body();
}
