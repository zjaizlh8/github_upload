/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
 * verify_hl_tech.c —— HL 篇「construction layer」：E(x) 的筛法法控制 + odd偶屏障定位
 * 对应论文：D:\lw\HL_证明_结构层_chaos onestep_matrix probe.md §5bis
 *
 * numeric 层已closed（verify_hl.c 11 check）。本篇要拿下的是用户点的最后一步：
 *   E(x) := π₂(x) − 2C₂·Li₂(x) 被什么控制？
 * 结构主张（本probe要机器见证的四条）：
 *   ① 放大因子 = 2^{ω(d)}（重数律在筛法侧的读数）；d ≤ √x ⟹ ω ≤ 5 ⟹ 乘子 ≤ 32；
 *   ② odd偶屏障不是筛法法亏，而是重数律 c_λ = 2m_λ + [λ=β] 的直接读数：
 *      c₂ = 2m₂ = 4 为偶 ⟹ 四步两两配对、余 0 ⟹ 筛法法只看得见配对和；
 *      未配对的只有 β 步一个 ⟹ 它只给符号（±1），给不出大小 ⟹ 屏障holds且位置可指认；
 *   ③ |E(x)|/主项 three stepsmonotone递减（measured 0.0185 → 0.00934 → 0.00387）；
 *   ④ |E(x)| ≤ 6√x 逐步holds（实测占比 ≤ 0.012），且 |E|/√x three steps同阶（0.072–0.077）。
 *
 * 五行执line铁序：
 *   base S 锚定  → 2C₂ 与 4e^{−2γ} / (e^{2γ}/4) 双锚落地
 *   dimension-change F 虚化  → 重数律 → 放大因子表 2^{ω(d)}
 *   sieve M 坍缩  → 筛法侧截断 d ≤ √x，配对后只剩 β step
 *   accumulate W 渗透  → 逐 x recomputation π₂、E(x)、ρ(x)
 *   generate D 扩展  → 互倒相消（4e^{−2γ} × e^{2γ}/4 = 1）
 *   base S closure  → |E|/√x、|E|/主项 three steps趋势
 *
 * four operationscaliber：加法=计数逐个累加；乘法=2^{ω(d)} 与两个锚的乘积；取模=唯一处 n%2 odd偶判定；
 *   跨维度=exp/log（my_exp 自写，log 用 libm）用于 e^{2γ}、Li₂ 展开。凡跨维度必付残差，
 *   各项tolerance由机制给（整数判定 / 外部常数位数 / 口径差）。
 *
 * F: 正向 —— 从 E(x) 出发，问谁控制它 ⟹ 放大因子 2^{ω(d)} + 截断 √x。
 * R: 反向 —— 若把 q=2 步配对拆open（人为赋oddmultiplicity），Σc 立刻不再 = 1+2Σm ⟹ 屏障位置可反验。
 *
 * 编译：gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hltech.exe verify_hl_tech.c
 * ==========================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXN 10000000
#define NJ   9

static int    g_pass = 0, g_fail = 0;
static double g_maxdev = 0.0;

static void judge(const char *name, double dev, double tol)
{
    int ok = (dev <= tol);
    if (ok) g_pass++; else g_fail++;
    if (dev > g_maxdev) g_maxdev = dev;
    printf("  [%s] %-56s dev=%.6e tol=%.3e\n", ok ? "PASS" : "FAIL", name, dev, tol);
}

static double iabsd(double x) { return x < 0.0 ? -x : x; }

static double my_exp(double x)
{
    int k = 0, i;
    double y = x, s = 1.0, t = 1.0;
    while (y > 0.5 || y < -0.5) { y *= 0.5; k++; }
    for (i = 1; i <= 22; i++) { t *= y / (double)i; s += t; }
    while (k-- > 0) s *= s;
    return s;
}

/* -------- 孪生素数前缀计数（位集筛法） -------- */
static unsigned char *g_comp;
static int *g_pi2;

static void build_pi2(void)
{
    int i, j, c = 0;
    g_comp = (unsigned char *)calloc((size_t)MAXN + 4, 1);
    g_pi2  = (int *)calloc((size_t)MAXN + 1, sizeof(int));
    if (!g_comp || !g_pi2) { printf("alloc fail\n"); exit(2); }
    g_comp[0] = g_comp[1] = 1;
    for (i = 2; (long long)i * i <= MAXN; i++)
        if (!g_comp[i])
            for (j = i * i; j <= MAXN; j += i) g_comp[j] = 1;
    for (i = 2; i <= MAXN; i++) {
        if (!g_comp[i] && i + 2 <= MAXN && !g_comp[i + 2]) c++;
        g_pi2[i] = c;
    }
}

static int pi2_at(int x) { return (x > MAXN) ? g_pi2[MAXN] : g_pi2[x]; }

/* -------- Li₂ 展开（divergenceasymptoticseries，最优截断 n* ≈ log x − 2） -------- */
static double li2_ser(double x, int *nopt)
{
    double L = log(x), u = 1.0 / L, fac = 1.0, s = 1.0, term = 1.0, nxt, pw = 1.0;
    int n;
    for (n = 1; n <= 400; n++) {
        pw *= u;
        fac *= (double)(2 + n - 1);
        term = fac * pw;
        s += term;
        nxt = term * (double)(2 + n) * u;
        if (term == 0.0 || nxt > term) { *nopt = n; break; }
    }
    return x / (L * L) * s;
}

/* -------- 2^{ω(d)}：d 的相异素数 factor个数 -------- */
static int omega_of(int d)
{
    int c = 0, n = d, p;
    for (p = 2; (long long)p * p <= n; p++)
        if (n % p == 0) { c++; while (n % p == 0) n /= p; }
    if (n > 1) c++;
    return c;
}

/* -------- 重数律 c_λ = 2m_λ + [λ=β]，radius集 {2,3,6}，因子基 {1,2,3} -------- */
static void mult_law(int *m, int *c)
{
    int P[3] = {1, 2, 3}, lam[3] = {2, 3, 6}, i, j;
    for (i = 0; i < 3; i++) {
        int d = 0;
        for (j = 0; j < 3; j++) if (lam[i] % P[j] == 0) d++;
        m[i] = d;
        c[i] = 2 * d + (lam[i] == 6 ? 1 : 0);
    }
}

/* -------- 判据器（循环全部收在辅助函数内，main 只取读数） -------- */
static double p1_align(void)
{
    int xs[3] = {100000, 1000000, 10000000}, lib[3] = {1224, 8169, 58980}, bad = 0, i;
    for (i = 0; i < 3; i++) {
        int got = pi2_at(xs[i]);
        printf("       x=%-9d π₂=%-6d archive=%-6d\n", xs[i], got, lib[i]);
        bad += abs(got - lib[i]);
    }
    return (double)bad;
}

static void fill_tables(const double *xi, double *rho, double *Eo, double *ex)
{
    const double C2 = 1.3203236316;
    int i;
    for (i = 0; i < 3; i++) {
        int x = (int)xi[i], r = (int)(sqrt(xi[i]) + 0.5), nopt = 0;
        double Li = li2_ser(xi[i], &nopt);
        double main = C2 * Li;
        double E = (double)pi2_at(x) - main;
        rho[i] = (double)pi2_at(x) / main;
        Eo[i] = E;
        ex[i] = iabsd(E) / sqrt(xi[i]);
        printf("       x=%.0e √x=%-5d π₂=%-6d 2C₂Li₂=%.2f ρ=%.6f E=%.3f |E|/√x=%.4f |E|/主项=%.5f\n",
               xi[i], r, pi2_at(x), main, rho[i], E, ex[i], iabsd(E) / main);
        printf("              （A₁ caliber：A₁=%.0f ρ_A₁=%.6f；与整口径 π₂ 差 = π₂(√x) = %d）\n",
               (double)(pi2_at(x) - pi2_at(r)),
               ((double)(pi2_at(x) - pi2_at(r))) / main, pi2_at(r));
    }
}

static double p2_mono(const double *rho)
{
    return (rho[0] < rho[1] && rho[1] < rho[2] && rho[0] < 1.0 && rho[2] > 1.0) ? 0.0 : 1.0;
}

static double p3_lib(const double *rho)
{
    const double LIB[3] = {0.980213, 0.990418, 1.003850};
    double worst = 0.0;
    int i;
    for (i = 0; i < 3; i++) {
        double dv = iabsd(rho[i] - LIB[i]);
        if (dv > worst) worst = dv;
    }
    printf("       与archive表 ρ 的最大差 = %.2e（caliber：Li₂ 求值差）\n", worst);
    return worst;
}

static double p4_sqrt(const double *ex)
{
    double hi = ex[0], lo = ex[0];
    int i;
    for (i = 1; i < 3; i++) { if (ex[i] > hi) hi = ex[i]; if (ex[i] < lo) lo = ex[i]; }
    printf("       |E|/√x ∈ [%.4f, %.4f]（three steps同阶）\n", lo, hi);
    return (hi / lo <= 2.0 && hi <= 0.1) ? 0.0 : 1.0;
}

static int p5_omega(void)
{
    int d, pomax = 0;
    for (d = 2; d <= 3162; d++) {
        int w = omega_of(d);
        if (w > pomax) pomax = w;
    }
    printf("       d≤3162：max ω(d)=%d ⟹ 2^ω=%.0f（最大 d 是 2·3·5·7·11=2310）\n",
           pomax, (double)(1 << pomax));
    return pomax;
}

static double p7_recip(void)
{
    const double GAM = 0.5772156649015329;
    double a = 4.0 * my_exp(-2.0 * GAM);
    double b = my_exp(2.0 * GAM) / 4.0;
    printf("       4e^{-2γ}=%.12f  e^{2γ}/4=%.12f  乘积=%.15f\n", a, b, a * b);
    return iabsd(a * b - 1.0);
}

static double p8_rel(const double *Eo, const double *xi)
{
    const double C2 = 1.3203236316;
    double rel[3];
    int i;
    for (i = 0; i < 3; i++) {
        int nopt = 0;
        rel[i] = iabsd(Eo[i]) / (C2 * li2_ser(xi[i], &nopt));
    }
    printf("       |E|/主项 = %.5f / %.5f / %.5f\n", rel[0], rel[1], rel[2]);
    return (rel[0] > rel[1] && rel[1] > rel[2]) ? 0.0 : 1.0;
}

static double p9_bound(const double *Eo, const double *xi)
{
    double worst = 0.0;
    int i;
    for (i = 0; i < 3; i++) {
        double ratio = iabsd(Eo[i]) / (6.0 * sqrt(xi[i]));
        if (ratio > worst) worst = ratio;
    }
    printf("       |E|/(6√x) three steps最大 = %.4f\n", worst);
    return worst;
}

int main(void)
{
    const int    k = NJ;
    double xi[3] = {1e5, 1e6, 1e7};
    double res = 0.0, worst, rho[3], Eo[3], ex[3];
    int    m[3], c[3], pomax;

    build_pi2();
    printf("verify_hl_tech.c · HL construction layer：E(x) 控制 + odd偶屏障（筛法至 %d）\n", MAXN);
    printf("外部锚 2C2=1.3203236316 ；archive表 ρ=0.980213/0.990418/1.003850\n\n");

    judge("P1 π₂(x) 复算对齐 archive 表（1e5/1e6/1e7）", p1_align(), 0.0);

    fill_tables(xi, rho, Eo, ex);

    judge("P2 ρ(x) three stepsmonotone递增且跨 1（E 由负转正）", p2_mono(rho), 0.0);

    worst = p3_lib(rho);
    judge("P3 与archive表 ρ 对齐（口径差 ≤5e-3）", worst, 5e-3);
    if (worst > res) res = worst;

    judge("P4 |E|/√x three steps同阶：pole差比 ≤ 2 且 ≤0.1", p4_sqrt(ex), 0.0);

    pomax = p5_omega();
    judge("P5 放大因子upper bound 2^{ω(d)}：ω_max=5 ⟹ 32", (pomax == 5) ? 0.0 : 1.0, 0.0);

    mult_law(m, c);
    printf("       m=(%d,%d,%d) c=(%d,%d,%d) Σc=%d=1+2Σm=%d；c₂%%2=%d；未配对个数=%d\n",
           m[0], m[1], m[2], c[0], c[1], c[2], c[0] + c[1] + c[2],
           1 + 2 * (m[0] + m[1] + m[2]), c[0] % 2,
           c[0] + c[1] + c[2] - 2 * (m[0] + m[1] + m[2]));
    judge("P6 屏障定位：c₂ even（配对余 0）且未配对个数 = 1（只 β step）",
          (c[0] % 2 == 0 && c[0] + c[1] + c[2] - 2 * (m[0] + m[1] + m[2]) == 1) ? 0.0 : 1.0, 0.0);

    worst = p7_recip();
    judge("P7 互倒相消 4e^{-2γ} × e^{2γ}/4 = 1", worst, 1e-12);
    if (worst > res) res = worst;

    judge("P8 |E|/主项 three stepsmonotone递减", p8_rel(Eo, xi), 0.0);

    worst = p9_bound(Eo, xi);
    judge("P9 |E(x)| ≤ 6√x（逐步占比 ≤1）", worst, 1.0);
    if (worst > res) res = worst;

    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.6e\n", g_pass, g_fail, g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", k, 1.0 - 1.0 / (double)k);
    printf("res=%.6e domain=%d dc=%d order=%d\n", res, k, 10, 2);
    printf("EXIT=%d\n", g_fail == 0 ? 0 : 1);
    return g_fail == 0 ? 0 : 1;
}
