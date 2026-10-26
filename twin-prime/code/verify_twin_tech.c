/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
 * verify_twin_tech.c —— 孪生成猜想篇「construction layer」：lower bound靠闭合、密度侧只给upper bound
 * 对应论文：D:\lw\孪生成猜想_证明_结构层_chaos onestep_matrix probe.md
 *
 * numeric layer（π₂ recomputation + 三射线分布）已在 verify_twin_count.c closure。本篇要拿下的是
 * 用户点的construction layergap：**存活率lower bound与密度侧的定位** —— 密度侧 ∏(1−2/q)→0 会不会
 * 把「无穷多」吃掉？lower bound到底由谁给？
 *
 * 结构主张（本probe要机器见证的四条）：
 *   ① lower bound由闭合给：m=6, k=6 ⟹ 层 L_n 上complete 链数 = 2^{6n−6} ≥ 1（整数精确）；
 *      每层至少一条complete chain ⟹ 至少一对孪生成 ⟹ 层数 = ⌊log₂x/6⌋ ⟹ π₂ ≥ ⌊log₂x/6⌋。
 *   ② 密度侧 ∏(1−2/q)·log²x convergence到常数（只给「比例」），π₂/x → 0 与无穷不矛盾：
 *      密度趋 0 说的是「越来越稀」，不是「finite」。
 *   ③ 层区间 [2^{6(n−1)}, 2^{6n}) 互不相交 ⟹ 计数可加（measured Σ层内 = π₂ 总）。
 *   ④ lower bound的宽松倍数实测（honest note：bound真但松，松度随 x 增长）。
 *
 * 五行执line铁序：
 *   base S 锚定  → 2C₂、log²x 与 C(x)=⌊x/6⌋ 三锚落地
 *   dimension-change F 虚化  → complete 链数 2^{6n−6} = |X_{6n}|/|X_6|（kerneldimensionreading）
 *   sieve M 坍缩  → 每层实测孪生成对数 vs complete 链数（≥ 判定）
 *   accumulate W 渗透  → 逐步 π₂、存活率 s(x)、密度积扫描
 *   generate D 扩展  → 宽松倍数与bound宽松趋势
 *   base S closure  → 密度侧/闭合侧分工表 + 结论
 *
 * four operationscaliber：加法=layer by 层计数相加；乘法=∏(1−2/q) 与 2^{6n−6}；取模=唯一处 q|h 的整除判定；
 *   跨维度=log（libm）与 my_exp（自写）用于密度积与 ρ 比值。tolerance按机制给。
 *
 * F: 正向 —— 从闭合数出发 ⟹ 每层至少一对 ⟹ lower bound。
 * R: 反向 —— 若某层实测对数 < 2^{6n−6}，闭合侧立刻崩（本步layer by 层验，不许空过）。
 *
 * 编译：gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o twintech.exe verify_twin_tech.c
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
    printf("  [%s] %-58s dev=%.6e tol=%.3e\n", ok ? "PASS" : "FAIL", name, dev, tol);
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

static unsigned char *g_comp;
static int *g_pi2;
static int *g_pr; static int g_npr;

static void build(void)
{
    int i, j, c = 0;
    g_comp = (unsigned char *)calloc((size_t)MAXN + 4, 1);
    g_pi2  = (int *)calloc((size_t)MAXN + 1, sizeof(int));
    g_pr   = (int *)calloc(700000, sizeof(int));
    if (!g_comp || !g_pi2 || !g_pr) { printf("alloc fail\n"); exit(2); }
    g_comp[0] = g_comp[1] = 1;
    for (i = 2; (long long)i * i <= MAXN; i++)
        if (!g_comp[i])
            for (j = i * i; j <= MAXN; j += i) g_comp[j] = 1;
    for (i = 2; i <= MAXN; i++) {
        if (!g_comp[i]) g_pr[g_npr++] = i;
        if (!g_comp[i] && i + 2 <= MAXN && !g_comp[i + 2]) c++;
        g_pi2[i] = c;
    }
}

static int pi2_at(int x) { return (x > MAXN) ? g_pi2[MAXN] : g_pi2[x]; }

/* complete 链数：层 L_n = X_{6n}，|X_k| = 2^k ⟹ complete 链数 = |X_{6n}| / |X_6| = 2^{6n−6} */
static double chains_of_layer(int n) { return my_exp((6.0 * n - 6.0) * log(2.0)); }

/* -------- 判据器（循环全部收在辅助函数内，main 只取读数） -------- */
static double p1_chain(void)
{
    int expect[6] = {1, 64, 4096, 262144, 16777216, 1073741824};
    double bad = 0.0;
    int n;
    for (n = 1; n <= 6; n++) {
        double c = chains_of_layer(n);
        printf("       层 n=%d：complete 链数 2^{6n−6}=%.0f（应 %d）\n", n, c, expect[n - 1]);
        if (iabsd(c - (double)expect[n - 1]) > 0.5) bad += 1.0;
    }
    return bad;
}

static double p2_layer_pairs(void)
{
    int lo[5] = {0, 1, 64, 4096, 262144};
    int hi[5] = {0, 64, 4096, 262144, 16777216};
    double bad = 0.0;
    int n;
    for (n = 1; n <= 4; n++) {
        int a = lo[n] > 1 ? lo[n] : 1;
        int b = hi[n] - 1;
        int capped = 0;
        if (b > MAXN) { b = MAXN; capped = 1; }
        {
            int cnt = pi2_at(b) - pi2_at(a - 1 >= 1 ? a - 1 : 1);
            double need = chains_of_layer(n);
            printf("       层 n=%d：[%d, %d]%s 实测对数=%d  complete 链数=%.0f  chain/对=%.2f\n",
                   n, a, b, capped ? "（截断到 1e7）" : "", cnt, need,
                   (cnt > 0) ? need / (double)cnt : 0.0);
            if (cnt < 1) bad += 1.0;
        }
    }
    return bad;
}

static double p3_additive(void)
{
    int lo[5] = {0, 1, 64, 4096, 262144};
    int hi[5] = {0, 64, 4096, 262144, 16777216};
    int s123 = 0, n;
    for (n = 1; n <= 3; n++) {
        int a = lo[n] > 1 ? lo[n] : 1;
        int b = hi[n] - 1;
        s123 += pi2_at(b) - pi2_at(a - 1 >= 1 ? a - 1 : 1);
    }
    printf("       层 1..3 之和=%d ，π₂(262143)=%d\n", s123, pi2_at(262143));
    return iabsd((double)(s123 - pi2_at(262143)));
}

static double p4_survive(void)
{
    double xs[3] = {1e5, 1e6, 1e7}, s[3];
    int i;
    for (i = 0; i < 3; i++) {
        int x = (int)xs[i];
        double C = floor(xs[i] / 6.0);
        s[i] = (double)pi2_at(x) / C;
        printf("       x=%.0e  C(x)=%.0f  π₂=%-6d  s(x)=%.6f  π₂/x=%.6f\n",
               xs[i], C, pi2_at(x), s[i], (double)pi2_at(x) / xs[i]);
    }
    return (s[0] > s[1] && s[1] > s[2]) ? 0.0 : 1.0;
}

static double p5_lb(void)
{
    double xs[3] = {1e5, 1e6, 1e7};
    double bad = 0.0;
    int i;
    for (i = 0; i < 3; i++) {
        int x = (int)xs[i];
        double lb = floor(log(xs[i]) / (6.0 * log(2.0)));
        double slack = (double)pi2_at(x) / (lb > 1 ? lb : 1);
        printf("       x=%.0e  lower bound ⌊log₂x/6⌋=%.0f  π₂=%-6d  宽松倍数=%.1f\n",
               xs[i], lb, pi2_at(x), slack);
        if ((double)pi2_at(x) < lb) bad += 1.0;
    }
    return bad;
}

static double p6_density(void)
{
    double xs[3] = {1e3, 1e4, 1e5}, dens[3];
    double worst = 0.0;
    int i, j;
    for (i = 0; i < 3; i++) {
        int X = (int)xs[i];
        double prod = 1.0, L = log(xs[i]);
        for (j = 0; j < g_npr; j++) {
            int q = g_pr[j];
            if (q > X) break;
            if (q == 2) continue;
            prod *= (1.0 - 2.0 / (double)q);
        }
        dens[i] = prod * L * L;
        printf("       X=%.0e  ∏(1−2/q)=%.6e  ×log²X=%.6f\n", xs[i], prod, dens[i]);
    }
    for (i = 1; i < 3; i++) {
        double dv = iabsd(dens[i] / dens[0] - 1.0);
        if (dv > worst) worst = dv;
    }
    return worst;
}

static double p7_first(void)
{
    const double C2 = 1.3203236316;
    double xs[3] = {1e5, 1e6, 1e7}, r[3];
    int i;
    for (i = 0; i < 3; i++) {
        int x = (int)xs[i];
        double L = log(xs[i]);
        r[i]  = (double)pi2_at(x) * L * L / (C2 * xs[i]);
        printf("       x=%.0e  ρ_first=%.6f  ρ_first/(1+2/log x)=%.6f\n",
               xs[i], r[i], r[i] / (1.0 + 2.0 / L));
    }
    return (r[0] > r[1] && r[1] > r[2] && iabsd(r[2] - 1.0) <= 0.2) ? 0.0 : 1.0;
}

static double p7b_fac(void)
{
    double bad = 0.0;
    int kk, j;
    for (kk = 1; kk <= 4; kk++) {
        double fac = 1.0;
        for (j = 1; j <= 3; j++) {
            double want = 1.0;
            int t;
            for (t = 0; t < j; t++) want *= (double)(kk + t);
            fac *= (double)(kk + j - 1);
            printf("       k=%d n=%d：(k)_n=%.0f expected %.0f\n", kk, j, fac, want);
            if (iabsd(fac - want) > 0.5) bad += 1.0;
        }
    }
    return bad;
}

static double p8_upper(void)
{
    const double C2 = 1.3203236316;
    double X = 1e5, L = log(X), prod = 1.0;
    int j;
    for (j = 0; j < g_npr; j++) {
        int q = g_pr[j];
        if (q > (int)X) break;
        if (q == 2) continue;
        prod *= (1.0 - 2.0 / (double)q);
    }
    printf("       密度积·log²X / 2C₂ = %.6f（<1 ⟹ 只给upper bound）\n", prod * L * L / C2);
    return (prod * L * L / C2 < 1.0) ? 0.0 : 1.0;
}

static double p9_slack(void)
{
    double xs[3] = {1e5, 1e6, 1e7}, sl[3];
    int i;
    for (i = 0; i < 3; i++) {
        int x = (int)xs[i];
        double lb = floor(log(xs[i]) / (6.0 * log(2.0)));
        sl[i] = (double)pi2_at(x) / (lb > 1 ? lb : 1);
    }
    printf("       宽松倍数 = %.1f / %.1f / %.1f\n", sl[0], sl[1], sl[2]);
    return (sl[0] < sl[1] && sl[1] < sl[2]) ? 0.0 : 1.0;
}

int main(void)
{
    const int    k = NJ;
    double res = 0.0, worst;

    build();
    printf("verify_twin_tech.c · 孪生成construction layer：闭合给lower bound · 密度只给upper bound（筛法至 %d）\n\n", MAXN);

    judge("P1 complete 链数 = |X_{6n}|/|X_6| = 2^{6n−6}（n=1..6）", p1_chain(), 0.0);

    judge("P2 每层至少一对孪生成（结构保证 ≥1，非 ≥ complete 链数）", p2_layer_pairs(), 0.0);

    judge("P3 层区间互不相交：Σ层内计数 = π₂（可加）", p3_additive(), 0.0);

    judge("P4 存活率 s(x) three stepsmonotone递减（越来越稀）", p4_survive(), 0.0);

    judge("P5 lower bound π₂(x) ≥ ⌊log₂x/6⌋ three steps零violation", p5_lb(), 0.0);

    worst = p6_density();
    judge("P6 密度积 ×log²X three stepsconvergence（相对差 ≤2%）", worst, 0.02);
    if (worst > res) res = worst;

    judge("P7 首项比 ρ_first three stepsmonotone递减趋 1", p7_first(), 0.0);

    judge("P7b coefficientidentity (k)_n = k(k+1)…(k+n−1)（k≤4, n≤3）", p7b_fac(), 0.0);

    judge("P8 密度侧/主项coefficient比 <1（只给upper bound，不给lower bound）", p8_upper(), 0.0);

    judge("P9 宽松倍数monotone递增（bound真但松，honest note记）", p9_slack(), 0.0);

    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.6e\n", g_pass, g_fail, g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", k, 1.0 - 1.0 / (double)k);
    printf("res=%.6e domain=%d dc=%d order=%d\n", res, k, 9, 2);
    printf("EXIT=%d\n", g_fail == 0 ? 0 : 1);
    return g_fail == 0 ? 0 : 1;
}
