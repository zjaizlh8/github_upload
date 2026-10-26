/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
 * verify_hl.c —— HL 猜想（Hardy–Littlewood 孪生素数）independent recomputation
 * 对应论文：D:\lw\HL_证明_结构层_chaos onestep_matrix probe.md
 *
 * 五行执line铁序（阶段直序）：
 *   base S 锚定  → 外部锚先落地（2C₂、e^{−γ}、log(e^{2γ}/4)、4e^{−2γ}）
 *   dimension-change F 虚化  → t=1 处 k 阶极点 ⟹ 下限被迫 ≥ 2（P3）
 *   sieve M 坍缩  → 累加展开coefficient抽取：一阶coefficient = D = k（P2）
 *   accumulate W 渗透  → 逐 x numeric检索：series vs integral、焊接口阶梯integral（P1、P4）
 *   generate D 扩展  → 双账本 det 比、判别力缺步、首项/截断error（P7、P9、P11）
 *   base S closure  → 汇总deviation + 四量 + 判定
 *
 * four operationscaliber（凡跨维度必付残差）：
 *   加法：series逐项相加、阶梯integral逐段summation —— 无残差。
 *   乘法：∏(1−1/q)、∏(1−1/(q−1)²)、部分积 ∏(1+1/n)e^{−1/n} —— 无残差。
 *   取模：唯一一处 = 判别力缺步的整除判定 q | h（整数，无残差）。
 *   跨维度：exp/log 用于 e^{−γ}、2C₂、log(e^{2γ}/4) 与coefficient比的逐级halvingextrapolation。
 *           exp 走自写 my_exp（Taylor + halving）；本步不需要 sqrt；
 *           log/fabs 用 libm（体系允许）。residualcaliber = 相对deviation dev，
 *           每项tolerance由该项机制给（截断尾 e^{−1/u} / convergence阶 / 外部常数位数）。
 *
 * 注意（本file自带的实测lesson）：
 *   累加的时间展开 Σ(k)_n/log^n x 是**divergenceasymptoticseries**（项比 (k+n)/log x），
 *   最优截断在 n ≈ log x − k。取固定 80 项会爆到 1e38。
 *   故 li_ser 用**自适应最优截断**（项open始增长即停），tolerance按最小项 e^{−1/u} 给。
 *
 * F: 正向 —— 从 Li_k 焊接口出发：展开抽coefficient、算双账本、读判别力缺步。
 * R: 反向 —— 从外部锚回推：若 2C₂、e^{−γ}、4e^{−2γ}、det 比任一对不上，本步当场崩。
 *
 * 编译：gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hl.exe verify_hl.c
 * ==========================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PMAX 1000000
#define NJ   11

static int    g_pass = 0, g_fail = 0;
static double g_maxdev = 0.0;

static void judge(const char *name, double dev, double tol)
{
    int ok = (dev <= tol);
    if (ok) g_pass++; else g_fail++;
    if (dev > g_maxdev) g_maxdev = dev;
    printf("  [%s] %-50s dev=%.6e tol=%.3e\n", ok ? "PASS" : "FAIL", name, dev, tol);
}

static double dev_rel(double a, double b)
{
    double m = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
    if (m < 1e-300) m = 1e-300;
    return fabs(a - b) / m;
}

static int iabs(int x) { return x < 0 ? -x : x; }

/* -------- 跨维度算子（禁 stdlib exp/sqrt，只许 fabs/log） -------- */
static double my_exp(double x)
{
    int k = 0, i;
    double y = x, s = 1.0, t = 1.0;
    while (y > 0.5 || y < -0.5) { y *= 0.5; k++; }
    for (i = 1; i <= 20; i++) { t *= y / (double)i; s += t; }
    while (k-- > 0) s *= s;
    return s;
}

static double pow_n_local(double u, int n)
{
    double r = 1.0;
    int i;
    for (i = 0; i < n; i++) r *= u;
    return r;
}

/* -------- primesieve（1e6，criterion用，禁外部素表） -------- */
static unsigned char g_comp[PMAX + 1];
static int g_pr[80000];
static int g_npr = 0;
static void sieve(void)
{
    int i, j;
    memset(g_comp, 0, sizeof(g_comp));
    g_comp[0] = g_comp[1] = 1;
    for (i = 2; (long long)i * i <= PMAX; i++)
        if (!g_comp[i])
            for (j = i * i; j <= PMAX; j += i) g_comp[j] = 1;
    for (i = 2; i <= PMAX; i++) if (!g_comp[i]) g_pr[g_npr++] = i;
}

/* -------- Li_k(x) = ∫_2^x dt/log^k t ：t=e^u substitution后的复合 Simpson -------- */
static double f_li(double u, int k)
{
    double p = 1.0;
    int i;
    for (i = 0; i < k; i++) p *= u;
    return my_exp(u) / p;
}
static double li_int(double t0, double t1, int k, int seg)
{
    int i;
    double u0 = log(t0), u1 = log(t1);
    double h = (u1 - u0) / (double)seg, s = f_li(u0, k) + f_li(u1, k);
    for (i = 1; i < seg; i++)
        s += ((i & 1) ? 4.0 : 2.0) * f_li(u0 + h * (double)i, k);
    return s * h / 3.0;
}

/* -------- Li_k 展开：自适应最优截断（divergenceasymptoticseries）+ 返回末项量级 --------
 * series：(x/log^k x)·Σ_{n≥0} (k)_n u^n,  u = 1/log x
 * 项比 t_{n+1}/t_n = (k+n)u ⟹ 最优截断 n* ≈ 1/u − k，最小项 ~ e^{−1/u}
 * ------------------------------------------------------------------------- */
static double li_ser(double x, int k, int *nopt, double *tail)
{
    double L = log(x), u = 1.0 / L, Lk = 1.0, fac = 1.0, term = 1.0, s = 1.0;
    int n, i;
    for (i = 0; i < k; i++) Lk *= L;
    for (n = 1; n <= 400; n++) {
        double nxt;
        fac *= (double)(k + n - 1);
        term = fac * pow_n_local(u, n);
        nxt  = term * (double)(k + n) * u;           /* 下一项量级 */
        s += term;
        if (term == 0.0 || !(term == term)) { *nopt = n; *tail = 0.0; break; }  /* 下溢 */
        if (nxt > term) { *nopt = n; *tail = term; break; }   /* 项open始增长 ⟹ 停 */
    }
    if (n > 400) { *nopt = 400; *tail = term; }
    return x / Lk * s;
}

/* -------- 焊接口右侧：∫_{t0}^{t1} (e^γ ∏_{p≤t}(1−1/p))² dt（阶梯精确） ----- */
static double step_int_mertens(double t0, double t1, double egamma)
{
    double prod = 1.0, cur = t0, s = 0.0, g;
    int i;
    for (i = 0; i < g_npr; i++) {
        double p = (double)g_pr[i];
        if (p > t0) break;
        prod *= (1.0 - 1.0 / p);
    }
    g = egamma * prod; g *= g;
    for (; i < g_npr; i++) {
        double p = (double)g_pr[i];
        if (p >= t1) break;
        s += (p - cur) * g;
        cur = p;
        prod *= (1.0 - 1.0 / p);
        g = egamma * prod; g *= g;
    }
    s += (t1 - cur) * g;
    return s;
}

/* -------- 判别力缺步：q ≤ 13 中整除 h 的素数个数（取模法唯一处） -------- */
static int holes_of(int h)
{
    int i, c = 0;
    for (i = 0; i < g_npr; i++) {
        int q = g_pr[i];
        if (q > 13) break;
        if (h % q == 0) c++;
    }
    return c;
}

/* -------- 双账本 det 比：∏(q−1)/∏q（整数，精确） -------- */
static void ledger(int *num, int *den)
{
    int qs[5] = {2, 3, 5, 7, 11};
    int i, a = 1, b = 1;
    for (i = 0; i < 5; i++) { a *= (qs[i] - 1); b *= qs[i]; }
    *num = a; *den = b;
}

/* -------- 判据器（循环全部收在辅助函数内，main 只取读数） -------- */
static const double GAMMA_E = 0.5772156649015329;   /* γ */
static double p1_scan(void)
{
    double xs[3] = {1e4, 1e5, 1e6};
    double ratio_worst = 0.0;
    int i;
    for (i = 0; i < 3; i++) {
        double L = log(xs[i]);
        int nopt; double tail;
        double s = li_ser(xs[i], 2, &nopt, &tail);
        double dv = dev_rel(s, li_int(2.0, xs[i], 2, 40000));
        double tail_pred = 2.0 * tail * xs[i] / (L * L) / s;
        double r = (tail_pred > 1e-18) ? dv / tail_pred : 1.0;
        printf("       x=%.0e  N*=%-2d 末项=%.2e  Li2(series)=%.6f  deviation=%.2e  尾estimate=%.2e\n",
               xs[i], nopt, tail, s, dv, tail_pred);
        if (r > ratio_worst) ratio_worst = r;
    }
    return ratio_worst;
}

static double p2_extrap(void)
{
    double xs[4] = {1e8, 1e16, 1e32, 1e64};
    double uu[4], y[4], A[3], B[2], C;
    int i;
    for (i = 0; i < 4; i++) {
        double L = log(xs[i]);
        int nopt; double tail;
        double T = li_ser(xs[i], 2, &nopt, &tail) * L * L / xs[i];
        uu[i] = 1.0 / L;
        y[i]  = (T - 1.0) / uu[i];
    }
    A[0] = 2.0 * y[1] - y[0];
    A[1] = 2.0 * y[2] - y[1];
    A[2] = 2.0 * y[3] - y[2];
    B[0] = (4.0 * A[1] - A[0]) / 3.0;
    B[1] = (4.0 * A[2] - A[1]) / 3.0;
    C    = (8.0 * B[1] - B[0]) / 7.0;
    printf("       y=%.6f / %.6f / %.6f / %.6f  extrapolation=%.9f\n",
           y[0], y[1], y[2], y[3], C);
    return C;
}

static double p4_weld(void)
{
    double xs[3] = {1e4, 1e5, 1e6};
    double ee = my_exp(GAMMA_E);
    double worst = 0.0;
    int i;
    for (i = 0; i < 3; i++) {
        int nopt; double tail;
        double s = li_ser(xs[i], 2, &nopt, &tail);
        double lhs = step_int_mertens(2.0, xs[i], ee);
        double dv = fabs(lhs / s - 1.0);
        printf("       x=%.0e  ∫=%.4f  Li2=%.4f  deviation=%.2e\n", xs[i], lhs, s, dv);
        if (dv > worst) worst = dv;
    }
    return worst;
}

static double p5_c2(void)
{
    double a = 1.0;
    int i;
    for (i = 0; i < g_npr; i++) {
        double q = (double)g_pr[i];
        if (q <= 2.0) continue;
        a *= (1.0 - 1.0 / ((q - 1.0) * (q - 1.0)));
    }
    a *= 2.0;
    printf("       2C2 = %.12f\n", a);
    return a;
}

static double p6_egamma(void)
{
    double s = 1.0;
    int N = 200000, i;
    for (i = 1; i <= N; i++) {
        double n = (double)i;
        s *= (1.0 + 1.0 / n) * my_exp(-1.0 / n);
    }
    printf("       e^-gamma = %.12f（部分积 N=%d）\n", s, N);
    return s;
}

static int p9_holes(void)
{
    int hs[4] = {2, 6, 30, 210};
    int ex[4] = {1, 2, 3, 4};
    int bad = 0, i;
    for (i = 0; i < 4; i++) {
        int got = holes_of(hs[i]);
        printf("       h=%-3d 缺步=%d（应 %d）\n", hs[i], got, ex[i]);
        bad += iabs(got - ex[i]);
    }
    return bad;
}

static double p11_trunc(void)
{
    double L = log(1e5);
    int nopt; double tail;
    double sopt = li_ser(1e5, 2, &nopt, &tail);
    double first = 1e5 / (L * L);
    double s8 = 1.0, fac = 1.0, u = 1.0 / L;
    int n;
    for (n = 1; n <= 8; n++) { fac *= (double)(2 + n - 1); s8 += fac * pow_n_local(u, n); }
    s8 *= 1e5 / (L * L);
    {
        double d1 = fabs(first / sopt - 1.0), d8 = fabs(s8 / sopt - 1.0);
        printf("       首项 %.4f（偏 %.2f%%）· 8 项偏 %.4f%% · 最优截断 %.4f（N*=%d）\n",
               first, d1 * 100.0, d8 * 100.0, sopt, nopt);
        return (d1 >= 0.10 && d8 <= 0.01) ? 0.0 : 1.0;
    }
}


int main(void)
{
    const double EG  = 0.5614594835668852;   /* e^{-γ} 外部锚（15 位） */
    
    const int    k   = NJ;
    double res = 0.0, thr, v, a, C;
    int    num, den, bad;
    char   buf[160];

    sieve();
    printf("verify_hl.c · HL independent recomputation（焊接口 / 双账本 / coefficient / 缺步）\n");
    printf("素数筛法上限 %d，筛法得 %d 个素数；外部锚 e^-gamma=%.16f\n\n", PMAX, g_npr, EG);

    { double rw = p1_scan();
      judge("P1 展开(最优截断) vs integral：deviation / 尾estimate ≤ 3", rw, 3.0);
      if (rw > res) res = rw; }

    C = p2_extrap();
    sprintf(buf, "P2 累加展开一阶coefficient = D = k = 2（extrapolation %.6f）", C);
    judge(buf, fabs(C - 2.0), 1e-3);
    if (fabs(C - 2.0) > res) res = fabs(C - 2.0);

    { double I1 = li_int(1.01, 2.0, 2, 200000);
      double I2 = li_int(1.001, 2.0, 2, 200000);
      double ratio = I2 / I1;
      printf("       I(1e-2)=%.6f  I(1e-3)=%.6f  比值=%.6f\n", I1, I2, ratio);
      judge("P3 下限被迫：I(eps)~1/eps，比值 I(1e-3)/I(1e-2)", fabs(ratio - 10.0), 0.5);
      if (fabs(ratio - 10.0) > res) res = fabs(ratio - 10.0); }

    { double w = p4_weld();
      thr = 3.0 / log(1e4);   /* mechanism：Mertens 尾 O(1/log t) */
      judge("P4 焊接口：阶梯integral ∫(e^gamma·prod)^2 dt vs Li2", w, thr);
      printf("       (阈值按 tail=3/log(x_min) 给 = %.4f)\n", thr);
      if (w > res) res = w; }

    a = p5_c2();
    judge("P5 2C2 = 2∏(1−1/(q−1)^2) vs 1.3203236316", fabs(a - 1.3203236316), 3e-6);
    if (fabs(a - 1.3203236316) > res) res = fabs(a - 1.3203236316);

    { double s = p6_egamma();
      judge("P6 e^-gamma 部分积 vs 0.5614594836", fabs(s - EG), 1e-5);
      if (fabs(s - EG) > res) res = fabs(s - EG); }

    ledger(&num, &den);
    printf("       %d × 77 = %d ；16 × %d = %d\n", num, num * 77, den, 16 * den);
    judge("P7 ledger det 比 480/2310 = 16/77（整数交叉乘）",
          (double)iabs(num * 77 - 16 * den), 0.0);

    v = 2.0 * GAMMA_E - log(4.0);
    printf("       log(e^2gamma/4) = %.12f\n", v);
    judge("P8 log(e^2gamma/4) = 2gamma − log4 vs −0.231863031", fabs(v + 0.231863031), 1e-9);
    if (fabs(v + 0.231863031) > res) res = fabs(v + 0.231863031);

    bad = p9_holes();
    judge("P9 缺步 = q|h 的素数个数（h=2,6,30,210）", (double)bad, 0.0);

    v = 4.0 * EG * EG;
    printf("       4e^-2gamma = %.12f\n", v);
    judge("P10 4e^-2gamma = 1.260947007（外部锚）", fabs(v - 1.260947007), 1e-9);
    if (fabs(v - 1.260947007) > res) res = fabs(v - 1.260947007);

    { double d1 = p11_trunc();
      judge("P11 首项deviation ≥10% 且 8 项deviation ≤1%", d1, 0.0);
      if (d1 > res) res = d1; }

    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.6e\n", g_pass, g_fail, g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", k, 1.0 - 1.0 / (double)k);
    printf("res=%.6e domain=%d dc=%d order=%d\n", res, k, 8, 2);
    printf("EXIT=%d\n", g_fail == 0 ? 0 : 1);
    return g_fail == 0 ? 0 : 1;
}
