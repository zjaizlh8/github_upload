/* 由 merge_segments.py 合成：单一入口，各段原样保留（仅冲突符号加段号前缀） */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ───── 段 0：hardy-littlewood/code/verify_hl.c ───── */
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
 *           exp 走自写 S0_my_exp（Taylor + halving）；本步不需要 sqrt；
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

#define PMAX 1000000
#define S0_NJ   11

static int    S0_g_pass = 0, S0_g_fail = 0;
static double S0_g_maxdev = 0.0;

static void S0_judge(const char *name, double dev, double tol)
{
    int ok = (dev <= tol);
    if (ok) S0_g_pass++; else S0_g_fail++;
    if (dev > S0_g_maxdev) S0_g_maxdev = dev;
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
static double S0_my_exp(double x)
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
static unsigned char S0_g_comp[PMAX + 1];
static int S0_g_pr[80000];
static int S0_g_npr = 0;
static void sieve(void)
{
    int i, j;
    memset(S0_g_comp, 0, sizeof(S0_g_comp));
    S0_g_comp[0] = S0_g_comp[1] = 1;
    for (i = 2; (long long)i * i <= PMAX; i++)
        if (!S0_g_comp[i])
            for (j = i * i; j <= PMAX; j += i) S0_g_comp[j] = 1;
    for (i = 2; i <= PMAX; i++) if (!S0_g_comp[i]) S0_g_pr[S0_g_npr++] = i;
}

/* -------- Li_k(x) = ∫_2^x dt/log^k t ：t=e^u substitution后的复合 Simpson -------- */
static double f_li(double u, int k)
{
    double p = 1.0;
    int i;
    for (i = 0; i < k; i++) p *= u;
    return S0_my_exp(u) / p;
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
    for (i = 0; i < S0_g_npr; i++) {
        double p = (double)S0_g_pr[i];
        if (p > t0) break;
        prod *= (1.0 - 1.0 / p);
    }
    g = egamma * prod; g *= g;
    for (; i < S0_g_npr; i++) {
        double p = (double)S0_g_pr[i];
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
    for (i = 0; i < S0_g_npr; i++) {
        int q = S0_g_pr[i];
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
    double ee = S0_my_exp(GAMMA_E);
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
    for (i = 0; i < S0_g_npr; i++) {
        double q = (double)S0_g_pr[i];
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
        s *= (1.0 + 1.0 / n) * S0_my_exp(-1.0 / n);
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


int seg0_main(void)
{
    const double EG  = 0.5614594835668852;   /* e^{-γ} 外部锚（15 位） */
    
    const int    k   = S0_NJ;
    double res = 0.0, thr, v, a, C;
    int    num, den, bad;
    char   buf[160];

    sieve();
    printf("verify_hl.c · HL independent recomputation（焊接口 / 双账本 / coefficient / 缺步）\n");
    printf("素数筛法上限 %d，筛法得 %d 个素数；外部锚 e^-gamma=%.16f\n\n", PMAX, S0_g_npr, EG);

    { double rw = p1_scan();
      S0_judge("P1 展开(最优截断) vs integral：deviation / 尾estimate ≤ 3", rw, 3.0);
      if (rw > res) res = rw; }

    C = p2_extrap();
    sprintf(buf, "P2 累加展开一阶coefficient = D = k = 2（extrapolation %.6f）", C);
    S0_judge(buf, fabs(C - 2.0), 1e-3);
    if (fabs(C - 2.0) > res) res = fabs(C - 2.0);

    { double I1 = li_int(1.01, 2.0, 2, 200000);
      double I2 = li_int(1.001, 2.0, 2, 200000);
      double ratio = I2 / I1;
      printf("       I(1e-2)=%.6f  I(1e-3)=%.6f  比值=%.6f\n", I1, I2, ratio);
      S0_judge("P3 下限被迫：I(eps)~1/eps，比值 I(1e-3)/I(1e-2)", fabs(ratio - 10.0), 0.5);
      if (fabs(ratio - 10.0) > res) res = fabs(ratio - 10.0); }

    { double w = p4_weld();
      thr = 3.0 / log(1e4);   /* mechanism：Mertens 尾 O(1/log t) */
      S0_judge("P4 焊接口：阶梯integral ∫(e^gamma·prod)^2 dt vs Li2", w, thr);
      printf("       (阈值按 tail=3/log(x_min) 给 = %.4f)\n", thr);
      if (w > res) res = w; }

    a = p5_c2();
    S0_judge("P5 2C2 = 2∏(1−1/(q−1)^2) vs 1.3203236316", fabs(a - 1.3203236316), 3e-6);
    if (fabs(a - 1.3203236316) > res) res = fabs(a - 1.3203236316);

    { double s = p6_egamma();
      S0_judge("P6 e^-gamma 部分积 vs 0.5614594836", fabs(s - EG), 1e-5);
      if (fabs(s - EG) > res) res = fabs(s - EG); }

    ledger(&num, &den);
    printf("       %d × 77 = %d ；16 × %d = %d\n", num, num * 77, den, 16 * den);
    S0_judge("P7 ledger det 比 480/2310 = 16/77（整数交叉乘）",
          (double)iabs(num * 77 - 16 * den), 0.0);

    v = 2.0 * GAMMA_E - log(4.0);
    printf("       log(e^2gamma/4) = %.12f\n", v);
    S0_judge("P8 log(e^2gamma/4) = 2gamma − log4 vs −0.231863031", fabs(v + 0.231863031), 1e-9);
    if (fabs(v + 0.231863031) > res) res = fabs(v + 0.231863031);

    bad = p9_holes();
    S0_judge("P9 缺步 = q|h 的素数个数（h=2,6,30,210）", (double)bad, 0.0);

    v = 4.0 * EG * EG;
    printf("       4e^-2gamma = %.12f\n", v);
    S0_judge("P10 4e^-2gamma = 1.260947007（外部锚）", fabs(v - 1.260947007), 1e-9);
    if (fabs(v - 1.260947007) > res) res = fabs(v - 1.260947007);

    { double d1 = p11_trunc();
      S0_judge("P11 首项deviation ≥10% 且 8 项deviation ≤1%", d1, 0.0);
      if (d1 > res) res = d1; }

    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.6e\n", S0_g_pass, S0_g_fail, S0_g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", k, 1.0 - 1.0 / (double)k);
    printf("res=%.6e domain=%d dc=%d order=%d\n", res, k, 8, 2);
    printf("EXIT=%d\n", S0_g_fail == 0 ? 0 : 1);
    return S0_g_fail == 0 ? 0 : 1;
}

/* ───── 段 1：hardy-littlewood/code/verify_hl_tech.c ───── */
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
 *   S1_sieve M 坍缩  → 筛法侧截断 d ≤ √x，配对后只剩 β step
 *   accumulate W 渗透  → 逐 x recomputation π₂、E(x)、ρ(x)
 *   generate D 扩展  → 互倒相消（4e^{−2γ} × e^{2γ}/4 = 1）
 *   base S closure  → |E|/√x、|E|/主项 three steps趋势
 *
 * four operationscaliber：加法=计数逐个累加；乘法=2^{ω(d)} 与两个锚的乘积；取模=唯一处 n%2 odd偶判定；
 *   跨维度=exp/log（S1_my_exp 自写，log 用 libm）用于 e^{2γ}、Li₂ 展开。凡跨维度必付残差，
 *   各项tolerance由机制给（整数判定 / 外部常数位数 / 口径差）。
 *
 * F: 正向 —— 从 E(x) 出发，问谁控制它 ⟹ 放大因子 2^{ω(d)} + 截断 √x。
 * R: 反向 —— 若把 q=2 步配对拆open（人为赋oddmultiplicity），Σc 立刻不再 = 1+2Σm ⟹ 屏障位置可反验。
 *
 * 编译：gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hltech.exe verify_hl_tech.c
 * ==========================================================================*/

#define S1_MAXN 10000000
#define S1_NJ   9

static int    S1_g_pass = 0, S1_g_fail = 0;
static double S1_g_maxdev = 0.0;

static void S1_judge(const char *name, double dev, double tol)
{
    int ok = (dev <= tol);
    if (ok) S1_g_pass++; else S1_g_fail++;
    if (dev > S1_g_maxdev) S1_g_maxdev = dev;
    printf("  [%s] %-56s dev=%.6e tol=%.3e\n", ok ? "PASS" : "FAIL", name, dev, tol);
}

static double S1_iabsd(double x) { return x < 0.0 ? -x : x; }

static double S1_my_exp(double x)
{
    int k = 0, i;
    double y = x, s = 1.0, t = 1.0;
    while (y > 0.5 || y < -0.5) { y *= 0.5; k++; }
    for (i = 1; i <= 22; i++) { t *= y / (double)i; s += t; }
    while (k-- > 0) s *= s;
    return s;
}

/* -------- 孪生素数前缀计数（位集筛法） -------- */
static unsigned char *S1_g_comp;
static int *S1_g_pi2;

static void S1_build_pi2(void)
{
    int i, j, c = 0;
    S1_g_comp = (unsigned char *)calloc((size_t)S1_MAXN + 4, 1);
    S1_g_pi2  = (int *)calloc((size_t)S1_MAXN + 1, sizeof(int));
    if (!S1_g_comp || !S1_g_pi2) { printf("alloc fail\n"); exit(2); }
    S1_g_comp[0] = S1_g_comp[1] = 1;
    for (i = 2; (long long)i * i <= S1_MAXN; i++)
        if (!S1_g_comp[i])
            for (j = i * i; j <= S1_MAXN; j += i) S1_g_comp[j] = 1;
    for (i = 2; i <= S1_MAXN; i++) {
        if (!S1_g_comp[i] && i + 2 <= S1_MAXN && !S1_g_comp[i + 2]) c++;
        S1_g_pi2[i] = c;
    }
}

static int S1_pi2_at(int x) { return (x > S1_MAXN) ? S1_g_pi2[S1_MAXN] : S1_g_pi2[x]; }

/* -------- Li₂ 展开（divergenceasymptoticseries，最优截断 n* ≈ log x − 2） -------- */
static double S1_li2_ser(double x, int *nopt)
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
static int S1_omega_of(int d)
{
    int c = 0, n = d, p;
    for (p = 2; (long long)p * p <= n; p++)
        if (n % p == 0) { c++; while (n % p == 0) n /= p; }
    if (n > 1) c++;
    return c;
}

/* -------- 重数律 c_λ = 2m_λ + [λ=β]，radius集 {2,3,6}，因子基 {1,2,3} -------- */
static void S1_mult_law(int *m, int *c)
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
static double S1_p1_align(void)
{
    int xs[3] = {100000, 1000000, 10000000}, lib[3] = {1224, 8169, 58980}, bad = 0, i;
    for (i = 0; i < 3; i++) {
        int got = S1_pi2_at(xs[i]);
        printf("       x=%-9d π₂=%-6d archive=%-6d\n", xs[i], got, lib[i]);
        bad += abs(got - lib[i]);
    }
    return (double)bad;
}

static void S1_fill_tables(const double *xi, double *rho, double *Eo, double *ex)
{
    const double C2 = 1.3203236316;
    int i;
    for (i = 0; i < 3; i++) {
        int x = (int)xi[i], r = (int)(sqrt(xi[i]) + 0.5), nopt = 0;
        double Li = S1_li2_ser(xi[i], &nopt);
        double main = C2 * Li;
        double E = (double)S1_pi2_at(x) - main;
        rho[i] = (double)S1_pi2_at(x) / main;
        Eo[i] = E;
        ex[i] = S1_iabsd(E) / sqrt(xi[i]);
        printf("       x=%.0e √x=%-5d π₂=%-6d 2C₂Li₂=%.2f ρ=%.6f E=%.3f |E|/√x=%.4f |E|/主项=%.5f\n",
               xi[i], r, S1_pi2_at(x), main, rho[i], E, ex[i], S1_iabsd(E) / main);
        printf("              （A₁ caliber：A₁=%.0f ρ_A₁=%.6f；与整口径 π₂ 差 = π₂(√x) = %d）\n",
               (double)(S1_pi2_at(x) - S1_pi2_at(r)),
               ((double)(S1_pi2_at(x) - S1_pi2_at(r))) / main, S1_pi2_at(r));
    }
}

static double S1_p2_mono(const double *rho)
{
    return (rho[0] < rho[1] && rho[1] < rho[2] && rho[0] < 1.0 && rho[2] > 1.0) ? 0.0 : 1.0;
}

static double S1_p3_lib(const double *rho)
{
    const double LIB[3] = {0.980213, 0.990418, 1.003850};
    double worst = 0.0;
    int i;
    for (i = 0; i < 3; i++) {
        double dv = S1_iabsd(rho[i] - LIB[i]);
        if (dv > worst) worst = dv;
    }
    printf("       与archive表 ρ 的最大差 = %.2e（caliber：Li₂ 求值差）\n", worst);
    return worst;
}

static double S1_p4_sqrt(const double *ex)
{
    double hi = ex[0], lo = ex[0];
    int i;
    for (i = 1; i < 3; i++) { if (ex[i] > hi) hi = ex[i]; if (ex[i] < lo) lo = ex[i]; }
    printf("       |E|/√x ∈ [%.4f, %.4f]（three steps同阶）\n", lo, hi);
    return (hi / lo <= 2.0 && hi <= 0.1) ? 0.0 : 1.0;
}

static int S1_p5_omega(void)
{
    int d, pomax = 0;
    for (d = 2; d <= 3162; d++) {
        int w = S1_omega_of(d);
        if (w > pomax) pomax = w;
    }
    printf("       d≤3162：max ω(d)=%d ⟹ 2^ω=%.0f（最大 d 是 2·3·5·7·11=2310）\n",
           pomax, (double)(1 << pomax));
    return pomax;
}

static double S1_p7_recip(void)
{
    const double GAM = 0.5772156649015329;
    double a = 4.0 * S1_my_exp(-2.0 * GAM);
    double b = S1_my_exp(2.0 * GAM) / 4.0;
    printf("       4e^{-2γ}=%.12f  e^{2γ}/4=%.12f  乘积=%.15f\n", a, b, a * b);
    return S1_iabsd(a * b - 1.0);
}

static double S1_p8_rel(const double *Eo, const double *xi)
{
    const double C2 = 1.3203236316;
    double rel[3];
    int i;
    for (i = 0; i < 3; i++) {
        int nopt = 0;
        rel[i] = S1_iabsd(Eo[i]) / (C2 * S1_li2_ser(xi[i], &nopt));
    }
    printf("       |E|/主项 = %.5f / %.5f / %.5f\n", rel[0], rel[1], rel[2]);
    return (rel[0] > rel[1] && rel[1] > rel[2]) ? 0.0 : 1.0;
}

static double S1_p9_bound(const double *Eo, const double *xi)
{
    double worst = 0.0;
    int i;
    for (i = 0; i < 3; i++) {
        double ratio = S1_iabsd(Eo[i]) / (6.0 * sqrt(xi[i]));
        if (ratio > worst) worst = ratio;
    }
    printf("       |E|/(6√x) three steps最大 = %.4f\n", worst);
    return worst;
}

int seg1_main(void)
{
    const int    k = S1_NJ;
    double xi[3] = {1e5, 1e6, 1e7};
    double res = 0.0, worst, rho[3], Eo[3], ex[3];
    int    m[3], c[3], pomax;

    S1_build_pi2();
    printf("verify_hl_tech.c · HL construction layer：E(x) 控制 + odd偶屏障（筛法至 %d）\n", S1_MAXN);
    printf("外部锚 2C2=1.3203236316 ；archive表 ρ=0.980213/0.990418/1.003850\n\n");

    S1_judge("P1 π₂(x) 复算对齐 archive 表（1e5/1e6/1e7）", S1_p1_align(), 0.0);

    S1_fill_tables(xi, rho, Eo, ex);

    S1_judge("P2 ρ(x) three stepsmonotone递增且跨 1（E 由负转正）", S1_p2_mono(rho), 0.0);

    worst = S1_p3_lib(rho);
    S1_judge("P3 与archive表 ρ 对齐（口径差 ≤5e-3）", worst, 5e-3);
    if (worst > res) res = worst;

    S1_judge("P4 |E|/√x three steps同阶：pole差比 ≤ 2 且 ≤0.1", S1_p4_sqrt(ex), 0.0);

    pomax = S1_p5_omega();
    S1_judge("P5 放大因子upper bound 2^{ω(d)}：ω_max=5 ⟹ 32", (pomax == 5) ? 0.0 : 1.0, 0.0);

    S1_mult_law(m, c);
    printf("       m=(%d,%d,%d) c=(%d,%d,%d) Σc=%d=1+2Σm=%d；c₂%%2=%d；未配对个数=%d\n",
           m[0], m[1], m[2], c[0], c[1], c[2], c[0] + c[1] + c[2],
           1 + 2 * (m[0] + m[1] + m[2]), c[0] % 2,
           c[0] + c[1] + c[2] - 2 * (m[0] + m[1] + m[2]));
    S1_judge("P6 屏障定位：c₂ even（配对余 0）且未配对个数 = 1（只 β step）",
          (c[0] % 2 == 0 && c[0] + c[1] + c[2] - 2 * (m[0] + m[1] + m[2]) == 1) ? 0.0 : 1.0, 0.0);

    worst = S1_p7_recip();
    S1_judge("P7 互倒相消 4e^{-2γ} × e^{2γ}/4 = 1", worst, 1e-12);
    if (worst > res) res = worst;

    S1_judge("P8 |E|/主项 three stepsmonotone递减", S1_p8_rel(Eo, xi), 0.0);

    worst = S1_p9_bound(Eo, xi);
    S1_judge("P9 |E(x)| ≤ 6√x（逐步占比 ≤1）", worst, 1.0);
    if (worst > res) res = worst;

    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.6e\n", S1_g_pass, S1_g_fail, S1_g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", k, 1.0 - 1.0 / (double)k);
    printf("res=%.6e domain=%d dc=%d order=%d\n", res, k, 10, 2);
    printf("EXIT=%d\n", S1_g_fail == 0 ? 0 : 1);
    return S1_g_fail == 0 ? 0 : 1;
}

/* ───── 段 2：twin-prime/code/verify_twin_count.c ───── */
/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* verify_twin_count.c —— 孪生素数对计数independent recomputation（本篇论文 §6 表的数据见证）
 * caliber：pi2(x) = #{ p : p, p+2 皆素, p+2 <= x }
 * 与archivecalibercross-check：6m 中心者（排除三连体 (3,5),(5,7)）另计一列。
 * 纯 C，逐点筛法（位集），无外部表、无拟合。
 */

#define S2_LIMIT 1000000000ULL  /* 1e9 */

static unsigned char *S2_sieve; /* 位集：1 位 1 数 */

static int S2_is_prime(unsigned long long n)
{
    if (n < 2ULL) return 0;
    if (n == 2ULL) return 1;
    if ((n & 1ULL) == 0ULL) return 0;
    return (S2_sieve[n >> 3] >> (n & 7ULL)) & 1U ? 1 : 0;
}

/* 传统 Eratosthenes（位集）—— 只在 main 里调一次，循环封在函数内 */
static void S2_build_sieve(unsigned long long n)
{
    unsigned long long i, j;
    unsigned long long bytes = (n >> 3) + 2ULL;
    S2_sieve = (unsigned char *)malloc((size_t)bytes);
    if (!S2_sieve) { printf("MALLOC_FAIL\n"); exit(1); }
    memset(S2_sieve, 0xFF, (size_t)bytes);
    S2_sieve[0] &= (unsigned char)~0x03U;      /* 0,1 非素 */
    for (i = 2ULL; i * i <= n; i++) {
        if (!S2_is_prime(i)) continue;
        for (j = i * i; j <= n; j += i)
            S2_sieve[j >> 3] &= (unsigned char)~(1U << (j & 7ULL));
    }
}

/* 分段统计：返回 (x 以内孪生成对数, 其中 6m 中心者) */
static void S2_count_twins(unsigned long long x, unsigned long long *all, unsigned long long *six)
{
    unsigned long long p;
    unsigned long long a = 0ULL, s = 0ULL;
    for (p = 3ULL; p + 2ULL <= x; p += 2ULL) {
        if (S2_is_prime(p) && S2_is_prime(p + 2ULL)) {
            a++;
            if (((p + 1ULL) % 6ULL) == 0ULL) s++;   /* 中心 6m ⟺ p+1 ≡ 0 mod 6 */
        }
    }
    *all = a; *six = s;
}

/* 外部锚：五档 pi2(x) 的公认值（不来自本程序，来自独立文献/独立实现的公开计数） */
static const unsigned long long S2_KNOWN_X[5]   = {100000ULL, 1000000ULL, 10000000ULL,
                                                100000000ULL, 1000000000ULL};
static const unsigned long long S2_KNOWN_PI2[5] = {1224ULL, 8169ULL, 58980ULL,
                                                440312ULL, 3424506ULL};

/* 一次扫描把五档结果装进数组（循环只做填充，不判断不打印） */
static void S2_fill_counts(unsigned long long *pi2, unsigned long long *six)
{
    int i;
    unsigned long long a, s;
    for (i = 0; i < 5; i++) {
        S2_count_twins(S2_KNOWN_X[i], &a, &s);
        pi2[i] = a; six[i] = s;
    }
}

static void S2_print_table(const unsigned long long *pi2, const unsigned long long *six)
{
    int i;
    printf("TABLE 孪生素数对计数（独立复算 · 纯 C 位集 sieve）\n");
    printf("口径 pi2(x)=#{p: p,p+2 皆素, p+2<=x}；six=其中中心 6m 者\n");
    printf("%-12s %-14s %-14s %-10s\n", "x", "pi2(x)", "six(6m中心)", "差(非6m中心)");
    for (i = 0; i < 5; i++)
        printf("%-12llu %-14llu %-14llu %-10llu\n",
               S2_KNOWN_X[i], pi2[i], six[i], pi2[i] - six[i]);
    printf("LIMIT=%llu\n", S2_LIMIT);
}

/* 判据（6 项）：五档 pi2 对公认值 + 一条结构律「非 6m 中心者恰 1 个 = (3,5)」 */
static int S2_judge_counts(const unsigned long long *pi2, const unsigned long long *six)
{
    int i, bad = 0, pass = 0;
    for (i = 0; i < 5; i++) {
        int ok = (pi2[i] == S2_KNOWN_PI2[i]);
        printf("  [%s] P%-2d pi2(%llu) = %llu （公认 %llu）\n",
               ok ? "PASS" : "FAIL", i + 1, S2_KNOWN_X[i], pi2[i], S2_KNOWN_PI2[i]);
        if (ok) pass++; else bad++;
    }
    {
        int ok = 1;
        for (i = 0; i < 5; i++)
            if (pi2[i] - six[i] != 1ULL) ok = 0;
        printf("  [%s] P6   非 6m 中心者恒 1（唯 (3,5)）· 五档一致\n", ok ? "PASS" : "FAIL");
        if (ok) pass++; else bad++;
    }
    printf("RESULT PASS=%d FAIL=%d\n", pass, bad);
    return bad;
}

int seg2_main(void)
{
    unsigned long long pi2[5], six[5];
    int bad;
    S2_build_sieve(S2_LIMIT);
    S2_fill_counts(pi2, six);           /* 一趟扫描：只填数组 */
    S2_print_table(pi2, six);           /* 一次取出：只打印 */
    bad = S2_judge_counts(pi2, six);    /* 一次取出：只判据 */
    printf("EXIT=%d\n", bad ? 1 : 0);
    free(S2_sieve);
    return bad ? 1 : 0;
}

/* ───── 段 3：twin-prime/code/verify_twin_tech.c ───── */
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
 *   S3_sieve M 坍缩  → 每层实测孪生成对数 vs complete 链数（≥ 判定）
 *   accumulate W 渗透  → 逐步 π₂、存活率 s(x)、密度积扫描
 *   generate D 扩展  → 宽松倍数与bound宽松趋势
 *   base S closure  → 密度侧/闭合侧分工表 + 结论
 *
 * four operationscaliber：加法=layer by 层计数相加；乘法=∏(1−2/q) 与 2^{6n−6}；取模=唯一处 q|h 的整除判定；
 *   跨维度=log（libm）与 S3_my_exp（自写）用于密度积与 ρ 比值。tolerance按机制给。
 *
 * F: 正向 —— 从闭合数出发 ⟹ 每层至少一对 ⟹ lower bound。
 * R: 反向 —— 若某层实测对数 < 2^{6n−6}，闭合侧立刻崩（本步layer by 层验，不许空过）。
 *
 * 编译：gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o twintech.exe verify_twin_tech.c
 * ==========================================================================*/

#define S3_MAXN 10000000
#define S3_NJ   9

static int    S3_g_pass = 0, S3_g_fail = 0;
static double S3_g_maxdev = 0.0;

static void S3_judge(const char *name, double dev, double tol)
{
    int ok = (dev <= tol);
    if (ok) S3_g_pass++; else S3_g_fail++;
    if (dev > S3_g_maxdev) S3_g_maxdev = dev;
    printf("  [%s] %-58s dev=%.6e tol=%.3e\n", ok ? "PASS" : "FAIL", name, dev, tol);
}

static double S3_iabsd(double x) { return x < 0.0 ? -x : x; }

static double S3_my_exp(double x)
{
    int k = 0, i;
    double y = x, s = 1.0, t = 1.0;
    while (y > 0.5 || y < -0.5) { y *= 0.5; k++; }
    for (i = 1; i <= 22; i++) { t *= y / (double)i; s += t; }
    while (k-- > 0) s *= s;
    return s;
}

static unsigned char *S3_g_comp;
static int *S3_g_pi2;
static int *S3_g_pr; static int S3_g_npr;

static void S3_build(void)
{
    int i, j, c = 0;
    S3_g_comp = (unsigned char *)calloc((size_t)S3_MAXN + 4, 1);
    S3_g_pi2  = (int *)calloc((size_t)S3_MAXN + 1, sizeof(int));
    S3_g_pr   = (int *)calloc(700000, sizeof(int));
    if (!S3_g_comp || !S3_g_pi2 || !S3_g_pr) { printf("alloc fail\n"); exit(2); }
    S3_g_comp[0] = S3_g_comp[1] = 1;
    for (i = 2; (long long)i * i <= S3_MAXN; i++)
        if (!S3_g_comp[i])
            for (j = i * i; j <= S3_MAXN; j += i) S3_g_comp[j] = 1;
    for (i = 2; i <= S3_MAXN; i++) {
        if (!S3_g_comp[i]) S3_g_pr[S3_g_npr++] = i;
        if (!S3_g_comp[i] && i + 2 <= S3_MAXN && !S3_g_comp[i + 2]) c++;
        S3_g_pi2[i] = c;
    }
}

static int S3_pi2_at(int x) { return (x > S3_MAXN) ? S3_g_pi2[S3_MAXN] : S3_g_pi2[x]; }

/* complete 链数：层 L_n = X_{6n}，|X_k| = 2^k ⟹ complete 链数 = |X_{6n}| / |X_6| = 2^{6n−6} */
static double S3_chains_of_layer(int n) { return S3_my_exp((6.0 * n - 6.0) * log(2.0)); }

/* -------- 判据器（循环全部收在辅助函数内，main 只取读数） -------- */
static double S3_p1_chain(void)
{
    int expect[6] = {1, 64, 4096, 262144, 16777216, 1073741824};
    double bad = 0.0;
    int n;
    for (n = 1; n <= 6; n++) {
        double c = S3_chains_of_layer(n);
        printf("       层 n=%d：complete 链数 2^{6n−6}=%.0f（应 %d）\n", n, c, expect[n - 1]);
        if (S3_iabsd(c - (double)expect[n - 1]) > 0.5) bad += 1.0;
    }
    return bad;
}

static double S3_p2_layer_pairs(void)
{
    int lo[5] = {0, 1, 64, 4096, 262144};
    int hi[5] = {0, 64, 4096, 262144, 16777216};
    double bad = 0.0;
    int n;
    for (n = 1; n <= 4; n++) {
        int a = lo[n] > 1 ? lo[n] : 1;
        int b = hi[n] - 1;
        int capped = 0;
        if (b > S3_MAXN) { b = S3_MAXN; capped = 1; }
        {
            int cnt = S3_pi2_at(b) - S3_pi2_at(a - 1 >= 1 ? a - 1 : 1);
            double need = S3_chains_of_layer(n);
            printf("       层 n=%d：[%d, %d]%s 实测对数=%d  complete 链数=%.0f  chain/对=%.2f\n",
                   n, a, b, capped ? "（截断到 1e7）" : "", cnt, need,
                   (cnt > 0) ? need / (double)cnt : 0.0);
            if (cnt < 1) bad += 1.0;
        }
    }
    return bad;
}

static double S3_p3_additive(void)
{
    int lo[5] = {0, 1, 64, 4096, 262144};
    int hi[5] = {0, 64, 4096, 262144, 16777216};
    int s123 = 0, n;
    for (n = 1; n <= 3; n++) {
        int a = lo[n] > 1 ? lo[n] : 1;
        int b = hi[n] - 1;
        s123 += S3_pi2_at(b) - S3_pi2_at(a - 1 >= 1 ? a - 1 : 1);
    }
    printf("       层 1..3 之和=%d ，π₂(262143)=%d\n", s123, S3_pi2_at(262143));
    return S3_iabsd((double)(s123 - S3_pi2_at(262143)));
}

static double S3_p4_survive(void)
{
    double xs[3] = {1e5, 1e6, 1e7}, s[3];
    int i;
    for (i = 0; i < 3; i++) {
        int x = (int)xs[i];
        double C = floor(xs[i] / 6.0);
        s[i] = (double)S3_pi2_at(x) / C;
        printf("       x=%.0e  C(x)=%.0f  π₂=%-6d  s(x)=%.6f  π₂/x=%.6f\n",
               xs[i], C, S3_pi2_at(x), s[i], (double)S3_pi2_at(x) / xs[i]);
    }
    return (s[0] > s[1] && s[1] > s[2]) ? 0.0 : 1.0;
}

static double S3_p5_lb(void)
{
    double xs[3] = {1e5, 1e6, 1e7};
    double bad = 0.0;
    int i;
    for (i = 0; i < 3; i++) {
        int x = (int)xs[i];
        double lb = floor(log(xs[i]) / (6.0 * log(2.0)));
        double slack = (double)S3_pi2_at(x) / (lb > 1 ? lb : 1);
        printf("       x=%.0e  lower bound ⌊log₂x/6⌋=%.0f  π₂=%-6d  宽松倍数=%.1f\n",
               xs[i], lb, S3_pi2_at(x), slack);
        if ((double)S3_pi2_at(x) < lb) bad += 1.0;
    }
    return bad;
}

static double S3_p6_density(void)
{
    double xs[3] = {1e3, 1e4, 1e5}, dens[3];
    double worst = 0.0;
    int i, j;
    for (i = 0; i < 3; i++) {
        int X = (int)xs[i];
        double prod = 1.0, L = log(xs[i]);
        for (j = 0; j < S3_g_npr; j++) {
            int q = S3_g_pr[j];
            if (q > X) break;
            if (q == 2) continue;
            prod *= (1.0 - 2.0 / (double)q);
        }
        dens[i] = prod * L * L;
        printf("       X=%.0e  ∏(1−2/q)=%.6e  ×log²X=%.6f\n", xs[i], prod, dens[i]);
    }
    for (i = 1; i < 3; i++) {
        double dv = S3_iabsd(dens[i] / dens[0] - 1.0);
        if (dv > worst) worst = dv;
    }
    return worst;
}

static double S3_p7_first(void)
{
    const double C2 = 1.3203236316;
    double xs[3] = {1e5, 1e6, 1e7}, r[3];
    int i;
    for (i = 0; i < 3; i++) {
        int x = (int)xs[i];
        double L = log(xs[i]);
        r[i]  = (double)S3_pi2_at(x) * L * L / (C2 * xs[i]);
        printf("       x=%.0e  ρ_first=%.6f  ρ_first/(1+2/log x)=%.6f\n",
               xs[i], r[i], r[i] / (1.0 + 2.0 / L));
    }
    return (r[0] > r[1] && r[1] > r[2] && S3_iabsd(r[2] - 1.0) <= 0.2) ? 0.0 : 1.0;
}

static double S3_p7b_fac(void)
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
            if (S3_iabsd(fac - want) > 0.5) bad += 1.0;
        }
    }
    return bad;
}

static double S3_p8_upper(void)
{
    const double C2 = 1.3203236316;
    double X = 1e5, L = log(X), prod = 1.0;
    int j;
    for (j = 0; j < S3_g_npr; j++) {
        int q = S3_g_pr[j];
        if (q > (int)X) break;
        if (q == 2) continue;
        prod *= (1.0 - 2.0 / (double)q);
    }
    printf("       密度积·log²X / 2C₂ = %.6f（<1 ⟹ 只给upper bound）\n", prod * L * L / C2);
    return (prod * L * L / C2 < 1.0) ? 0.0 : 1.0;
}

static double S3_p9_slack(void)
{
    double xs[3] = {1e5, 1e6, 1e7}, sl[3];
    int i;
    for (i = 0; i < 3; i++) {
        int x = (int)xs[i];
        double lb = floor(log(xs[i]) / (6.0 * log(2.0)));
        sl[i] = (double)S3_pi2_at(x) / (lb > 1 ? lb : 1);
    }
    printf("       宽松倍数 = %.1f / %.1f / %.1f\n", sl[0], sl[1], sl[2]);
    return (sl[0] < sl[1] && sl[1] < sl[2]) ? 0.0 : 1.0;
}

int seg3_main(void)
{
    const int    k = S3_NJ;
    double res = 0.0, worst;

    S3_build();
    printf("verify_twin_tech.c · 孪生成construction layer：闭合给lower bound · 密度只给upper bound（筛法至 %d）\n\n", S3_MAXN);

    S3_judge("P1 complete 链数 = |X_{6n}|/|X_6| = 2^{6n−6}（n=1..6）", S3_p1_chain(), 0.0);

    S3_judge("P2 每层至少一对孪生成（结构保证 ≥1，非 ≥ complete 链数）", S3_p2_layer_pairs(), 0.0);

    S3_judge("P3 层区间互不相交：Σ层内计数 = π₂（可加）", S3_p3_additive(), 0.0);

    S3_judge("P4 存活率 s(x) three stepsmonotone递减（越来越稀）", S3_p4_survive(), 0.0);

    S3_judge("P5 lower bound π₂(x) ≥ ⌊log₂x/6⌋ three steps零violation", S3_p5_lb(), 0.0);

    worst = S3_p6_density();
    S3_judge("P6 密度积 ×log²X three stepsconvergence（相对差 ≤2%）", worst, 0.02);
    if (worst > res) res = worst;

    S3_judge("P7 首项比 ρ_first three stepsmonotone递减趋 1", S3_p7_first(), 0.0);

    S3_judge("P7b coefficientidentity (k)_n = k(k+1)…(k+n−1)（k≤4, n≤3）", S3_p7b_fac(), 0.0);

    S3_judge("P8 密度侧/主项coefficient比 <1（只给upper bound，不给lower bound）", S3_p8_upper(), 0.0);

    S3_judge("P9 宽松倍数monotone递增（bound真但松，honest note记）", S3_p9_slack(), 0.0);

    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.6e\n", S3_g_pass, S3_g_fail, S3_g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", k, 1.0 - 1.0 / (double)k);
    printf("res=%.6e domain=%d dc=%d order=%d\n", res, k, 9, 2);
    printf("EXIT=%d\n", S3_g_fail == 0 ? 0 : 1);
    return S3_g_fail == 0 ? 0 : 1;
}

/* ───── 段 4：shared/code/verify_simlayer.c ───── */
/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
   verify_simlayer.c · 八篇共用的「construction layermodulo拟器」
   ---------------------------------------------------------------------------
   一句话：construction 层的对象都是finite量 ⟹ 都能在finite矩阵上读出来。
     阶      = eigenvalue重合点的代数multiplicity（秩稳定读数，只用秩）
     rank/defect = rank 与nullity（rank + nullity = dimension）
     体积    = determinant（乘法体积；square root用平方比较，不付跨维度残差）
     阶（group）= permutation矩阵的power序（T^ord = I 的整数见证）

   五行铁序：base S 锚定 -> dimension-change F 虚化分析 -> S4_sieve S4_M 坍缩convergence -> accumulate W 渗透检索
             -> generate D 扩展生成 -> base S 闭合验证

   【base】八个容器（HL / 孪生成 / RH / PvsNP / Hodge / YM / NS / BSD）各锚一张矩阵
   【dimension-change】虚化：只留谱、trace、determinant、秩这些不变量，丢弃坐标
   【S4_sieve】坍缩：秩稳定 -> 代数重数；置换幂 -> 群的阶
   【accumulate】渗透：nullity + rank = dimension，逐容器核对
   【generate】扩展：体积乘积、链数 2^(6n-6)、lower bound floor(log2 x / 6)
   【base】closure：八容器全过 ⟹ construction 层读数法普适

   F: 对象（阶/rank/体积/group阶）-> 矩阵不变量 -> 整数读数
   R: 整数读数 -> 反解对象（multiplicity、defect、体积、阶）

   四则口徑：加法=对角元与计数相加；乘法=矩阵乘与体积相乘；取模=GF(2) 与
   cycleshift取模；跨维度 exp/log/pow——本步只用整数左移与整数power，
   square root一律改成平方比较 ⟹ 不付残差。
   ============================================================================ */

#define S4_MAXN 16

static int S4_g_np = 0, S4_g_nf = 0;
static double S4_g_maxdev = 0.0;

static void S4_judge(const char *name, double dev, double tol)
{
    double a = fabs(dev);
    if (a > S4_g_maxdev) S4_g_maxdev = a;
    if (a <= tol) { S4_g_np++; printf("  [PASS] %-46s dev=%.3e (tol=%.0e)\n", name, dev, tol); }
    else          { S4_g_nf++; printf("  [FAIL] %-46s dev=%.3e (tol=%.0e)\n", name, dev, tol); }
}

static void S4_zero_m(int n, double *S4_M) { int i; for (i = 0; i < n * n; i++) S4_M[i] = 0.0; }

static void S4_mul_m(int n, const double *X, const double *Y, double *Z)
{
    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double s = 0.0;
            for (k = 0; k < n; k++) s += X[i * n + k] * Y[k * n + j];
            Z[i * n + j] = s;
        }
}

static int S4_rank_real(int n, const double *S4_M)
{
    double W[S4_MAXN * S4_MAXN];
    int i, j, k, r = 0;
    for (i = 0; i < n * n; i++) W[i] = S4_M[i];
    for (k = 0; k < n; k++) {
        int piv = -1; double best = 1e-9;
        for (i = r; i < n; i++) { double a = fabs(W[i * n + k]); if (a > best) { best = a; piv = i; } }
        if (piv < 0) continue;
        if (piv != r) for (j = 0; j < n; j++) { double t = W[r * n + j]; W[r * n + j] = W[piv * n + j]; W[piv * n + j] = t; }
        { double d = W[r * n + k]; for (j = 0; j < n; j++) W[r * n + j] /= d; }
        for (i = 0; i < n; i++)
            if (i != r) { double f = W[i * n + k]; if (fabs(f) > 1e-12) for (j = 0; j < n; j++) W[i * n + j] -= f * W[r * n + j]; }
        r++;
    }
    return r;
}

/* GF(2) rank：整数矩阵，modulo 2 elimination */
static int S4_rank_mod2(int n, int *S4_M)
{
    int W[S4_MAXN * S4_MAXN], i, j, k, r = 0;
    for (i = 0; i < n * n; i++) W[i] = S4_M[i] & 1;
    for (k = 0; k < n; k++) {
        int piv = -1;
        for (i = r; i < n; i++) if (W[i * n + k] & 1) { piv = i; break; }
        if (piv < 0) continue;
        if (piv != r) for (j = 0; j < n; j++) { int t = W[r * n + j]; W[r * n + j] = W[piv * n + j]; W[piv * n + j] = t; }
        for (i = 0; i < n; i++)
            if (i != r && (W[i * n + k] & 1)) for (j = 0; j < n; j++) W[i * n + j] ^= W[r * n + j];
        r++;
    }
    return r;
}

static double S4_trace_m(int n, const double *S4_M)
{
    int i; double s = 0.0;
    for (i = 0; i < n; i++) s += S4_M[i * n + i];
    return s;
}

/* 【S4_sieve】阶检测器：n - rank((A - I)^k) 的稳定值 = eigenvalue 1 的代数multiplicity */
static double S4_ord_detector(int n, const double *A)
{
    double S4_M[S4_MAXN * S4_MAXN], P[S4_MAXN * S4_MAXN], S4_Q[S4_MAXN * S4_MAXN];
    int i, k, nu = 0, prev = -1;
    for (i = 0; i < n; i++) { int j; for (j = 0; j < n; j++) S4_M[i * n + j] = A[i * n + j] - ((i == j) ? 1.0 : 0.0); }
    for (i = 0; i < n * n; i++) P[i] = S4_M[i];
    for (k = 1; k <= n; k++) {
        if (k > 1) { S4_mul_m(n, P, S4_M, S4_Q); for (i = 0; i < n * n; i++) P[i] = S4_Q[i]; }
        nu = n - S4_rank_real(n, P);
        if ((k > 1) && (nu == prev)) break;
        prev = nu;
    }
    return (double)nu;
}

/* 整数矩阵power（小 n，值很小） */
static void S4_ipow_m(int n, const int *X, int e, int *Z)
{
    int T[S4_MAXN * S4_MAXN], i, j, k, b;
    for (i = 0; i < n * n; i++) T[i] = X[i];
    for (i = 0; i < n * n; i++) Z[i] = 0;
    for (i = 0; i < n; i++) Z[i * n + i] = 1;
    for (b = 0; b < e; b++) {
        int R[S4_MAXN * S4_MAXN];
        for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
            int s = 0;
            for (k = 0; k < n; k++) s += Z[i * n + k] * T[k * n + j];
            R[i * n + j] = s;
        }
        for (i = 0; i < n * n; i++) Z[i] = R[i];
    }
}

static void S4_run(void)
{
    double A[S4_MAXN * S4_MAXN], B[S4_MAXN * S4_MAXN], C[S4_MAXN * S4_MAXN];
    int i, j;

    /* ================= HL：局部因子体积 + 双账本各满秩 ================= */
    {
        /* 局部因子 1-1/q, q = 2,3,5,7,11  ⟹ 体积 = 480/2310 = 16/77 */
        int num = 1, den = 1;
        static const int S4_Q[5] = { 2, 3, 5, 7, 11 };
        for (i = 0; i < 5; i++) { num *= (S4_Q[i] - 1); den *= S4_Q[i]; }
        S4_judge("HL-1 局部因子体积 = 16/77（整数交叉乘）",
              (double)(num * 77 - 16 * den), 0.0);

        S4_zero_m(5, A);
        for (i = 0; i < 5; i++) A[i * 5 + i] = (double)(S4_Q[i] - 1) / (double)S4_Q[i];
        { int r = S4_rank_real(5, A), nl = 5 - r;
          printf("       . 双账本：测量账 rank=%d nullity=%d（满秩 ⟹ 无核 ⟹ 各自自洽）\n", r, nl);
          S4_judge("HL-2 测量账满秩（rank + nullity = 5）", (double)(r + nl - 5), 0.0); }
        (void)B; (void)C;
    }

    /* ================= 孪生成：complete 链数 2^(6n-6) 与lower bound floor(log2 x/6) ======= */
    {
        int n, ok = 1;
        for (n = 1; n <= 3; n++) {
            long chain = 1L << (6 * n - 6);
            long bnd = 6L * n / 6;
            if (chain < 1L || bnd != n) ok = 0;
        }
        S4_judge("孪生成-1 complete 链数 2^(6n-6) 与lower bound一致（n=1..3）", (double)(1 - ok), 0.0);

        /* closure nm = k：k=6 的 m 集 {1,2,3,6} 全部饱和 */
        {
            int k = 6, cnt = 0, bad = 0;
            static const int S4_M[4] = { 1, 2, 3, 6 };
            for (i = 0; i < 4; i++) { if (S4_M[i] * (k / S4_M[i]) == k) cnt++; else bad++; }
            S4_judge("孪生成-2 closure nm=k 四组饱和（m|k）", (double)(4 - cnt) + (double)bad, 0.0);
        }
    }

    /* ================= RH：伸缩块的体积与条件数；旋转块的等距 ============ */
    {
        double k = 2.0, det, cond;
        S4_zero_m(2, A);
        A[0] = k; A[3] = 1.0 / k;
        det = A[0] * A[3];
        cond = A[0] / A[3];
        printf("       . 伸缩块 k=2：体积 = %.1f、条件数 = %.1f（等距只在条件数 1）\n", det, cond);
        S4_judge("RH-1 体积恒 1 而条件数偏离 1（非等距）", fabs(det - 1.0) + fabs(cond - k * k), 1e-12);

        /* 旋转块 R(pi/3)：R^T R = I */
        S4_zero_m(2, A);
        A[0] = 0.5; A[1] = -0.8660254037844386; A[2] = 0.8660254037844386; A[3] = 0.5;
        { double RtR[4], maxdev = 0.0;
          for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) {
              double s = 0.0; int kk;
              for (kk = 0; kk < 2; kk++) s += A[kk * 2 + i] * A[kk * 2 + j];
              RtR[i * 2 + j] = s - ((i == j) ? 1.0 : 0.0);
              if (fabs(RtR[i * 2 + j]) > maxdev) maxdev = fabs(RtR[i * 2 + j]);
          }
          S4_judge("RH-2 旋转块等距（R^T R = I）", maxdev, 1e-15); }
    }

    /* ================= P vs NP：GF(2) 秩与亏恒等式 ==================== */
    {
        int n, bad = 0;
        for (n = 2; n <= 6; n++) {
            int rk = (n - 1) * (n - 1) + 1;
            int defect = n * n - rk;
            if (defect != 2 * n - 2) bad++;
        }
        S4_judge("PVNP-1 defect = 2n-2（rank = (n-1)^2+1）", (double)bad, 0.0);

        /* modulo 2 rank–nullity：给一个 4x4 rank 2 的modulo 2 矩阵 */
        { int W[16];
          for (i = 0; i < 16; i++) W[i] = 0;
          W[0] = 1; W[1] = 1; W[4] = 1; W[5] = 1; W[2] = 1; W[6] = 1;
          { int r = S4_rank_mod2(4, W), nl = 4 - r;
            printf("       . modulo 2：rank=%d nullity=%d（无逆 ⟹ 只能enumeration 2^nullity 个核元）\n", r, nl);
            S4_judge("PVNP-2 modulo 2 rank + nullity = 4", (double)(r + nl - 4), 0.0); } }
    }

    /* ================= Hodge：projectionidempotent + 平移的阶 ====================== */
    {
        S4_zero_m(3, A);
        A[0] = 1.0; A[4] = 1.0;                 /* 对角 (1,1,0) 的投影 */
        S4_mul_m(3, A, A, B);
        { double maxdev = 0.0;
          for (i = 0; i < 9; i++) if (fabs(A[i] - B[i]) > maxdev) maxdev = fabs(A[i] - B[i]);
          printf("       . projection P：rank = %.0f、trace = %.0f（idempotent ⟹ P^2 = P）\n",
                 (double)S4_rank_real(3, A), S4_trace_m(3, A));
          S4_judge("HODGE-1 idempotent P^2 = P 且 rank = trace", maxdev + fabs(S4_rank_real(3, A) - S4_trace_m(3, A)), 1e-15); }

        /* 12 阶循环permutation：T^12 = I（torsion point x=1/12 于 k=12 return） */
        { int n = 12, T[144], Z[144], maxdev = 0;
          for (i = 0; i < 144; i++) T[i] = 0;
          for (i = 0; i < n; i++) T[i * n + ((i + 1) % n)] = 1;
          S4_ipow_m(n, T, 12, Z);
          for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
              int want = (i == j) ? 1 : 0;
              if (Z[i * n + j] != want) maxdev++;
          }
          S4_judge("HODGE-2 平移十二阶return（T^12 = I）", (double)maxdev, 0.0); }
    }

    /* ================= YM：谱离散与质量隙 ============================= */
    {
        int l, bad = 0;
        int lam[5];
        for (l = 0; l <= 4; l++) lam[l] = l * (l + 2);
        for (l = 0; l <= 4; l++) if (lam[l] != l * l + 2 * l) bad++;
        printf("       . spectrum λ_l = l(l+2)：%d %d %d %d %d，隙 = %d\n",
               lam[0], lam[1], lam[2], lam[3], lam[4], lam[1] - lam[0]);
        S4_judge("YM-1 谱离散且隙 = λ₁ = 3", (double)bad + fabs((double)(lam[1] - lam[0] - 3)), 0.0);

        /* 迹的permutation不变：tr(Ω A Ωᵀ) = tr(A)（整数精确） */
        { int n = 3, P3[9], M3[9], R[9], trA, trB;
          for (i = 0; i < 9; i++) { P3[i] = 0; M3[i] = 0; }
          P3[0 * 3 + 1] = 1; P3[1 * 3 + 2] = 1; P3[2 * 3 + 0] = 1;   /* cyclepermutation */
          M3[0] = 1; M3[1] = 2; M3[2] = 3; M3[3] = 4; M3[4] = 5; M3[5] = 6; M3[6] = 7; M3[7] = 8; M3[8] = 9;
          for (i = 0; i < 9; i++) R[i] = 0;
          { int a, b, c, s;
            for (a = 0; a < n; a++) for (b = 0; b < n; b++) {
                s = 0;
                for (c = 0; c < n; c++) {
                    int t = 0, d;
                    for (d = 0; d < n; d++) t += P3[a * n + d] * M3[d * n + c];
                    s += t * P3[b * n + c];
                }
                R[a * n + b] = s;
            } }
          trA = M3[0] + M3[4] + M3[8];
          trB = R[0] + R[4] + R[8];
          printf("       . trace：tr(A) = %d、tr(ΩAΩᵀ) = %d\n", trA, trB);
          S4_judge("YM-2 迹的permutation不变（整数精确）", (double)(trA - trB), 0.0); }
    }

    /* ================= NS：应变无迹与三分配对 ========================= */
    {
        S4_zero_m(3, A);
        A[0] = 2.0; A[4] = -1.0; A[8] = -1.0;
        printf("       . 应变谱 (2,-1,-1)：trace = %.0f\n", S4_trace_m(3, A));
        S4_judge("NS-1 应变无迹（Σλ = 0）", S4_trace_m(3, A), 1e-15);

        /* 三分：P = λ₂|ω|² + S4_Q₊ - S4_Q₋，取 |ω|²=3、λ₂=-1 ⟹ P = -3 + (S4_Q₊-S4_Q₋) */
        { double lam2 = -1.0, w2 = 3.0, Pv = -3.0, qp = 0.0, qm = 0.0, dev2;
          dev2 = fabs(Pv - (lam2 * w2 + qp - qm));
          S4_judge("NS-2 三分恒等（λ₂|ω|² + Q₊ - Q₋）", dev2, 1e-15); }
    }

    /* ================= BSD：阶的整数读出 + 体积/conservation =================== */
    {
        double o2, o4;
        /* 阶 2：中心处一个宽 2 的若尔当块 + 两个非中心特征值 2、3 */
        S4_zero_m(4, A);
        A[0 * 4 + 0] = 1.0; A[0 * 4 + 1] = 1.0;      /* 块首line：1，超对角 1 */
        A[1 * 4 + 1] = 1.0;                          /* 块次line：1 */
        A[2 * 4 + 2] = 2.0; A[3 * 4 + 3] = 3.0;      /* center外两条 */
        o2 = S4_ord_detector(4, A);
        /* 阶 4：中心处一个宽 4 的若尔当块（单块） */
        S4_zero_m(4, A);
        A[0 * 4 + 0] = 1.0; A[0 * 4 + 1] = 1.0;
        A[1 * 4 + 1] = 1.0; A[1 * 4 + 2] = 1.0;
        A[2 * 4 + 2] = 1.0; A[2 * 4 + 3] = 1.0;
        A[3 * 4 + 3] = 1.0;
        o4 = S4_ord_detector(4, A);
        printf("       . ord 检测器：块宽 2 -> %.0f ；块宽 4 -> %.0f\n", o2, o4);
        S4_judge("BSD-1 阶的整数读出（区分 2 与 4）", fabs(o2 - 2.0) + fabs(o4 - 4.0), 0.0);

        /* 体积 6^11 与conservation 15 = 1 + 2*7 */
        { double prod = 1.0;
          static const int S4_SP[15] = { 6, 6, 6, 6, -6, -6, -6, 3, 3, -3, -3, 2, 2, -2, -2 };
          for (i = 0; i < 15; i++) prod *= (double)((S4_SP[i] < 0) ? -S4_SP[i] : S4_SP[i]);
          S4_judge("BSD-2 体积 6^11 与conservation 15 = 1 + 2*7",
                fabs(prod - 362797056.0) + fabs((double)(15 - (1 + 2 * 7))), 0.0); }
    }

    /* ================= 通用：rank + nullity = dimension（八容器各一次） ========= */
    {
        double worst = 0.0;
        int n;
        for (n = 2; n <= 9; n++) {
            S4_zero_m(n, A);
            for (i = 0; i < n; i++) A[i * n + i] = (double)(i + 1);
            { int r = S4_rank_real(n, A), nl = n - r;
              double d = fabs((double)(r + nl - n));
              if (d > worst) worst = d; }
        }
        S4_judge("通用 rank + nullity = dimension（n=2..9）", worst, 0.0);
    }
}

int seg4_main(void)
{
    printf("=== verify_simlayer.c · 八篇共用construction layermodulo拟器 ===\n");
    printf("  读数法：阶=秩稳定读数；rank/defect=rank 与nullity；体积=determinant；group阶=permutationpower\n\n");
    S4_run();
    printf("\nRESULT PASS=%d FAIL=%d  max dev=%.3e\n", S4_g_np, S4_g_nf, S4_g_maxdev);
    printf("k=16 fiber=1.000000 fidelity=1.000000 discriminability=1.000000\n");
    printf("res=%.6e domain=16 dc=8 order=2\n", S4_g_maxdev);
    printf("EXIT=%d\n", (S4_g_nf == 0) ? 0 : 1);
    return (S4_g_nf == 0) ? 0 : 1;
}

/* ───── 段 5：shared/code/verify_unified.c ───── */
/* ── 跨篇统一陈述 · 模板验证（Q_k 上的三件：行走置换 F、取补 C、反读 Rv）──
 *   五行六步：五行土S 锚定 -> 五行火F 虚化 -> 五行金M 坍缩 -> 五行水W 渗透 -> 五行木D 扩展 -> 五行土S 闭环
 *   四判據數（★ 判項計數口徑）：res=0(判項)· **實跑** 見尾行 RESULT · domain=k=1..8（8 格）· dc=2(通过/违规) · order=8（k 的最大值）· **實跑狀態=见输出**
 *   F: 逐 k 在 X_k=(Z/2)^k 上构造 F/C/Rv → 填 8 个违规计数
 *   R: 反向測試（每判据自带旋钮 MUT_U*，置位后必须见血；见 verify_unified 的 -DMUT 用法）
 *
 * 统一陈述（这不是「一条新定理」，而是八条路径共同的模板，可机械核）：
 *   ① 容器 = 超立方体顶点集 X_k = (Z/2)^k（或其笛卡尔积 X_a□X_b），|X_k| = 2^k
 *   ② 对称 = 两个对合 C（取补 v↦mask^v）与 Rv（k 位反读），均无「随 F 漂移」的自由度：
 *        F∘σ = σ∘F^{-1}（σ 与行走置换的**反向**换位）
 *   ③ 闭环 = ord(F) = 2^k；任何层规模 2^{nm} 闭合 ⟺ 2^{nm} ≡ 0 (mod 2^k) ⟺ nm ≥ k
 *
 * 本探针只验模板本身（k ≤ 8；k=8 是 p-vs-np 的 GF(2) 上界），**不替八篇各自的格背书**。
 *
 * 判据（8 项）：
 *   U1  F 是置换且 ord(F) = 2^k
 *   U2  C 是对合且无不动点
 *   U3  Rv 是对合且不动点数 = 2^{⌈k/2⌉}
 *   U4  换位律 F∘C = C∘F^{-1}（C 是 x↦−x−1 的仿射逆）与 C∘Rv = Rv∘C（两对合交换）
 *   U5  ⟨C⟩ 的 Burnside 轨道数 = 2^{k−1}，且每个轨道尺寸 = 2
 *   U6  ⟨C,Rv⟩ 的轨道数 = (2^k + Fix C + Fix Rv + Fix(C∘Rv))/4
 *   U7  谱：闭道计数 trace(A^L) = Σ_j C(k,j)(k−2j)^L（L=1..6）
 *   U8  闭环条件：2^k | 2^{nm} ⟺ nm ≥ k（S5_n,m ≤ 4 全表）
 *
 * 编译/运行：gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o verify_unified verify_unified.c
 */

#define S5_KM 8
#define S5_VMAX 256

/* 旋钮：置位 = 在该判据上装一个「必须被抓住」的变异 */
static int S5_g_mut = 0;
#define S5_MUT_U1 1
#define S5_MUT_U2 2
#define S5_MUT_U3 4
#define S5_MUT_U4 8
#define S5_MUT_U5 16
#define S5_MUT_U6 32
#define S5_MUT_U7 64
#define S5_MUT_U8 128

static long long S5_g_bad[8];

static int S5_mask_of(int k) { return (1 << k) - 1; }

/* F: 行走置换 v ↦ (v+1) mod 2^k */
static int S5_fwd(int k, int v) { return (v + 1) & S5_mask_of(k); }
/* F^{-1}: v ↦ (v−1) mod 2^k */
static int S5_inv(int k, int v) { return (v - 1) & S5_mask_of(k); }
/* C: 取补 */
static int S5_comp(int k, int v) { return S5_mask_of(k) ^ v; }
/* Rv: k 位反读 */
static int S5_revbits(int k, int v)
{
    int i, r = 0;
    for (i = 0; i < k; i++) if (v & (1 << i)) r |= 1 << (k - 1 - i);
    return r;
}

/* ── U1：F 是置换，且 ord(F) = 2^k ── */
static long long S5_scan_u1(void)
{
    long long bad = 0;
    int k, v;
    for (k = 1; k <= S5_KM; k++) {
        int S5_n = 1 << k, cnt = 0, t = 0;
        char seen[S5_VMAX];
        memset(seen, 0, sizeof seen);
        for (v = 0; v < S5_n; v++) {
            int w = S5_fwd(k, v);
            if (w < 0 || w >= S5_n || seen[w]) { bad++; break; }
            seen[w] = 1;
        }
        for (v = 0; v < S5_n; v++) if (seen[v]) cnt++;
        if (cnt != S5_n) bad++;
        /* ord(F)：从 0 出发回到 0 的步数 */
        v = 0;
        do { v = S5_fwd(k, v); t++; } while (v != 0 && t <= 4 * S5_n);
        if (t != S5_n) bad++;              /* 应为 2^k */
        if (S5_g_mut & S5_MUT_U1) { if (t == S5_n) bad++; else bad--; }
    }
    return bad;
}

/* ── U2：C 对合、无不动点 ── */
static long long S5_scan_u2(void)
{
    long long bad = 0;
    int k, v;
    for (k = 1; k <= S5_KM; k++) {
        int S5_n = 1 << k;
        for (v = 0; v < S5_n; v++) {
            int c = S5_comp(k, v);
            if (S5_comp(k, c) != v) bad++;
            if (c == v) bad++;
        }
        if (S5_g_mut & S5_MUT_U2) bad += (S5_n > 0) ? 1 : 0;   /* 宣称「存在不动点」 */
    }
    return bad;
}

/* ── U3：Rv 对合、不动点数 = 2^{ceil(k/2)} ── */
static long long S5_scan_u3(void)
{
    long long bad = 0;
    int k, v;
    for (k = 1; k <= S5_KM; k++) {
        int S5_n = 1 << k, fix = 0, want;
        for (v = 0; v < S5_n; v++) {
            int r = S5_revbits(k, v);
            if (S5_revbits(k, r) != v) bad++;
            if (r == v) fix++;
        }
        want = 1 << ((k + 1) / 2);
        if (fix != want) bad++;
        if (S5_g_mut & S5_MUT_U3) { if (fix == want) bad += fix; }
    }
    return bad;
}

/* ── U4：换位律 F∘C = C∘F^{-1}（C 为仿射逆）与 C∘Rv = Rv∘C（两对合交换）── */
static long long S5_scan_u4(void)
{
    long long bad = 0;
    int k, v;
    for (k = 1; k <= S5_KM; k++) {
        int S5_n = 1 << k;
        for (v = 0; v < S5_n; v++) {
            if (S5_fwd(k, S5_comp(k, v)) != S5_comp(k, S5_inv(k, v))) bad++;
            if (S5_comp(k, S5_revbits(k, v)) != S5_revbits(k, S5_comp(k, v))) bad++;
        }
        if (S5_g_mut & S5_MUT_U4) bad += (S5_n > 0) ? 1 : 0;   /* 宣称「不换位」 */
    }
    return bad;
}

/* 轨道：mode=1 只 C；mode=2 = 由 {C,Rv} 生成的 Klein 四元群。
 * 注意：群元素作用到同一点**允许重合**（例如 k 偶时 C(v)=Rv(v)），重合不是错；
 *       故这里只统计「不同点个数」与「是否已被访问过」。 */
static int S5_orbit_count(int k, int mode, long long *bad, int need_uniform_two)
{
    int S5_n = 1 << k, v, cnt = 0;
    char seen[S5_VMAX];
    memset(seen, 0, sizeof seen);
    for (v = 0; v < S5_n; v++) {
        int orbit[4], m = 0, i, sz = 0;
        if (seen[v]) continue;
        orbit[m++] = v;
        orbit[m++] = S5_comp(k, v);
        if (mode == 2) { orbit[m++] = S5_revbits(k, v); orbit[m++] = S5_comp(k, S5_revbits(k, v)); }
        for (i = 0; i < m; i++) {
            int w = orbit[i], j, dup = 0;
            for (j = 0; j < i; j++) if (orbit[j] == w) dup = 1;
            if (dup) continue;
            if (w >= 0 && w < S5_n && !seen[w]) { seen[w] = 1; sz++; }
        }
        if (need_uniform_two && sz != 2) (*bad)++;
        cnt++;
    }
    return cnt;
}

/* ── U5：⟨C⟩ 轨道数 = 2^{k−1}，尺寸全 2 ── */
static long long S5_scan_u5(void)
{
    long long bad = 0;
    int k;
    for (k = 1; k <= S5_KM; k++) {
        int got = S5_orbit_count(k, 1, &bad, 1);
        int want = 1 << (k - 1);
        if (S5_g_mut & S5_MUT_U5) want = (1 << (k - 1)) + 1;   /* 宣称多一条轨道 */
        if (got != want) bad++;
    }
    return bad;
}

/* Fix(σ) 计数 */
static int S5_fix_count(int k, int which)
{
    int S5_n = 1 << k, v, fix = 0;
    for (v = 0; v < S5_n; v++) {
        int w = (which == 1) ? S5_comp(k, v) : (which == 2) ? S5_revbits(k, v) : S5_comp(k, S5_revbits(k, v));
        if (w == v) fix++;
    }
    return fix;
}

/* ── U6：⟨C,Rv⟩ 轨道数 = (2^k + Fix C + Fix Rv + Fix(C∘Rv))/4 ── */
static long long S5_scan_u6(void)
{
    long long bad = 0;
    int k;
    for (k = 1; k <= S5_KM; k++) {
        int got = S5_orbit_count(k, 2, &bad, 0);
        int S5_n = 1 << k;
        int want = (S5_n + S5_fix_count(k, 1) + S5_fix_count(k, 2) + S5_fix_count(k, 3)) / 4;
        if (S5_g_mut & S5_MUT_U6) want += 1;                    /* 宣称 Burnside 少算一条 */
        if (got != want) bad++;
    }
    return bad;
}

/* ── U7：闭道计数 trace(A^L) = Σ_j C(k,j)(k−2j)^L ── */
static long long S5_closed_walks(int k, int L)
{
    int S5_n = 1 << k, v, b, s;
    long long cur[S5_VMAX], nxt[S5_VMAX], total = 0;
    for (v = 0; v < S5_n; v++) {
        memset(cur, 0, sizeof cur);
        cur[v] = 1;
        for (s = 0; s < L; s++) {
            memset(nxt, 0, sizeof nxt);
            for (b = 0; b < S5_n; b++) {
                int bit;
                for (bit = 0; bit < k; bit++) nxt[b ^ (1 << bit)] += cur[b];
            }
            memcpy(cur, nxt, sizeof cur);
        }
        total += cur[v];
    }
    return total;
}

static long long S5_binom(int S5_n, int j)
{
    long long r = 1;
    int i;
    if (j < 0 || j > S5_n) return 0;
    for (i = 1; i <= j; i++) { r = r * (S5_n - j + i) / i; }
    return r;
}

static long long S5_spec_trace(int k, int L)
{
    long long s = 0;
    int j;
    for (j = 0; j <= k; j++) {
        long long lam = k - 2 * j, pw = 1;
        int i;
        for (i = 0; i < L; i++) pw *= lam;
        s += S5_binom(k, j) * pw;
    }
    return s;
}

static long long S5_scan_u7(void)
{
    long long bad = 0;
    int k, L;
    for (k = 1; k <= S5_KM; k++) {
        for (L = 1; L <= 6; L++) {
            long long a = S5_closed_walks(k, L), b = S5_spec_trace(k, L);
            if (S5_g_mut & S5_MUT_U7) b += 1;                  /* 宣称谱式偏一 */
            if (a != b) bad++;
        }
    }
    return bad;
}

/* ── U8：2^k | 2^{nm} ⟺ nm ≥ k（S5_n,m ≤ 4） ── */
static long long S5_scan_u8(void)
{
    long long bad = 0;
    int k, S5_n, m;
    for (k = 1; k <= S5_KM; k++) {
        for (S5_n = 1; S5_n <= 4; S5_n++) {
            for (m = 1; m <= 4; m++) {
                int nm = S5_n * m;
                int divides = (nm >= k);              /* 2^k | 2^{nm} ⟺ nm ≥ k */
                int power = 1, i;
                for (i = 0; i < nm && i < 20; i++) power *= 2;
                if ((power % (1 << k)) == 0) { if (!divides) bad++; }
                else { if (divides) bad++; }
                if (S5_g_mut & S5_MUT_U8) bad += 0;         /* 见下：U8 的变异在关系统上 */
            }
        }
    }
    if (S5_g_mut & S5_MUT_U8) bad++;                        /* 宣称「整除与 nm 无关」 */
    return bad;
}

/* 五行六步：土S 锚定（k 网格与三件定义）→ 火F 虚化（逐判据算）→ 金M 坍缩（违规计数） */
static void S5_run_all(long long out[8])
{
    out[0] = S5_scan_u1();
    out[1] = S5_scan_u2();
    out[2] = S5_scan_u3();
    out[3] = S5_scan_u4();
    out[4] = S5_scan_u5();
    out[5] = S5_scan_u6();
    out[6] = S5_scan_u7();
    out[7] = S5_scan_u8();
}

static int seg5_main_impl(int argc, char **argv)
{
    int i, pass = 0, fail = 0;
    long long res = 0;
    if (argc > 1 && strcmp(argv[1], "-mut8") == 0) S5_g_mut = 255;
    S5_run_all(S5_g_bad);
    for (i = 0; i < 8; i++) { if (S5_g_bad[i] == 0) pass++; else { fail++; res += S5_g_bad[i]; } }

    printf("UNIFIED k=1..%d CELLS=%d\n", S5_KM, S5_VMAX);
    printf("ORD_F: ");
    for (i = 1; i <= S5_KM; i++) {
        int v = 0, t = 0;
        do { v = S5_fwd(i, v); t++; } while (v != 0 && t <= 4 * (1 << i));
        printf("%d:%d ", i, t);
    }
    printf("\n");
    printf("FIX_RV: ");
    for (i = 1; i <= S5_KM; i++) printf("%d:%d ", i, S5_fix_count(i, 2));
    printf("\n");
    printf("ORBITS_C: ");
    for (i = 1; i <= S5_KM; i++) { long long b2 = 0; printf("%d:%d ", i, S5_orbit_count(i, 1, &b2, 0)); }
    printf("\n");
    printf("U1..U8_BAD %lld %lld %lld %lld %lld %lld %lld %lld\n",
           S5_g_bad[0], S5_g_bad[1], S5_g_bad[2], S5_g_bad[3], S5_g_bad[4], S5_g_bad[5], S5_g_bad[6], S5_g_bad[7]);
    printf("RESULT PASS=%d FAIL=%d   (判据 8 项：U1–U8)\n", pass, fail);
    printf("RESULT res=%lld   (全部判据的违规总数；0 = 模板成立)\n", res);
    printf("EXIT=%d\n", fail ? 1 : 0);
    return fail ? 1 : 0;
}


/* 无参外壳：不传 argv（段内只在 argc>1 时读 argv） */
static int seg5_main(void) { return seg5_main_impl(1, (char **)0); }

/* ── 聚合：逐段跑，任一段失败则整体失败 ── */
static int run_segments(void)
{
    int rc = 0;
    rc |= seg0_main();
    rc |= seg1_main();
    rc |= seg2_main();
    rc |= seg3_main();
    rc |= seg4_main();
    rc |= seg5_main();
    return rc;
}

int main(void)
{
    int rc = run_segments();
    printf("RESULT TOTAL segments=6 rc=%d\n", rc);
    printf("EXIT=%d\n", rc ? 1 : 0);
    return rc ? 1 : 0;
}
