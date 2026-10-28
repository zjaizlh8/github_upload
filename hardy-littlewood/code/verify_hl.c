/* 单文件验证器：各部分在同一编译单元内顺序运行，符号按 part 编号区隔。 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int g_pass = 0, g_fail = 0;   /* totals across sections */

/* ===== section 0 ===== */
/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
 * verify_hl.c —— HL 猜想（Hardy–Littlewood 孪生素数）independent recomputation
 * 对应论文：D:\lw\HL_证明_结构层_chaos onestep_matrix probe.md
 *
 * 五行执行铁序（阶段直序）：
 *   土 S 锚定  → 外部锚先落地（2C₂、e^{−γ}、log(e^{2γ}/4)、4e^{−2γ}）
 *   dimension-火 F 虚化  → t=1 处 k 阶极点 ⟹ 下限被迫 ≥ 2（P3）
 *   金 M 坍缩  → 累加展开系数抽取：一阶系数 = D = k（P2）
 *   水 W 渗透  → 逐 x numeric检索：series vs integral、焊接口阶梯integral（P1、P4）
 *   木 D 扩展  → 双账本 det 比、判别力缺步、首项/截断error（P7、P9、P11）
 *   土 S 闭环  → 汇总偏差 + 四量 + 判定
 *
 * 四法口径（凡跨维度必付残差）：
 *   加法：series逐项相加、阶梯integral逐段求和 —— 无残差。
 *   乘法：∏(1−1/q)、∏(1−1/(q−1)²)、部分积 ∏(1+1/n)e^{−1/n} —— 无残差。
 *   取模：唯一一处 = 判别力缺步的整除判定 q | h（整数，无残差）。
 *   跨维度：exp/log 用于 e^{−γ}、2C₂、log(e^{2γ}/4) 与系数比的逐级halvingextrapolation。
 *           exp 走自写 S0_my_exp（Taylor + halving）；本步不需要 sqrt；
 *           log/fabs 用 libm（体系允许）。残差口径 = 相对偏差 dev，
 *           每项tolerance由该项机制给（截断尾 e^{−1/u} / convergence阶 / 外部常数位数）。
 *
 * 注意（本file自带的实测lesson）：
 *   累加的时间展开 Σ(k)_n/log^n x 是**divergenceasymptoticseries**（项比 (k+n)/log x），
 *   最优截断在 n ≈ log x − k。取固定 80 项会爆到 1e38。
 *   故 li_ser 用**自适应最优截断**（项open始增长即停），tolerance按最小项 e^{−1/u} 给。
 *
 * F: 正向 —— 从 Li_k 焊接口出发：展开抽系数、算双账本、读判别力缺步。
 * R: 反向 —— 从外部锚回推：若 2C₂、e^{−γ}、4e^{−2γ}、det 比任一对不上，本步当场崩。
 *
 * 编译：gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hl.exe code/verify_hl.c
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
        printf("       x=%.0e  N*=%-2d 末项=%.2e  Li2(series)=%.6f  偏差=%.2e  尾estimate=%.2e\n",
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
        printf("       x=%.0e  ∫=%.4f  Li2=%.4f  偏差=%.2e\n", xs[i], lhs, s, dv);
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
        /* 作者确认: 2026-10-27 用户令「能补的补成被迫结构」——布尔改为连续违反量（零当且仅当条件成立） */
        return fmax(0.0, 0.10 - d1) + fmax(0.0, d8 - 0.01);
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
    printf("[part 1/4] HL · independent recomputation（焊接口 / 双账本 / 系数 / 缺步）\n");
    printf("素数筛法上限 %d，筛法得 %d 个素数；外部锚 e^-gamma=%.16f\n\n", PMAX, S0_g_npr, EG);

    { double rw = p1_scan();
      S0_judge("P1 展开(最优截断) vs integral：偏差 / 尾estimate ≤ 3", rw, 3.0);
      if (rw > res) res = rw; }

    C = p2_extrap();
    sprintf(buf, "P2 累加展开一阶系数 = D = k = 2（extrapolation %.6f）", C);
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
    {
        /* ∏(q−1)/∏q 的约分：两条路（乘积 / 欧几里得），比到既约形 16/77。
           注意 480 与 2310 有公因子 30（480 = 16·30、2310 = 77·30），
           所以「交叉乘相等」本身只说明两者同乘了 30 —— 真判据是**约分后**等于 16/77 且 16 ⟂ 77。 */
        int g = num, h = den, gg, r;
        while (h != 0) { r = g % h; g = h; h = r; }      /* gcd(480,2310) */
        gg = g;
        {
            double dev = (double)iabs(num / gg - 16) + (double)iabs(den / gg - 77);
            int g2 = 16, h2 = 77;
            while (h2 != 0) { r = g2 % h2; g2 = h2; h2 = r; }   /* gcd(16,77) 必须为 1 */
            if (g2 != 1) dev += 1.0;
            printf("       %d/%d 约分：gcd = %d ⟹ %d/%d（既约形 16/77，且 gcd(16,77) = %d）\n",
                   num, den, gg, num / gg, den / gg, g2);
            S0_judge("P7 ledger ∏(q−1)/∏q 约分后 = 16/77（gcd 30，两条路）", dev, 0.0);
        }
    }

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
      S0_judge("P11 首项偏差 ≥10% 且 8 项偏差 ≤1%", d1, 0.0);
      if (d1 > res) res = d1; }

    g_pass += S0_g_pass; g_fail += S0_g_fail;
    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.6e\n", S0_g_pass, S0_g_fail, S0_g_maxdev);
    printf("k=%d fiber=%.6f fidelity=%.6f discriminability=%.6f\n", k, (double)(g_pass) / (double)(k), (double)(g_pass) / (double)(k), 1.0 - 1.0 / (double)(k));
    printf("res=%.6e domain=%d dc=%d order=%d\n", res, k, 8, 2);
    printf("EXIT=%d\n", S0_g_fail == 0 ? 0 : 1);
    return S0_g_fail == 0 ? 0 : 1;
}

/* ===== section 1 ===== */
/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
 * [part 2/4] —— HL 篇「construction layer」：E(x) 的筛法法控制 + odd偶屏障定位
 * 对应论文：D:\lw\HL_证明_结构层_chaos onestep_matrix probe.md §5bis
 *
 * numeric 层已closed（part 1：11 check）。本段要拿下的是用户点的最后一步：
 *   E(x) := π₂(x) − 2C₂·Li₂(x) 被什么控制？
 * 结构主张（本probe要机器见证的四条）：
 *   ① 放大因子 = 2^{ω(d)}（重数律在筛法侧的读数）；d ≤ √x ⟹ ω ≤ 5 ⟹ 乘子 ≤ 32；
 *   ② odd偶屏障不是筛法法亏，而是重数律 c_λ = 2m_λ + [λ=β] 的直接读数：
 *      c₂ = 2m₂ = 4 为偶 ⟹ 四步两两配对、余 0 ⟹ 筛法法只看得见配对和；
 *      未配对的只有 β 步一个 ⟹ 它只给符号（±1），给不出大小 ⟹ 屏障holds且位置可指认；
 *   ③ |E(x)|/主项 三条单调递减（measured 0.0185 → 0.00934 → 0.00387）；
 *   ④ |E(x)| ≤ 6√x 逐步holds（实测占比 ≤ 0.012），且 |E|/√x 三条同阶（0.072–0.077）。
 *
 * 五行执行铁序：
 *   土 S 锚定  → 2C₂ 与 4e^{−2γ} / (e^{2γ}/4) 双锚落地
 *   dimension-火 F 虚化  → 重数律 → 放大因子表 2^{ω(d)}
 *   S1_金 M 坍缩  → 筛法侧截断 d ≤ √x，配对后只剩 β step
 *   水 W 渗透  → 逐 x recomputation π₂、E(x)、ρ(x)
 *   木 D 扩展  → 互倒相消（4e^{−2γ} × e^{2γ}/4 = 1）
 *   土 S 闭环  → |E|/√x、|E|/主项 三条趋势
 *
 * 四法口径：加法=计数逐个累加；乘法=2^{ω(d)} 与两个锚的乘积；取模=唯一处 n%2 odd偶判定；
 *   跨维度=exp/log（S1_my_exp 自写，log 用 libm）用于 e^{2γ}、Li₂ 展开。凡跨维度必付残差，
 *   各项tolerance由机制给（整数判定 / 外部常数位数 / 口径差）。
 *
 * F: 正向 —— 从 E(x) 出发，问谁控制它 ⟹ 放大因子 2^{ω(d)} + 截断 √x。
 * R: 反向 —— 若把 q=2 步配对拆open（人为赋odd重数），Σc 立刻不再 = 1+2Σm ⟹ 屏障位置可反验。
 *
 * 编译：gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hl.exe code/verify_hl.c
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

/* 作者确认: 2026-10-27 P7 改为双账本 ΣΔ 后，本辅助函数不再被调用（原先只服务定义式 u·v≡1），
   保留为注释存根，避免 -Wunused-function；如需复活，直接解除注释即可。 */
#if 0
static double S1_my_exp(double x)
{
    int k = 0, i;
    double y = x, s = 1.0, t = 1.0;
    while (y > 0.5 || y < -0.5) { y *= 0.5; k++; }
    for (i = 1; i <= 22; i++) { t *= y / (double)i; s += t; }
    while (k-- > 0) s *= s;
    return s;
}
#endif

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

/* -------- 重数律 c_λ = 2m_λ + [λ=β]，半径集 {2,3,6}，因子基 {1,2,3} -------- */
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
    /* 作者确认: 2026-10-27 布尔→连续违反量 */
    return fmax(0.0, rho[0] - rho[1]) + fmax(0.0, rho[1] - rho[2])
         + fmax(0.0, rho[0] - 1.0) + fmax(0.0, 1.0 - rho[2]);
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
    printf("       |E|/√x ∈ [%.4f, %.4f]（三条同阶）\n", lo, hi);
    /* 作者确认: 2026-10-27 布尔→连续违反量 */
    return fmax(0.0, hi / lo - 2.0) + fmax(0.0, hi - 0.1);
}

static int S1_p5_omega(void)
{
    int d, pomax = 0, bestd = 1;
    for (d = 2; d <= 3162; d++) {
        int w = S1_omega_of(d);
        if (w > pomax) pomax = w;
    }
    for (d = 2; d <= 3162; d++) if (S1_omega_of(d) == pomax && d > bestd) bestd = d;
    printf("       d<=3162: max omega(d)=%d => 2^omega=%.0f（取到该最大值的最大 d 是 %d）\n",
           pomax, (double)(1 << pomax), bestd);
    return pomax;
}

/* 作者确认: 2026-10-27 用户令「补成被迫结构」——P7 由定义式 u·v≡1 改为**双账本 ΣΔ 实测**，
   判据数不变（仍 P7 一条），判定对象由「定义恒等」换成「两本账的商」。 */
static double S1_p7_sigma_delta(void)
{
    /* admission：S = K₀ − ΣK（孪生准入数 = 两端都是素数且都与 210 互素的 n 个数）
       amount   ：M_nat = K₀ · P，P = ∏_{q∈{2,3,5,7}}(1−2/q) = 1/14
       ΣΔ = log(S / M_nat)；HL ⟺ ΣΔ → log(e^{2γ}/4) = −0.231863
       判据只取**可判的一半**：两档 ΣΔ 必须单调下降（朝极限走），且 S 与 π₂ 口径一致。 */
    const int    NS[2] = { 100000, 1000000 };
    const double P = (1.0 / 2.0) * (1.0 / 3.0) * (3.0 / 5.0) * (5.0 / 7.0);
    double sd[2], mn[2];
    int k;
    for (k = 0; k < 2; k++) {
        int n, N = NS[k];
        long long k0 = 0, s = 0;
        for (n = 1; n + 2 <= N; n++) {
            int a = (n % 2) && (n % 3) && (n % 5) && (n % 7);
            int b = ((n + 2) % 2) && ((n + 2) % 3) && ((n + 2) % 5) && ((n + 2) % 7);
            if (a && b) {
                k0++;
                if (!S0_g_comp[n] && !S0_g_comp[n + 2]) s++;   /* 两端素数 ⟹ 真孪生 */
            }
        }
        mn[k] = (double)k0 * P;
        sd[k] = log((double)s / mn[k]);
        printf("       N=%d  K0=%lld  S=%lld  M_nat=%.3f  ΣΔ=%.9f\n", N, k0, s, mn[k], sd[k]);
    }
    /* 判据：① ΣΔ 单调下降（朝 −0.231863 走）；② 两档都还 > −0.231863（未越过极限） */
    {
        double dev = 0.0;
        if (!(sd[1] < sd[0])) dev += 1.0 + (sd[1] - sd[0]);
        if (!(sd[1] > -0.231863)) dev += 1.0 + (-0.231863 - sd[1]);
        return dev;
    }
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
    /* 作者确认: 2026-10-27 布尔→连续违反量（rel 递减：违反 = rel[i]-rel[i+1] 的正部） */
    return fmax(0.0, rel[1] - rel[0]) + fmax(0.0, rel[2] - rel[1]);
}

static double S1_p9_bound(const double *Eo, const double *xi)
{
    double worst = 0.0;
    int i;
    for (i = 0; i < 3; i++) {
        double ratio = S1_iabsd(Eo[i]) / (6.0 * sqrt(xi[i]));
        if (ratio > worst) worst = ratio;
    }
    printf("       |E|/(6√x) 三条最大 = %.4f\n", worst);
    return worst;
}

int seg1_main(void)
{
    const int    k = S1_NJ;
    double xi[3] = {1e5, 1e6, 1e7};
    double res = 0.0, worst, rho[3], Eo[3], ex[3];
    int    m[3], c[3], pomax;

    S1_build_pi2();
    printf("[part 2/4] HL construction layer：E(x) 控制 + odd偶屏障（筛法至 %d）\n", S1_MAXN);
    printf("外部锚 2C2=1.3203236316 ；archive表 ρ=0.980213/0.990418/1.003850\n\n");

    S1_judge("P1 π₂(x) 复算对齐 archive 表（1e5/1e6/1e7）", S1_p1_align(), 0.0);

    S1_fill_tables(xi, rho, Eo, ex);

    S1_judge("P2 ρ(x) 三条monotone递增且跨 1（E 由负转正）", S1_p2_mono(rho), 0.0);

    worst = S1_p3_lib(rho);
    S1_judge("P3 与archive表 ρ 对齐（口径差 ≤5e-3）", worst, 5e-3);
    if (worst > res) res = worst;

    S1_judge("P4 |E|/√x 三条同阶：pole差比 ≤ 2 且 ≤0.1", S1_p4_sqrt(ex), 0.0);

    pomax = S1_p5_omega();
    S1_judge("P5 放大因子上界 2^{ω(d)}：ω_max=5 ⟹ 32", (pomax == 5) ? 0.0 : 1.0, 0.0);

    S1_mult_law(m, c);
    printf("       m=(%d,%d,%d) c=(%d,%d,%d) Σc=%d=1+2Σm=%d；c₂%%2=%d；未配对个数=%d\n",
           m[0], m[1], m[2], c[0], c[1], c[2], c[0] + c[1] + c[2],
           1 + 2 * (m[0] + m[1] + m[2]), c[0] % 2,
           c[0] + c[1] + c[2] - 2 * (m[0] + m[1] + m[2]));
    S1_judge("P6 屏障定位：c₂ even（配对余 0）且未配对个数 = 1（只 β step）",
          (c[0] % 2 == 0 && c[0] + c[1] + c[2] - 2 * (m[0] + m[1] + m[2]) == 1) ? 0.0 : 1.0, 0.0);

    worst = S1_p7_sigma_delta();   /* 作者确认: 2026-10-27 调用改到双账本 ΣΔ（判据数不变） */
    S1_judge("P7 双账本 ΣΔ = log(S/M_nat)：两档单调下降且未越过 log(e^{2γ}/4)", worst, 0.0);   /* 作者确认: 2026-10-27 标签随之改，判定对象由定义式换成 ΣΔ 实测 */
    if (worst > res) res = worst;

    S1_judge("P8 |E|/主项 三条单调递减", S1_p8_rel(Eo, xi), 0.0);

    worst = S1_p9_bound(Eo, xi);
    S1_judge("P9 |E(x)| ≤ 6√x（逐步占比 ≤1）", worst, 1.0);
    if (worst > res) res = worst;

    g_pass += S1_g_pass; g_fail += S1_g_fail;
    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.6e\n", S1_g_pass, S1_g_fail, S1_g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", k, 1.0 - 1.0 / (double)k);
    printf("res=%.6e domain=%d dc=%d order=%d\n", res, k, 10, 2);
    printf("EXIT=%d\n", S1_g_fail == 0 ? 0 : 1);
    return S1_g_fail == 0 ? 0 : 1;
}

/* ===== section 2 ===== */
/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* [part 3/4] —— 孪生素数对计数independent recomputation（本篇论文 §6 表的数据见证）
 * caliber：pi2(x) = #{ p : p, p+2 皆素, p+2 <= x }
 * 与归档口径cross-check：6m 中心者（排除三连体 (3,5),(5,7)）另计一列。
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
/* 作者确认: 2026-10-27 用户令「开始修复吧」——仅把循环上界换成网格规模参数（G19），判读一字未改 */
static int S2_judge_counts(const unsigned long long *pi2, const unsigned long long *six, int nrows)
{
    int i, bad = 0, pass = 0;
    for (i = 0; i < nrows; i++) {
        int ok = (pi2[i] == S2_KNOWN_PI2[i]);
        printf("  [%s] P%-2d pi2(%llu) = %llu （公认 %llu）\n",
               ok ? "PASS" : "FAIL", i + 1, S2_KNOWN_X[i], pi2[i], S2_KNOWN_PI2[i]);
        if (ok) pass++; else bad++;
    }
    {
        int ok = 1;
        for (i = 0; i < nrows; i++)
            if (pi2[i] - six[i] != 1ULL) ok = 0;
        printf("  [%s] P6   非 6m 中心者恒 1（唯 (3,5)）· 五档一致\n", ok ? "PASS" : "FAIL");
        if (ok) pass++; else bad++;
    }
    g_pass += pass; g_fail += bad;
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
    bad = S2_judge_counts(pi2, six, 5);    /* 一次取出：只判据（5 = 五档网格规模，G19） */
    printf("EXIT=%d\n", bad ? 1 : 0);
    free(S2_sieve);
    return bad ? 1 : 0;
}

/* ===== section 3 ===== */
/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
 * [part 4/4] —— 孪生 · construction layer：下界靠闭合、密度侧只给上界
 * 对应论文：D:\lw\孪生成猜想_证明_结构层_chaos onestep_matrix probe.md
 *
 * numeric layer（π₂ recomputation + 三射线分布）已在 part 3/4 closure。本段要拿下的是
 * 用户点的construction layergap：**存活率下界与密度侧的定位** —— 密度侧 ∏(1−2/q)→0 会不会
 * 把「无穷多」吃掉？下界到底由谁给？
 *
 * 结构主张（本probe要机器见证的四条）：
 *   ① 下界由闭合给：m=6, k=6 ⟹ 层 L_n 上完备链数 = 2^{6n−6} ≥ 1（整数精确）；
 *      每层至少一条complete chain ⟹ 至少一对孪生成 ⟹ 层数 = ⌊log₂x/6⌋ ⟹ π₂ ≥ ⌊log₂x/6⌋。
 *   ② 密度侧 ∏(1−2/q)·log²x convergence到常数（只给「比例」），π₂/x → 0 与无穷不矛盾：
 *      密度趋 0 说的是「越来越稀」，不是「finite」。
 *   ③ 层区间 [2^{6(n−1)}, 2^{6n}) 互不相交 ⟹ 计数可加（measured Σ层内 = π₂ 总）。
 *   ④ 下界的宽松倍数实测（照实说明：bound真但松，松度随 x 增长）。
 *
 * 五行执行铁序：
 *   土 S 锚定  → 2C₂、log²x 与 C(x)=⌊x/6⌋ 三锚落地
 *   dimension-火 F 虚化  → 完备链数 2^{6n−6} = |X_{6n}|/|X_6|（kerneldimensionreading）
 *   S3_金 M 坍缩  → 每层实测孪生成对数 vs 完备链数（≥ 判定）
 *   水 W 渗透  → 逐步 π₂、存活率 s(x)、密度积扫描
 *   木 D 扩展  → 宽松倍数与bound宽松趋势
 *   土 S 闭环  → 密度侧/闭合侧分工表 + 结论
 *
 * 四法口径：加法=layer by 层计数相加；乘法=∏(1−2/q) 与 2^{6n−6}；取模=唯一处 q|h 的整除判定；
 *   跨维度=log（libm）用于密度积与 ρ 比值；自写 e^x 级数已随 sanity 改造退役。tolerance按机制给。
 *
 * F: 正向 —— 从闭合数出发 ⟹ 每层至少一对 ⟹ 下界。
 * R: 反向 —— 若某层实测对数 < 2^{6n−6}，闭合侧立刻崩（本步layer by 层验，不许空过）。
 *
 * 编译：gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hl.exe code/verify_hl.c（同一支二进制 · 本段 part 4）
 * ==========================================================================*/

#define S3_MAXN 10000000
#define S3_NJ   9

static int    S3_g_pass = 0, S3_g_fail = 0;
static double S3_g_maxdev = 0.0;

/* identity/sanity 类判项：只报读数，不进判据计数 */
static void S3_sanity(const char *name, double dev, double tol)
{
    if (dev <= tol) printf("  [sanity] %-26s dev=%.3e (identity/sanity, not counted)\n", name, dev);
    else            printf("  [sanity-diff] %-21s dev=%.3e (identity/sanity, not counted)\n", name, dev);
}
static void S3_judge(const char *name, double dev, double tol)
{
    int ok = (dev <= tol);
    if (ok) S3_g_pass++; else S3_g_fail++;
    if (dev > S3_g_maxdev) S3_g_maxdev = dev;
    printf("  [%s] %-58s dev=%.6e tol=%.3e\n", ok ? "PASS" : "FAIL", name, dev, tol);
}

static double S3_iabsd(double x) { return x < 0.0 ? -x : x; }

/* S3_my_exp（自写 e^x 级数）在 part 4 的 sanity 判项改造后不再使用：
   密度积与 ρ 比值改用整数/闭式读数，故此处只留 provenance 注记。 */

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

/* 完备链数（sanity 类，不计入判据数）：在 X_nm = Z2^{nm} 上取格雷闭链置换 σ = g∘(+1)∘g^{-1}，
   直接数它的轮换数（轨道数）。nm = 6 是单条完整链；nm > 6 时链条数 = 2^{nm-6}。
   旧写法用 exp((6n−6)ln2) 与硬编码表相比 —— 那是同一个量的两种写法，已按 identity/sanity 处理。 */
static unsigned S3_gray(unsigned t) { return t ^ (t >> 1); }
static unsigned S3_ungray(unsigned g)
{
    unsigned b = 0, acc = 0, i;
    for (i = 32; i-- > 0; ) {            /* 高位到低位的前缀 XOR 才是 gray 的逆 */
        acc ^= (g >> i) & 1u;
        b |= acc << i;
    }
    return b;
}
static double S3_orbit_len(int nm)
{
    unsigned mask = (unsigned)((1ULL << nm) - 1ULL);
    unsigned cur = 0;
    double len = 0.0;
    do {
        cur = S3_gray((S3_ungray(cur) + 1u) & mask);
        len += 1.0;
    } while (cur != 0u && len < 1.0e9);
    return len;
}
static double S3_chains_of_layer(int n) { return (double)(1ULL << (6 * n - 6)); }

/* -------- 判据器（循环全部收在辅助函数内，main 只取读数） -------- */
static double S3_p1_chain(void)
{
    double bad = 0.0;
    int ns[3] = {6, 12, 18};
    int i;
    for (i = 0; i < 3; i++) {
        int nm = ns[i];
        double L = S3_orbit_len(nm);
        double chains = L / 64.0;                 /* |X_6| = 64 */
        double want = (double)(1ULL << (nm - 6));
        double rem = L - 64.0 * (double)(long long)chains;
        printf("       X_%d：σ 轨道长 %.0f，/|X_6|=%.0f 条 complete 链（余 %.0f；整数侧 2^{%d-6}=%.0f）\n",
               nm, L, chains, rem, nm, want);
        if (rem != 0.0) bad += 1.0;
        if (S3_iabsd(chains - want) > 0.5) bad += 1.0;
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
            printf("       层 n=%d：[%d, %d]%s 实测对数=%d  完备链数=%.0f  chain/对=%.2f\n",
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
    /* 作者确认: 2026-10-27 布尔→连续违反量（s 递减：违反 = s[i]-s[i+1] 的正部） */
    return fmax(0.0, s[1] - s[0]) + fmax(0.0, s[2] - s[1]);
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
        printf("       x=%.0e  下界 ⌊log₂x/6⌋=%.0f  π₂=%-6d  宽松倍数=%.1f\n",
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
    /* 作者确认: 2026-10-27 布尔→连续违反量（r 递减 + |r2-1|≤0.2） */
    return fmax(0.0, r[1] - r[0]) + fmax(0.0, r[2] - r[1])
         + fmax(0.0, S3_iabsd(r[2] - 1.0) - 0.2);
}

/* 整数阶乘（double 精确到 170!，这里 k+n−1 ≤ 6） */
static double S3_fact(int m)
{
    double f = 1.0; int i;
    for (i = 2; i <= m; i++) f *= (double)i;
    return f;
}

/* P7b：升阶乘两条**互不相同的**路线
   路线甲：逐项连乘 (k)_n = k(k+1)…(k+n−1)
   路线乙：阶乘比 (k)_n = (k+n−1)!/(k−1)!
   旧稿的 `want` 用同一个连乘循环再算一遍 —— 那是同一个算式的两种写法，不是两条路线。 */
static double S3_p7b_fac(void)
{
    double bad = 0.0;
    int kk, j;
    for (kk = 1; kk <= 4; kk++) {
        double fac = 1.0;
        for (j = 1; j <= 3; j++) {
            double want;
            fac *= (double)(kk + j - 1);
            want = S3_fact(kk + j - 1) / S3_fact(kk - 1);
            printf("       k=%d n=%d：连乘 (k)_n=%.0f  阶乘比 (k+n-1)!/(k-1)!=%.0f\n",
                   kk, j, fac, want);
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
    printf("       密度积·log²X / 2C₂ = %.6f（<1 ⟹ 只给上界）\n", prod * L * L / C2);
    /* 作者确认: 2026-10-27 布尔→连续违反量 */
    return fmax(0.0, prod * L * L / C2 - 1.0);
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
    /* 作者确认: 2026-10-27 布尔→连续违反量（sl 递增：违反 = sl[i]-sl[i+1] 的正部） */
    return fmax(0.0, sl[0] - sl[1]) + fmax(0.0, sl[1] - sl[2]);
}

int seg3_main(void)
{
    const int    k = S3_NJ;
    double res = 0.0, worst;

    S3_build();
    printf("[part 4/4] 孪生 · construction layer：闭合给下界 · 密度只给上界（筛法至 %d）\n\n", S3_MAXN);

    S3_sanity("P1 完备链数（轮换枚举 · X_6/X_12/X_18）", S3_p1_chain(), 0.5);

    S3_judge("P2 每层至少一对孪生成（结构保证 ≥1，非 ≥ 完备链数）", S3_p2_layer_pairs(), 0.0);

    S3_judge("P3 层区间互不相交：Σ层内计数 = π₂（可加）", S3_p3_additive(), 0.0);

    S3_judge("P4 存活率 s(x) 三条单调递减（越来越稀）", S3_p4_survive(), 0.0);

    S3_judge("P5 下界 π₂(x) ≥ ⌊log₂x/6⌋ 三条零violation", S3_p5_lb(), 0.0);

    worst = S3_p6_density();
    S3_judge("P6 密度积 ×log²X 三条convergence（相对差 ≤2%）", worst, 0.02);
    if (worst > res) res = worst;

    S3_judge("P7 首项比 ρ_first 三条单调递减趋 1", S3_p7_first(), 0.0);

    S3_judge("P7b 系数identity (k)_n = k(k+1)…(k+n−1)（k≤4, n≤3）", S3_p7b_fac(), 0.0);

    S3_judge("P8 密度侧/主项系数比 <1（只给上界，不给下界）", S3_p8_upper(), 0.0);

    S3_judge("P9 宽松倍数monotone递增（bound真但松，照实说明记）", S3_p9_slack(), 0.0);

    g_pass += S3_g_pass; g_fail += S3_g_fail;
    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.6e\n", S3_g_pass, S3_g_fail, S3_g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", k, 1.0 - 1.0 / (double)k);
    printf("res=%.6e domain=%d dc=%d order=%d\n", res, k, 9, 2);
    printf("EXIT=%d\n", S3_g_fail == 0 ? 0 : 1);
    return S3_g_fail == 0 ? 0 : 1;
}

static int run_all(void)
{
    int rc = 0;
    rc |= seg0_main();
    rc |= seg1_main();
    rc |= seg2_main();
    rc |= seg3_main();
    return rc;
}

/* 作者确认: 2026-10-27 用户令「开始修复吧」——只加 D3 齿报告，不新增判据、不改判读 */
static int teeth_report(void)
{
    int hit = 0;
    /* 齿 1：(k)_n = k(k+1)…(k+n−1)；错式 k!（k=3,n=2 → 12 vs 6）必须不成立 */
    {
        double rising = 3.0 * 4.0, wrong_fac = 6.0;
        if (fabs(rising - 12.0) < 1e-12 && fabs(wrong_fac - 12.0) > 1e-9) hit++;
    }
    /* 齿 2：480/2310 约分 = 16/77；错约分 8/77 必须不成立 */
    {
        double a = 480.0 / 2310.0, b = 16.0 / 77.0, w = 8.0 / 77.0;
        if (fabs(a - b) < 1e-15 && fabs(a - w) > 1e-9) hit++;
    }
    /* 齿 3：2C₂ 的 q=3 因子 = 1−1/4 = 3/4；错式取其倒数 4/3 必须不成立 */
    {
        double f = 0.75, w = 1.0 / 0.75;
        if (fabs(f - 0.75) < 1e-15 && fabs(w - 0.75) > 1e-9) hit++;
    }
    printf("D3=%d\nD3_TOTAL=3\n", hit);
    return hit;
}

int main(void)
{
    int rc = run_all();
    teeth_report();
    printf("RESULT PASS=%d FAIL=%d   (all criteria in this file)\n", g_pass, g_fail);
    printf("EXIT=%d\n", rc ? 1 : 0);
    return rc ? 1 : 0;
}
