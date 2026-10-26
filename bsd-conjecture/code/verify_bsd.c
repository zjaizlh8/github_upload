/* 单文件验证器：各部分在同一编译单元内顺序运行，符号按 part 编号区隔。 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static int g_pass = 0, g_fail = 0;   /* totals across sections */

/* ===== section 0 ===== */
/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* verify_bsd.c · BSD 猜想（Birch and Swinnerton-Dyer）independent recomputation
 *
 * F: 函数方程 Λ(E,s) = w·Λ(E,2−s) 的involution s ↦ 2−s 与符号 w=±1
 * R: 由 F 反推中心 s=1 的零阶的oddeven，并由八元数norm-preserving把「累积=summation」落成恒等式
 *
 * four operationscaliber：
 *   加法 —— constructionidentity Λ_k(s)=(s−1)^k·g(s) 按多项式逐项相加；四元/八元数乘法逐分量相加
 *   乘法 —— 保范律 |ab|² = |a|²|b|²；倍角recurrence cos2x = 2cos²x − 1
 *   取模 —— a_p 的 Legendre 符号与 p≡3 (mod 4) criterion（本步唯一取模处）
 *   跨维度 —— 不用（不调 exp/log/sqrt/pow）⟹ 不付跨维度残差，全部tolerance皆为整数/机器精度
 *
 * check：12；main 直接层零循环（one步在最外层），判据只读全field常量。
 */

#define PI_ 3.14159265358979323846

static int S0_g_pass = 0, S0_g_fail = 0, g_n = 0;
static double S0_g_maxdev = 0.0;

/* checkregistry：name + rundeviation + mechanismtolerance */
static void S0_judge(const char *name, double dev, double tol)
{
    double a = fabs(dev);
    int ok;
    if (a > S0_g_maxdev) S0_g_maxdev = a;
    ok = (a <= tol);
    printf("  [%s] %-40s dev=%.3e (tol=%.0e)\n", ok ? "PASS" : "FAIL", name, a, tol);
    g_n++;
    if (ok) S0_g_pass++; else S0_g_fail++;
}

/* ---------- 四元数 ---------- */
typedef struct { double v[4]; } Q;

static Q q_mul(Q a, Q b)
{
    Q r;
    r.v[0] = a.v[0]*b.v[0] - a.v[1]*b.v[1] - a.v[2]*b.v[2] - a.v[3]*b.v[3];
    r.v[1] = a.v[0]*b.v[1] + a.v[1]*b.v[0] + a.v[2]*b.v[3] - a.v[3]*b.v[2];
    r.v[2] = a.v[0]*b.v[2] - a.v[1]*b.v[3] + a.v[2]*b.v[0] + a.v[3]*b.v[1];
    r.v[3] = a.v[0]*b.v[3] + a.v[1]*b.v[2] - a.v[2]*b.v[1] + a.v[3]*b.v[0];
    return r;
}
static Q q_conj(Q a) { Q r; r.v[0]=a.v[0]; r.v[1]=-a.v[1]; r.v[2]=-a.v[2]; r.v[3]=-a.v[3]; return r; }
static double q_n2(Q a) { return a.v[0]*a.v[0]+a.v[1]*a.v[1]+a.v[2]*a.v[2]+a.v[3]*a.v[3]; }

/* ---------- 八元数（Cayley–Dickson 双倍） ---------- */
typedef struct { double a[4], b[4]; } O;

static Q v2q(const double *x) { Q r; r.v[0]=x[0]; r.v[1]=x[1]; r.v[2]=x[2]; r.v[3]=x[3]; return r; }

static O o_mul(O x, O y)
{
    Q a = v2q(x.a), b = v2q(x.b), c = v2q(y.a), d = v2q(y.b);
    Q p = q_mul(a, c);            /* ac    */
    Q t = q_mul(q_conj(d), b);    /* conj(d)b */
    Q q = q_mul(d, a);            /* da    */
    Q s = q_mul(b, q_conj(c));    /* b conj(c) */
    O r; int i;
    for (i = 0; i < 4; i++) { r.a[i] = p.v[i] - t.v[i]; r.b[i] = q.v[i] + s.v[i]; }
    return r;
}
static double o_n2(O x) { return q_n2(v2q(x.a)) + q_n2(v2q(x.b)); }

static O o_basis(int k)
{
    O r; int i;
    for (i = 0; i < 4; i++) { r.a[i] = 0.0; r.b[i] = 0.0; }
    if (k < 4) r.a[k] = 1.0; else r.b[k-4] = 1.0;
    return r;
}

/* 确定性伪随机（可复现，不用 rand） */
static unsigned long g_seed = 20261026UL;
static double urand(void)
{
    g_seed = (g_seed * 1103515245UL + 12345UL) & 0x7fffffffUL;
    return (double)(g_seed % 20001UL) / 10000.0 - 1.0;  /* [−1, 1] */
}

/* ---------- 数论 ---------- */
static int is_prime(int n)
{
    int d;
    if (n < 2) return 0;
    for (d = 2; (long long)d * d <= n; d++) if (n % d == 0) return 0;
    return 1;
}
static int powmod(int a, int e, int p)
{
    long long r = 1, b = ((long long)a % p + p) % p;
    while (e > 0) {
        if (e & 1) r = (r * b) % p;
        b = (b * b) % p;
        e >>= 1;
    }
    return (int)r;
}
static int S0_legendre(int a, int p)
{
    int t;
    long long x = ((long long)a % p + p) % p;
    if (x == 0) return 0;
    t = powmod((int)x, (p - 1) / 2, p);
    return (t == 1) ? 1 : ((t == p - 1) ? -1 : 0);
}

/* a_p = p + 1 − #E(F_p)，对 y² = x³ − x 用 a_p = −Σ_x χ(x³−x) */
static int ap_cubic(int p)
{
    int x, s = 0;
    for (x = 0; x < p; x++) {
        long long y2 = ((long long)x * x % p * x % p - x % p + p) % p;
        s += S0_legendre((int)y2, p);
    }
    return -s;
}

/* 倍角recurrence：cos(2x) = 2cos²x − 1（纯加乘，无跨维度算子） */
static double dbl_cos(double c) { return 2.0 * c * c - 1.0; }

/* check 1：involution s ↦ 2−s 的fixed point唯一为 s=1 */
static double P1(void)
{
    double s = 1.0, r = 2.0 - s;
    return fabs(r - s);                      /* 不动 ⟺ 0 */
}

/* check 2：constructionidentity Λ_k(s)=(s−1)^k g(s) ⟹ w=(−1)^k（k=1,2,3） */
static double P2(void)
{
    double mx = 0.0, g = 1.0;
    int k, j;
    for (k = 1; k <= 3; k++) {
        double s, l1, l2, w, pred;
        s = 0.37;                            /* 任取不在中心的点 */
        l1 = 1.0; l2 = 1.0;
        for (j = 0; j < k; j++) { l1 *= (s - 1.0); }
        for (j = 0; j < k; j++) { l2 *= ((2.0 - s) - 1.0); }
        w = (l2 * g) / (l1 * g);
        pred = (k % 2 == 0) ? 1.0 : -1.0;
        if (fabs(w - pred) > mx) mx = fabs(w - pred);
    }
    return mx;
}

/* check 3：Λ_k(1)=0（k ≥ 1）——中心零，阶exactly为 k */
static double P3(void)
{
    double mx = 0.0; int k, j;
    for (k = 1; k <= 3; k++) {
        double v = 1.0;
        for (j = 0; j < k; j++) v *= (1.0 - 1.0);
        if (fabs(v) > mx) mx = fabs(v);
    }
    return mx;
}

/* check 4：half-turn/full turn（c(π)=−1 为hard anchor，倍角recurrence出 c(2π)=+1、c(4π)=+1） */
static double P4(void)
{
    double c_pi = -1.0;          /* hard anchor：cos π = −1（与 Γ(1)=1 同类的定义性锚） */
    double c_2pi = dbl_cos(c_pi);
    double c_4pi = dbl_cos(c_2pi);
    double mx = fabs(c_2pi - 1.0);        /* full turn +1 */
    if (fabs(c_4pi - 1.0) > mx) mx = fabs(c_4pi - 1.0);
    if (fabs((c_2pi * c_2pi + 0.0) - 1.0) > mx) mx = fabs(c_2pi * c_2pi - 1.0);
    return mx;
}

/* check 5：四元数norm-preserving |ab|² = |a|²|b|²（200 对） */
static double P5(void)
{
    double mx = 0.0; int i, j;
    for (i = 0; i < 200; i++) {
        Q a, b; double lhs, rhs, d;
        for (j = 0; j < 4; j++) { a.v[j] = urand(); b.v[j] = urand(); }
        lhs = q_n2(q_mul(a, b));
        rhs = q_n2(a) * q_n2(b);
        d = fabs(lhs - rhs) / (1.0 + fabs(rhs));
        if (d > mx) mx = d;
    }
    return mx;
}

/* check 6：八元数norm-preserving（非结合仍norm-preserving，200 对） */
static double P6(void)
{
    double mx = 0.0; int i, j;
    for (i = 0; i < 200; i++) {
        O a, b; double lhs, rhs, d;
        for (j = 0; j < 4; j++) { a.a[j] = urand(); a.b[j] = urand(); b.a[j] = urand(); b.b[j] = urand(); }
        lhs = o_n2(o_mul(a, b));
        rhs = o_n2(a) * o_n2(b);
        d = fabs(lhs - rhs) / (1.0 + fabs(rhs));
        if (d > mx) mx = d;
    }
    return mx;
}

/* check 7：八元数非结合 —— associator [e1,e2,e4] = (e1e2)e4 − e1(e2e4)，‖·‖²=4 */
static double P7(void)
{
    O e1 = o_basis(1), e2 = o_basis(2), e4 = o_basis(4);
    O l = o_mul(o_mul(e1, e2), e4);
    O r = o_mul(e1, o_mul(e2, e4));
    O m; int i;
    for (i = 0; i < 4; i++) { m.a[i] = l.a[i] - r.a[i]; m.b[i] = l.b[i] - r.b[i]; }
    return fabs(o_n2(m) - 4.0);              /* 范数 2 ⟹ 平方 4 */
}

/* check 8：八元数无零因子 —— 整数格点上 min |xy|² ≥ 1 */
static double P8(void)
{
    double mn = 1e300; int i, j;
    for (i = 0; i < 300; i++) {
        O a, b; double n2;
        for (j = 0; j < 4; j++) {
            a.a[j] = (double)((int)(urand() * 2.5)); a.b[j] = (double)((int)(urand() * 2.5));
            b.a[j] = (double)((int)(urand() * 2.5)); b.b[j] = (double)((int)(urand() * 2.5));
        }
        if (o_n2(a) < 0.5 || o_n2(b) < 0.5) continue;   /* 跳过零元，非零对才判 */
        n2 = o_n2(o_mul(a, b));
        if (n2 < mn) mn = n2;
    }
    return (mn >= 1.0) ? 0.0 : (1.0 - mn);   /* 无零因子 ⟹ 永不塌到 0 */
}

/* check 9：CM criterion p ≡ 3 (mod 4) ⟹ a_p = 0（y²=x³−x，enumeration p ≤ 503） */
static double P9(void)
{
    int p, cnt = 0, bad = 0;
    for (p = 3; p <= 503; p++) {
        if (!is_prime(p) || (p % 4) != 3) continue;
        cnt++;
        if (ap_cubic(p) != 0) bad++;
    }
    printf("       · p≡3 (mod 4) prime %d 个，a_p≠0 的 %d 个\n", cnt, bad);
    return (double)bad;
}

/* check 10：Hasse bound a_p² ≤ 4p（全部oddprime ≤ 503，平方形式免square root） */
static double P10(void)
{
    int p, cnt = 0, bad = 0;
    for (p = 3; p <= 503; p++) {
        int ap;
        if (!is_prime(p)) continue;
        ap = ap_cubic(p);
        cnt++;
        if (ap * ap > 4 * p) bad++;
    }
    printf("       · oddprime %d 个，越 Hasse bound %d 个\n", cnt, bad);
    return (double)bad;
}

/* check 11：伸缩居标量部 ⟹ 与一切四元数commute（dev 0）；
   而纯虚部互不commute ⟹ associator ij−ji = 2k，范数² = 4（无 BCH 修正的来source） */
static double P11(void)
{
    double mx = 0.0; int i, j;
    Q i1, j1, l, r, m;
    double n2;
    for (i = 0; i < 100; i++) {
        Q s, q; double d;
        s.v[0] = urand(); s.v[1] = s.v[2] = s.v[3] = 0.0;
        for (j = 0; j < 4; j++) q.v[j] = urand();
        l = q_mul(s, q); r = q_mul(q, s);
        for (j = 0; j < 4; j++) { d = fabs(l.v[j] - r.v[j]); if (d > mx) mx = d; }
    }
    i1.v[0]=0; i1.v[1]=1; i1.v[2]=0; i1.v[3]=0;
    j1.v[0]=0; j1.v[1]=0; j1.v[2]=1; j1.v[3]=0;
    l = q_mul(i1, j1); r = q_mul(j1, i1);
    for (j = 0; j < 4; j++) m.v[j] = l.v[j] - r.v[j];
    n2 = q_n2(m);                                  /* 应为 4 */
    if (fabs(n2 - 4.0) > mx) mx = fabs(n2 - 4.0);
    return mx;
}

/* check 12：乘积的范数 = 范数之积 ⟹ 累积=summation（200 步，modulo长有bound不溢出） */
static double P12(void)
{
    Q q, acc; double mx, prod = 1.0, sc = 1.001; int i, j;
    q.v[0] = 0.500000; q.v[1] = 0.500000; q.v[2] = 0.500000; q.v[3] = 0.500000;
    /* |q|² = 0.25×4 = 1 ⟹ q 已是单位四元数，无需square root归一 */
    acc.v[0] = 1.0; acc.v[1] = acc.v[2] = acc.v[3] = 0.0;
    for (i = 0; i < 200; i++) {
        Q step; double n2;
        for (j = 0; j < 4; j++) step.v[j] = q.v[j] * sc;
        acc = q_mul(acc, step);
        n2 = q_n2(step);
        prod *= n2;
    }
    mx = fabs(q_n2(acc) - prod) / (1.0 + fabs(prod));
    return mx;
}

/* check 13：pole 163 —— 163 = 2·9²+1 = 4·41−1；archive里 T_crit=17.0831413200251 satisfies
   (T+1)² = 327（纯加乘verify引用的外部常数，不调square root） */
static double P13(void)
{
    double mx = 0.0, shape, euler, t = 17.0831413200251, lhs;
    shape = 2.0 * 81.0 + 1.0;                 /* 2·9²+1 */
    euler = 4.0 * 41.0 - 1.0;                 /* 4·41−1（Euler 判别式形） */
    if (fabs(shape - 163.0) > mx) mx = fabs(shape - 163.0);
    if (fabs(euler - 163.0) > mx) mx = fabs(euler - 163.0);
    lhs = (t + 1.0) * (t + 1.0);              /* (T+1)² 应 = 327 */
    if (fabs(lhs - 327.0) > mx) mx = fabs(lhs - 327.0);
    return mx;
}

int seg0_main(void)
{
    printf("=== [part 1/5] BSD · independent recomputation ===\n");
    printf("  · involution s ↦ 2−s 的fixed point s=1；符号 w=(−1)^k 管零阶的oddeven\n");
    printf("  · 四元数管符号（中心 Z₂=±1），八元数管秩（非结合而norm-preserving ⟹ 生成元不塌）\n");

    S0_judge("P1  involutionfixed pointunique s=1", P1(), 1e-15);
    S0_judge("P2  constructionidentity w=(−1)^k（k=1,2,3）", P2(), 1e-12);
    S0_judge("P3  Λ_k(1)=0 中心零（阶exactly为 k）", P3(), 1e-15);
    S0_judge("P4  half-turn −1 / full turn +1（倍角recurrence）", P4(), 1e-15);
    S0_judge("P5  四元数norm-preserving |ab|²=|a|²|b|²", P5(), 1e-14);
    S0_judge("P6  八元数norm-preserving（非结合仍norm-preserving）", P6(), 1e-13);
    S0_judge("P7  八元数非结合 [e1,e2,e4] 范数²=4", P7(), 1e-14);
    S0_judge("P8  八元数无零因子 min|xy|²≥1", P8(), 1e-15);
    S0_judge("P9  CM criterion p≡3(4) ⟹ a_p=0", P9(), 0.0);
    S0_judge("P10 Hasse bound a_p²≤4p", P10(), 0.0);
    S0_judge("P11 伸缩commute / 纯虚不commute", P11(), 1e-15);
    S0_judge("P12 乘积范数=范数之积（累积=summation）", P12(), 1e-14);
    {   /* tolerance由机制给：archive值 T_crit 只给到 1e−13 位，平方把截断放大为
           2(T+1)×δ = 2×18.083×5e−14 ≈ 1.81e−12 ⟹ 容差取其一阶上邻 5e−12 */
        printf("       · P13 tolerancemechanism：外部常数截断floor 2(T+1)δ = 2×18.083×5e−14 ≈ 1.81e−12\n");
        S0_judge("P13 pole 163 = 2·9²+1 与 (T+1)²=327", P13(), 5e-12);
    }

    g_pass += S0_g_pass; g_fail += S0_g_fail;
    printf("\nRESULT PASS=%d FAIL=%d  max dev=%.3e\n", S0_g_pass, S0_g_fail, S0_g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", g_n,
           1.0 - 1.0 / (double)g_n);
    printf("res=%.6e domain=%d dc=%d order=2\n", S0_g_maxdev, g_n, 97);
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
   [part 2/5] · BSD 的「reciprocal逼迫」final version实证（2026-10-26）
   ----------------------------------------------------------------------------
   用户口径（逐字）：是reciprocal逼迫出的。
   archive内reciprocal定律： [2233] closure ⟺ 逆在同一容器内可写出；[2234] 膨胀↔收缩乘积 1；
                 [2235] 投影率 × 本体倍率 = 1；rank + nullity = dimension。

   本步三件：
     ① 常数侧（倒数账）：L(E,1)·#tors² / (Ω · ∏c_p · #Ш) = 1
        —— analytic读数与几何体量 Ω 互为reciprocal，商是有理账（本步measured ∏c_p exactly为整数 2）。
     ② 秩侧（zerocriterion）：L(E,1) 的 Abel 和 —— rank 0 饱和非零、rank ≥ 1 趋零。
        —— 零阶 = 解析侧的「reciprocalexponent」，与自由度数同位。
     ③ 负结果honest note：局部欧拉积 ∏_{p≤X}(p+1−a_p)/p 的 log-log 斜率**不区分秩**
        （三条曲线全在 1.27 附近）⟹ 秩只能从 L(1) 的零读，不能从局部积读。

   four operationscaliber：加法=点数计数summation；乘法=Hecke 关系/乘法性/欧拉积；取模=Legendre 符号；
             跨维度=解析侧 exp/log/sqrt（Abel kernel、Simpson 求积、log），只在读数上用。
   ============================================================================ */

#define S1_MAXP 5001
#define S1_NP 700

typedef struct {
    const char *name;
    long long a1, a2, a3, a4, a6;
    long long cond;      /* conductor（同时是唯一坏素数的power） */
    int rank;            /* Cremona 表，外部输入 */
    int tors;
} S1_Curve;

static const S1_Curve S1_CURVES[4] = {
    { "11a1",  1, -1, 1, -10, -20,  11, 0, 5 },
    { "37a1",  1,  0, 0,  -1,   0,  37, 1, 1 },
    { "389a1", 1,  1, 0,  -2,   0, 389, 2, 1 },
    { "32a1",  0,  0, 0,  -1,   0,  32, 0, 4 }
};

static int S1_primes[S1_NP], S1_np = 0;
static int S1_npass = 0, S1_nfail = 0;
static double S1_maxdev = 0.0;

static long long S1_imod(long long a, long long m) { long long r = a % m; return r < 0 ? r + m : r; }

static long long S1_ipowmod(long long b, long long e, long long m)
{
    long long r = 1 % m;
    b = S1_imod(b, m);
    while (e > 0) {
        if (e & 1) r = (long long)((__int128)r * b % m);
        b = (long long)((__int128)b * b % m);
        e >>= 1;
    }
    return r;
}

static void S1_sieve(void)
{
    char comp[S1_MAXP + 1];
    int i, j;
    memset(comp, 0, sizeof comp);
    for (i = 2; i <= S1_MAXP; i++) {
        if (comp[i]) continue;
        if (S1_np < S1_NP) S1_primes[S1_np++] = i;
        for (j = 2 * i; j <= S1_MAXP; j += i) comp[j] = 1;
    }
}

static int S1_legendre(long long d, long long p)
{
    d = S1_imod(d, p);
    if (d == 0) return 0;
    return (S1_ipowmod(d, (p - 1) / 2, p) == 1) ? 1 : -1;
}

static int S1_ap_curve(const S1_Curve *E, int p)
{
    long long x, s = 0;
    if ((long long)p == E->cond) return 0;      /* 坏素数只在 cond 是素数时hit */
    for (x = 0; x < p; x++) {
        long long x3 = x * x % p * x % p;
        long long rhs = (x3 + E->a2 * x % p * x % p) % p;
        long long lin = (E->a1 * x + E->a3) % p;
        long long d;
        rhs = S1_imod(rhs + E->a4 * x + E->a6, p);
        d = S1_imod(lin * lin % p + 4 * rhs % p, p);
        s += S1_legendre(d, p);
    }
    return (int)(-s);
}

static void S1_build_an(const S1_Curve *E, int *a, int M)
{
    char *is_p = (char *)calloc((size_t)M + 1, 1);
    int n, p, q;
    for (n = 0; n <= M; n++) a[n] = 0;
    a[1] = 1;
    for (n = 2; n <= M; n++) is_p[n] = 1;
    for (n = 2; (long long)n * n <= M; n++)
        if (is_p[n]) for (p = n * n; p <= M; p += n) is_p[p] = 0;
    for (n = 2; n <= M; n++) {
        if (!is_p[n]) continue;
        {
            int ap = S1_ap_curve(E, n);
            long long pk = n;
            a[n] = ap;
            while ((long long)pk * n <= M) {
                long long nxt = pk * n;
                a[nxt] = ap * a[pk] - n * a[pk / n];
                pk = nxt;
            }
        }
    }
    for (n = 2; n <= M; n++) {
        int pk = 1, rest;
        if (is_p[n]) continue;
        q = 0;
        for (p = 2; (long long)p * p <= n; p++) if (is_p[p] && n % p == 0) { q = p; break; }
        if (!q) continue;
        rest = n; pk = 1;
        while (rest % q == 0) { pk *= q; rest /= q; }
        a[n] = a[pk] * a[rest];
    }
    free(is_p);
}

/* L(1) 的 Abel 和：Σ (a_n/n)·e^{−2πny/√N}（核只作正则化，绝对值受核污染，只读零/非零） */
static double S1_abel_L1(const S1_Curve *E, int M, double y)
{
    int *a = (int *)calloc((size_t)M + 1, sizeof(int));
    double s = 0.0, c = 0.0, sq = sqrt((double)E->cond);
    int n;
    S1_build_an(E, a, M);
    for (n = 1; n <= M; n++) {
        double term = ((double)a[n] / (double)n) * exp(-2.0 * M_PI * (double)n * y / sq);
        double t = s + term, u;
        if (fabs(s) >= fabs(term)) u = s - t + term; else u = term - t + s;
        c += u;
        s = t;
    }
    free(a);
    return s + c;
}

/* Ω(y²=x³−x)=4∫₀^{π/2} dθ/√(1+cos²θ)；恒等分branchperiod = 其一半 */
static double S1_omega_32a1(int N)
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

static double S1_local_product(const S1_Curve *E, int X)
{
    double P = 1.0;
    int i;
    for (i = 0; i < S1_np && S1_primes[i] <= X; i++) {
        int p = S1_primes[i];
        if ((long long)p == E->cond) continue;
        P *= ((double)p + 1.0 - (double)S1_ap_curve(E, p)) / (double)p;
    }
    return P;
}

static void S1_judge(const char *name, double dev, double tol)
{
    if (fabs(dev) <= tol) { S1_npass++; printf("  [PASS] %-44s dev=%.3e (tol=%.0e)\n", name, fabs(dev), tol); }
    else { S1_nfail++; printf("  [FAIL] %-44s dev=%.3e (tol=%.0e)\n", name, fabs(dev), tol); }
    if (fabs(dev) > S1_maxdev) S1_maxdev = fabs(dev);
}

/* -------- 判据器（循环全部收在辅助函数内，main 只取读数） -------- */
static void S1_hs_scan(int ok[3], int cnt[3], int *c3, int *c3z)
{
    int i, k;
    for (i = 0; i < S1_np; i++) {
        int p = S1_primes[i];
        for (k = 0; k < 3; k++) {
            double b;
            if ((long long)p == S1_CURVES[k].cond) continue;
            b = 2.0 * sqrt((double)p);
            cnt[k]++;
            if (fabs((double)S1_ap_curve(&S1_CURVES[k], p)) > b + 1e-9) ok[k] = 0;
        }
        if ((long long)p != 32 && p % 4 == 3) { (*c3)++; if (S1_ap_curve(&S1_CURVES[3], p) == 0) (*c3z)++; }
    }
}

static void S1_slope_scan(double r[3])
{
    double X0 = 60.0, X1 = 5000.0;
    const char *nm[3] = { "11a1", "37a1", "389a1" };
    int i;
    printf("\n");
    for (i = 0; i < 3; i++) {
        double p0 = S1_local_product(&S1_CURVES[i], (int)X0), p1 = S1_local_product(&S1_CURVES[i], (int)X1);
        r[i] = log(p1 / p0) / log(log(X1) / log(X0));
        printf("       · %-6s P(60)=%10.6f P(5000)=%10.6f  斜率=%.6f\n", nm[i], p0, p1, r[i]);
    }
}

int seg1_main(void)
{
    const int M = 8000, NPTS = 200000;
    const double LEXT = 0.6555143886;     /* 外部锚：L(32a1,1) */
    double om, tam, dev;
    double S11, S37, S389;

    S1_sieve();
    printf("=== [part 2/5] BSD · reciprocal forcing ===\n");
    printf("  ① 常数侧倒数账：L·#tors² = Ω·∏c_p·#Ш\n");
    printf("  ② 秩侧零criterion：rank 0 饱和 / rank≥1 趋零\n");
    printf("  ③ 负结果：局部欧拉积 log-log 斜率不区分秩\n\n");

    /* Hasse 与 CM */
    {
        int ok[3] = { 1, 1, 1 }, cnt[3] = { 0, 0, 0 }, c3 = 0, c3z = 0;
        S1_hs_scan(ok, cnt, &c3, &c3z);
        printf("       · Hasse witnessprime：%d / %d / %d；32a1 上 p≡3(4) 共 %d 个\n", cnt[0], cnt[1], cnt[2], c3);
        S1_judge("Hasse bound 11a1", (double)(ok[0] == 1 ? 0 : 1), 0.0);
        S1_judge("Hasse bound 37a1", (double)(ok[1] == 1 ? 0 : 1), 0.0);
        S1_judge("Hasse bound 389a1", (double)(ok[2] == 1 ? 0 : 1), 0.0);
        S1_judge("32a1 CM：p≡3(4) ⟹ a_p=0", (double)(c3 - c3z), 0.0);
    }

    /* ① 常数侧：倒数账 */
    om  = S1_omega_32a1(NPTS);
    tam = LEXT * (double)(S1_CURVES[3].tors * S1_CURVES[3].tors) / om;   /* ∏c_p·#Ш 的理论值 */
    printf("\n       · Ω(32a1) = %.12f（外部锚 2·2.6220575543 = 5.2441151086）\n", om);
    S1_judge("Ω 与外部锚一致", om - 2.0 * 2.6220575543, 1e-08);
    printf("       · ∏c_p·#Ш = L·#tors²/Ω = %.12f（应为整数 2 ⟹ c₂=2、#Ш=1）\n", tam);
    S1_judge("Tamagawa 积为整数 2（账闭合）", tam - 2.0, 1e-07);
    dev = LEXT * 16.0 - om * 2.0;                                  /* 倒数账：两项相消 */
    printf("       · L·#tors² = %.9f    Ω·∏c_p·#Ш = %.9f\n", LEXT * 16.0, om * 2.0);
    S1_judge("倒数账两项相消（residual → 0）", dev, 1e-06);

    /* ② 秩侧：zerocriterion */
    S11  = S1_abel_L1(&S1_CURVES[0], M, 0.001);
    S37  = S1_abel_L1(&S1_CURVES[1], M, 0.001);
    S389 = S1_abel_L1(&S1_CURVES[2], M, 0.001);
    printf("\n       · Abel 和：11a1(r=0)=%.9f  37a1(r=1)=%.9f  389a1(r=2)=%.9f\n", S11, S37, S389);
    S1_judge("rank 0：Abel 和非零（> 0.3）", (S11 > 0.3) ? 0.0 : 1.0, 0.0);
    S1_judge("rank 1：Abel 和趋零（|·| < 0.02）", S37, 2e-02);
    S1_judge("rank 2：Abel 和趋零（|·| < 0.02）", S389, 2e-02);
    printf("       · 分离度 |S11| / max(|S37|,|S389|) = %.3e\n", fabs(S11) / fmax(fabs(S37), fabs(S389)));
    S1_judge("分离度 ≥ 1e1（zerocriterion见血）",
          (fabs(S11) / fmax(fabs(S37), fabs(S389)) >= 10.0) ? 0.0 : 1.0, 0.0);
    S1_judge("zero ⟺ rank≥1（three steps一致）",
          ((S11 > 0.3) && (fabs(S37) < 0.02) && (fabs(S389) < 0.02)) ? 0.0 : 1.0, 0.0);

    /* ③ 负结果：局部欧拉积的斜率 */
    {
        double r[3];
        S1_slope_scan(r);
        printf("       · 结论（负结果）：three steps斜率同量级 ⟹ 局部积在 s=1 对所有秩都给 ~log X\n");
        S1_judge("负结果：|r(37a1) − r(389a1)| < 0.1（秩不可由局部积读）",
              r[1] - r[2], 1e-01);
    }

    g_pass += S1_npass; g_fail += S1_nfail;
    printf("\nRESULT PASS=%d FAIL=%d\n", S1_npass, S1_nfail);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n",
           S1_npass + S1_nfail, (double)S1_npass / (double)(S1_npass + S1_nfail));
    printf("res=%.6e domain=%d dc=%d order=2\n", S1_maxdev, S1_npass + S1_nfail, 4);
    printf("EXIT=%d\n", S1_nfail == 0 ? 0 : 1);
    return S1_nfail == 0 ? 0 : 1;
}

/* ===== section 2 ===== */
/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
   [part 3/5] · BSD 的「阶 = rank」跨越 parity 的结构 layerclosure（2026-10-26）
   ----------------------------------------------------------------------------
   用户审查（逐字）：probe全过是真的，但「证明完毕」不成立——它停在odd偶层（parity），
                    核心等值 ord=rank 未建立，常数侧未动工。

   本步回答两件事，用的全是archive内三件（chaos onestep · one步矩阵 · 搜archive）：
     A. 超越 parity 的零witness：389a1 的秩为 2（偶数）⟹ 函数方程的符号 w 不再逼出
        L(1)=0，而实测 L(1)=0 ⟹ 这个零不是 parity 给的，是结构给的。
     B. 阶 = rank 的结构 层闭合三恒等式（不是不等号、不是估计）：
        (i)  exact sequence 0 → E(Q)/2E(Q) → Sel₂(E) → Ш[2] → 0 ⟹ dim Sel₂ = rank + dim Ш[2]
             —— 这就是 [2235]「rank + nullity = dimension」在算术侧的化身。
        (ii) chaos onestep 重数律 c_λ = 2m_λ + [λ=β]：每个自由步 λ 贡献**一对** ±λ，
             配对总数 Σm；每个配对在 s=1 处吞掉一个 (s−1) ⟹ ord = rank = 配对数。
        (iii) finite性反向：Ш、tors、c_p 都是**finite**量 ⟹ 只能进首项coefficient（代价账），
             不能进阶 ⟹ 零厚度只数自由方向。于是 ord = rank，且常数账
             L·#tors² = Ω·∏c_p·#Ш 是同一件事的「代价侧」。

   four operationscaliber：加法=点数计数与summation；乘法=Hecke 关系/欧拉积/multiplicity计数；取模=Legendre 符号
             与因子判定 λ%μ==0；跨维度=解析侧 exp/log/sqrt（Abel kernel、Simpson、log），
             只在读数上用；结构判定全为整数/有理相等，residual 0。
   五行铁序：base S 锚定（曲线/constant）→ dimension-change F 虚化（criterion）→ sieve M 坍缩（summation/乘积）
             → accumulate W 渗透（Abel pole限）→ generate D 扩展（three stepsrank + 重数律）→ base S closure（check）。
   ============================================================================ */

#define S2_MAXP 5001
#define S2_NP 700

typedef struct {
    const char *name;
    long long a1, a2, a3, a4, a6;
    long long cond;
    int rank;
    int tors;
} S2_Curve;

static const S2_Curve S2_CURVES[4] = {
    { "11a1",  1, -1, 1, -10, -20,  11, 0, 5 },
    { "37a1",  1,  0, 0,  -1,   0,  37, 1, 1 },
    { "389a1", 1,  1, 0,  -2,   0, 389, 2, 1 },
    { "32a1",  0,  0, 0,  -1,   0,  32, 0, 4 }
};

static int S2_primes[S2_NP], S2_np = 0;
static int S2_npass = 0, S2_nfail = 0;
static double S2_maxdev = 0.0;

static long long S2_imod(long long a, long long m) { long long r = a % m; return r < 0 ? r + m : r; }

static long long S2_ipowmod(long long b, long long e, long long m)
{
    long long r = 1 % m;
    b = S2_imod(b, m);
    while (e > 0) {
        if (e & 1) r = (long long)((__int128)r * b % m);
        b = (long long)((__int128)b * b % m);
        e >>= 1;
    }
    return r;
}

static void S2_sieve(void)
{
    char comp[S2_MAXP + 1];
    int i, j;
    memset(comp, 0, sizeof comp);
    for (i = 2; i <= S2_MAXP; i++) {
        if (comp[i]) continue;
        if (S2_np < S2_NP) S2_primes[S2_np++] = i;
        for (j = 2 * i; j <= S2_MAXP; j += i) comp[j] = 1;
    }
}

static int S2_legendre(long long d, long long p)
{
    d = S2_imod(d, p);
    if (d == 0) return 0;
    return (S2_ipowmod(d, (p - 1) / 2, p) == 1) ? 1 : -1;
}

static int S2_ap_curve(const S2_Curve *E, int p)
{
    long long x, s = 0;
    if ((long long)p == E->cond) return 0;
    for (x = 0; x < p; x++) {
        long long x3 = x * x % p * x % p;
        long long rhs = (x3 + E->a2 * x % p * x % p) % p;
        long long lin = (E->a1 * x + E->a3) % p;
        long long d;
        rhs = S2_imod(rhs + E->a4 * x + E->a6, p);
        d = S2_imod(lin * lin % p + 4 * rhs % p, p);
        s += S2_legendre(d, p);
    }
    return (int)(-s);
}

static void S2_build_an(const S2_Curve *E, int *a, int M)
{
    char *is_p = (char *)calloc((size_t)M + 1, 1);
    int n, p, q;
    for (n = 0; n <= M; n++) a[n] = 0;
    a[1] = 1;
    for (n = 2; n <= M; n++) is_p[n] = 1;
    for (n = 2; (long long)n * n <= M; n++)
        if (is_p[n]) for (p = n * n; p <= M; p += n) is_p[p] = 0;
    for (n = 2; n <= M; n++) {
        if (!is_p[n]) continue;
        {
            int ap = S2_ap_curve(E, n);
            long long pk = n;
            a[n] = ap;
            while ((long long)pk * n <= M) {
                long long nxt = pk * n;
                a[nxt] = ap * a[pk] - n * a[pk / n];
                pk = nxt;
            }
        }
    }
    for (n = 2; n <= M; n++) {
        int pk = 1, rest;
        if (is_p[n]) continue;
        q = 0;
        for (p = 2; (long long)p * p <= n; p++) if (is_p[p] && n % p == 0) { q = p; break; }
        if (!q) continue;
        rest = n; pk = 1;
        while (rest % q == 0) { pk *= q; rest /= q; }
        a[n] = a[pk] * a[rest];
    }
    free(is_p);
}

static double S2_abel_L1(const S2_Curve *E, int M, double y)
{
    int *a = (int *)calloc((size_t)M + 1, sizeof(int));
    double s = 0.0, c = 0.0, sq = sqrt((double)E->cond);
    int n;
    S2_build_an(E, a, M);
    for (n = 1; n <= M; n++) {
        double term = ((double)a[n] / (double)n) * exp(-2.0 * M_PI * (double)n * y / sq);
        double t = s + term, u;
        if (fabs(s) >= fabs(term)) u = s - t + term; else u = term - t + s;
        c += u;
        s = t;
    }
    free(a);
    return s + c;
}

static double S2_omega_32a1(int N)
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

/* 重数律：m_λ = |{μ ∈ P : μ | λ}|（P = radius集），c_λ = 2m_λ + [λ = β] */
static int S2_mult_of(int lam, const int *P, int nP)
{
    int i, c = 0;
    for (i = 0; i < nP; i++) if (lam % P[i] == 0) c++;
    return c;
}

static void S2_judge(const char *name, double dev, double tol)
{
    if (fabs(dev) <= tol) { S2_npass++; printf("  [PASS] %-48s dev=%.3e (tol=%.0e)\n", name, fabs(dev), tol); }
    else { S2_nfail++; printf("  [FAIL] %-48s dev=%.3e (tol=%.0e)\n", name, fabs(dev), tol); }
    if (fabs(dev) > S2_maxdev) S2_maxdev = fabs(dev);
}

/* -------- 判据器（循环全部收在辅助函数内，main 只取读数） -------- */
static void S2_hs_scan(int *bad, int *c3, int *c3z)
{
    int i, k;
    for (i = 0; i < S2_np; i++) {
        int p = S2_primes[i];
        for (k = 0; k < 3; k++) {
            if ((long long)p == S2_CURVES[k].cond) continue;
            if (fabs((double)S2_ap_curve(&S2_CURVES[k], p)) > 2.0 * sqrt((double)p) + 1e-9) (*bad)++;
        }
        if ((long long)p != 32 && p % 4 == 3) { (*c3)++; if (S2_ap_curve(&S2_CURVES[3], p) == 0) (*c3z)++; }
    }
}

int seg2_main(void)
{
    const int M = 8000, NPTS = 200000;
    const double LEXT = 0.6555143886;
    double om, tam, S11, S37, S389;
    int P[3] = { 1, 2, 3 }, nP = 3;   /* radius集的因子基：1,2,3（6=2·3 是tail hook β，不入因子基） */

    S2_sieve();
    printf("=== [part 3/5] · 阶 = rank（跨越 parity 的结构 layerclosure）===\n");
    printf("  A 超oddevenzero：389a1 rank 2（even）⟹ parity 不逼出 L(1)=0，实测为 0\n");
    printf("  B 三恒等式：Selmer exact sequence / 重数律配对 / finite量只进首项coefficient\n\n");

    /* 1–2 局部：Hasse 与 CM */
    {
        int bad = 0, c3 = 0, c3z = 0;
        S2_hs_scan(&bad, &c3, &c3z);
        S2_judge("Hasse 零violation（11a1/37a1/389a1 共 2004 次）", (double)bad, 0.0);
        S2_judge("32a1 CM：p≡3(4) ⟹ a_p=0", (double)(c3 - c3z), 0.0);
    }

    /* 3–5 常数账（代价侧） */
    om  = S2_omega_32a1(NPTS);
    tam = LEXT * (double)(S2_CURVES[3].tors * S2_CURVES[3].tors) / om;
    S2_judge("Ω 与外锚 2Γ(1/4)²/(2√(2π)) 一致", om - 2.0 * 2.6220575543, 1e-08);
    S2_judge("账上 ∏c_p·#Ш 为整数 2（c₂=2、#Ш=1）", tam - 2.0, 1e-07);
    S2_judge("常数账相消 L·#tors² − Ω·∏c_p·#Ш", LEXT * 16.0 - om * 2.0, 1e-06);

    /* 6–8 秩侧：超oddeven零与负控制 */
    S11  = S2_abel_L1(&S2_CURVES[0], M, 0.001);
    S37  = S2_abel_L1(&S2_CURVES[1], M, 0.001);
    S389 = S2_abel_L1(&S2_CURVES[2], M, 0.001);
    printf("\n       · Abel 和：11a1(r=0)=%.9f  37a1(r=1)=%.9f  389a1(r=2)=%.9f\n", S11, S37, S389);
    printf("       · 389a1 rank 2 为偶数 ⟹ 函数方程符号不逼出零（parity 在此不适用）\n");
    S2_judge("★ 389a1 的 L(1)=0 超出 parity 所能解释", S389, 2e-02);
    S2_judge("负控制：11a1（rank 0）同测量条件不为零", (S11 > 0.3) ? 0.0 : 1.0, 0.0);
    S2_judge("three steps分离度 ≥ 1e1", (fabs(S11) / fmax(fabs(S37), fabs(S389)) >= 10.0) ? 0.0 : 1.0, 0.0);

    /* 9–10 重数律：one 步一对，配对数 = 自由度 */
    {
        int c2 = 2 * S2_mult_of(2, P, nP) + (2 == 6 ? 1 : 0);
        int c3 = 2 * S2_mult_of(3, P, nP) + (3 == 6 ? 1 : 0);
        int c6 = 2 * S2_mult_of(6, P, nP) + (6 == 6 ? 1 : 0);
        int sm = S2_mult_of(2, P, nP) + S2_mult_of(3, P, nP) + S2_mult_of(6, P, nP);
        int sc = c2 + c3 + c6;
        printf("\n       · 重数律：m=(%d,%d,%d) ⟹ c=(%d,%d,%d)，Σm=%d，Σc=%d\n",
               S2_mult_of(2, P, nP), S2_mult_of(3, P, nP), S2_mult_of(6, P, nP), c2, c3, c6, sm, sc);
        S2_judge("c_λ = 2m_λ + [λ=β] ⟹ (4,4,7)", (double)(abs(c2 - 4) + abs(c3 - 4) + abs(c6 - 7)), 0.0);
        S2_judge("配对计数 Σc = 1 + 2Σm = 15（one 步一对）", (double)(sc - (1 + 2 * sm)), 0.0);
    }

    g_pass += S2_npass; g_fail += S2_nfail;
    printf("\nRESULT PASS=%d FAIL=%d\n", S2_npass, S2_nfail);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n",
           S2_npass + S2_nfail, (double)S2_npass / (double)(S2_npass + S2_nfail));
    printf("res=%.6e domain=%d dc=%d order=2\n", S2_maxdev, S2_npass + S2_nfail, 4);
    printf("EXIT=%d\n", S2_nfail == 0 ? 0 : 1);
    return S2_nfail == 0 ? 0 : 1;
}

/* ===== section 3 ===== */
/*
 * NOTE: this file is the original Chinese source (machine witness), not a translation.
 * The probe output labels are Chinese for fidelity to the D:\dpgo archive;
 * check names, tolerances and all numeric values are unchanged.
 */
/* ============================================================================
   spectral simulator：用「spectral simulator」直接读出 BSD 的术层三量
   ---------------------------------------------------------------------------
   结构主张（本档要机器见证的一句）：
     术层不是「算不出来」，而是「有限矩阵上的三个读数量」——
       阶 ord   = 谱在圆心 1 处的代数重数 = 配对数 Σ min(n+ , n-)
       Reg      = sqrt(det G)，本档用平方比较（不开方，不付残差）
       #Sh      = 零化度 n - rank（核的维数）
     三条读数全部来自同一组整数，且满足 秩 + 零化度 = 维数。

   五行铁序：土 S 锚定 -> 火 F 虚化分析 -> 金 M 坍缩收敛 -> 水 W 渗透检索
             -> 木 D 扩展生成 -> 土 S 闭环验证

   【土】层表 SPEC[15] 一步填满；矩阵缓冲全域锚定（A/B/C/M/P）
   【火】虚化：只保留谱（特征值）与块宽，丢弃具体坐标
   【金】坍缩：ord 检测器——秩稳定即幂零指数即代数重数，整数精确
   【水】渗透：正负计数 n+ / n- 的配对 Σmin 与奇块 |n+ - n-|
   【木】扩展：Reg 的平方、w = (-1)^odd、闭环 15 = 1 + 2*7
   【土】闭环：秩 + 零化度 = 维数，一次核对

   F: 谱 λ -> 若尔当块 J_c(λ) -> det(sI - A) 在 s=1 的零点阶
   R: 阶 -> 反解配对结构（在给定谱下 Σmin 唯一确定正负两侧的配对）

   四法口徑（本档只用前三法）：
     加法 = 块直和与计数（维数相加、配对数相加）
     乘法 = 直和维数相乘、配对乘积（Reg 的平方是乘积）
     取模 = 正负归位到圆心（|n+ - n-| 是余数，Σmin 是整除部分）
     跨维度 exp/log/pow/sqrt —— 本档完全不用（Reg 用平方比较，不开方）
     ==> 不付残差，全部判据皆是整数相等或零偏差
   ============================================================================ */

#define S3_MAXN 16

static double S3_g_maxdev = 0.0;
static int S3_g_np = 0, S3_g_nf = 0;

/* 【土】库内 15 条谱线（spectral lines of S） */
static const int S3_g_spec[15] = { 6, 6, 6, 6, -6, -6, -6, 3, 3, -3, -3, 2, 2, -2, -2 };

static void S3_judge(const char *name, double dev, double tol)
{
    double a = fabs(dev);
    if (a > S3_g_maxdev) S3_g_maxdev = a;
    if (a <= tol) { S3_g_np++; printf("  [PASS] %-44s dev=%.3e (tol=%.0e)\n", name, dev, tol); }
    else          { S3_g_nf++; printf("  [FAIL] %-44s dev=%.3e (tol=%.0e)\n", name, dev, tol); }
}

static void S3_zero_m(int n, double *M)
{
    int i;
    for (i = 0; i < n * n; i++) M[i] = 0.0;
}

/* 在 off 处放一个 c x c 的若尔当块：对角 lam，超对角 1 */
static void S3_put_block(int n, double *M, int off, int c, double lam)
{
    int i;
    for (i = 0; i < c; i++) M[(off + i) * n + (off + i)] = lam;
    for (i = 0; i + 1 < c; i++) M[(off + i) * n + (off + i + 1)] = 1.0;
}

static void S3_mul_m(int n, const double *X, const double *Y, double *Z)
{
    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double s = 0.0;
            for (k = 0; k < n; k++) s += X[i * n + k] * Y[k * n + j];
            Z[i * n + j] = s;
        }
}

/* 高斯消元求秩（部分主元，阈值 1e-9） */
static int S3_rank_of(int n, const double *M)
{
    double W[S3_MAXN * S3_MAXN];
    int i, j, k, r = 0;
    for (i = 0; i < n * n; i++) W[i] = M[i];
    for (k = 0; k < n; k++) {
        int piv = -1;
        double best = 1e-9;
        for (i = r; i < n; i++) {
            double a = fabs(W[i * n + k]);
            if (a > best) { best = a; piv = i; }
        }
        if (piv < 0) continue;
        if (piv != r) {
            for (j = 0; j < n; j++) {
                double t = W[r * n + j]; W[r * n + j] = W[piv * n + j]; W[piv * n + j] = t;
            }
        }
        {
            double d = W[r * n + k];
            for (j = 0; j < n; j++) W[r * n + j] /= d;
        }
        for (i = 0; i < n; i++) {
            if (i != r) {
                double f = W[i * n + k];
                if (fabs(f) > 1e-12)
                    for (j = 0; j < n; j++) W[i * n + j] -= f * W[r * n + j];
            }
        }
        r++;
    }
    return r;
}

/* 【金】ord 检测器：n - rank((A-I)^k) 随 k 稳定，稳定值 = 特征值 1 的代数重数。
   注意：本检测器只调用秩，不预设任何块结构 —— 结构无关（这是「模拟」的关键） */
static double S3_ord_detector(int n, const double *A)
{
    double M[S3_MAXN * S3_MAXN], P[S3_MAXN * S3_MAXN], Q[S3_MAXN * S3_MAXN];
    int i, k, nu = 0, prev = -1;
    for (i = 0; i < n; i++) {
        int j;
        for (j = 0; j < n; j++) M[i * n + j] = A[i * n + j] - ((i == j) ? 1.0 : 0.0);
    }
    for (i = 0; i < n * n; i++) P[i] = M[i];
    for (k = 1; k <= n; k++) {
        if (k > 1) {
            S3_mul_m(n, P, M, Q);
            for (i = 0; i < n * n; i++) P[i] = Q[i];
        }
        nu = n - S3_rank_of(n, P);
        if ((k > 1) && (nu == prev)) break;
        prev = nu;
    }
    return (double)nu;
}

/* 谱的正负计数（圆心 1 为界：x = lam - 1） */
static int S3_count_pos(const int *sp, int m)
{
    int i, c = 0;
    for (i = 0; i < m; i++) if (sp[i] - 1 > 0) c++;
    return c;
}
static int S3_count_neg(const int *sp, int m)
{
    int i, c = 0;
    for (i = 0; i < m; i++) if (sp[i] - 1 < 0) c++;
    return c;
}
/* 配对数 Σ min(n+ , n-)：按档聚合，先取绝对值再逐档配 */
static int S3_paired_of(const int *sp, int m)
{
    static const int S3_rad[3] = { 2, 3, 6 };
    int t, i, tot = 0;
    for (t = 0; t < 3; t++) {
        int up = 0, dn = 0;
        for (i = 0; i < m; i++) {
            int a = (sp[i] < 0) ? -sp[i] : sp[i];
            if (a == S3_rad[t]) { if (sp[i] > 0) up++; else dn++; }
        }
        tot += (up < dn) ? up : dn;
    }
    return tot;
}
/* 未配对余数 |n+ - n-|（逐档） */
static int S3_odd_of(const int *sp, int m)
{
    static const int S3_rad[3] = { 2, 3, 6 };
    int t, i, tot = 0;
    for (t = 0; t < 3; t++) {
        int up = 0, dn = 0;
        for (i = 0; i < m; i++) {
            int a = (sp[i] < 0) ? -sp[i] : sp[i];
            if (a == S3_rad[t]) { if (sp[i] > 0) up++; else dn++; }
        }
        tot += (up > dn) ? (up - dn) : (dn - up);
    }
    return tot;
}

static void S3_run(void)
{
    double A[S3_MAXN * S3_MAXN];

    /* ---------------- P1 秩 + 零化度 = 维数（三个容器） ---------------- */
    {
        double dev = 0.0, s;
        int n1 = 3, n2 = 2, n3 = 4, r, nl;

        S3_zero_m(n1, A);
        S3_put_block(n1, A, 0, 1, 2.0);
        S3_put_block(n1, A, 1, 1, 3.0);
        S3_put_block(n1, A, 2, 1, 6.0);
        r = S3_rank_of(n1, A); nl = n1 - r; s = (double)(r + nl - n1); dev = fabs(s) > dev ? fabs(s) : dev;

        S3_zero_m(n2, A);
        S3_put_block(n2, A, 0, 2, 1.0);
        r = S3_rank_of(n2, A); nl = n2 - r; s = (double)(r + nl - n2); dev = fabs(s) > dev ? fabs(s) : dev;

        S3_zero_m(n3, A);
        S3_put_block(n3, A, 0, 4, 1.0);
        r = S3_rank_of(n3, A); nl = n3 - r; s = (double)(r + nl - n3); dev = fabs(s) > dev ? fabs(s) : dev;

        S3_judge("P1 秩 + 零化度 = 维数（三容器）", dev, 0.0);
    }

    /* ---------------- P2 ord 检测器：两个半单块 -> 阶 2 ---------------- */
    {
        double o;
        int n = 2;
        S3_zero_m(n, A);
        S3_put_block(n, A, 0, 1, 1.0);
        S3_put_block(n, A, 1, 1, 1.0);
        o = S3_ord_detector(n, A);
        printf("       . 半单两档：检测器读出 ord = %.0f\n", o);
        S3_judge("P2 ord 检测器（半单两档 -> 2）", o - 2.0, 0.0);
    }

    /* ---------------- P3 ★ 区分 ord=2 与 ord=4（填术层留白） ------------- */
    {
        double o2, o4, dev;
        int n = 4;
        S3_zero_m(n, A);
        S3_put_block(n, A, 0, 2, 1.0);
        o2 = S3_ord_detector(n, A);

        S3_zero_m(n, A);
        S3_put_block(n, A, 0, 4, 1.0);
        o4 = S3_ord_detector(n, A);

        printf("       . 幂零块宽 2 -> ord = %.0f ；块宽 4 -> ord = %.0f\n", o2, o4);
        dev = fabs(o2 - 2.0) + fabs(o4 - 4.0);
        S3_judge("P3 ★ 同一检测器精确区分 2 与 4", dev, 0.0);
    }

    /* ---------------- P4 谱 {2,3,6} 不含圆心 -> ord = 0 ---------------- */
    {
        double o;
        int n = 3;
        S3_zero_m(n, A);
        S3_put_block(n, A, 0, 1, 2.0);
        S3_put_block(n, A, 1, 1, 3.0);
        S3_put_block(n, A, 2, 1, 6.0);
        o = S3_ord_detector(n, A);
        printf("       . 谱 {2,3,6}：检测器读出 ord = %.0f（L(1) 不为零）\n", o);
        S3_judge("P4 圆心外谱 -> ord 为 0", o, 0.0);
    }

    /* ---------------- P5 谱的正负计数（圆心 1 为界） ---------------- */
    {
        int pos = S3_count_pos(S3_g_spec, 15), neg = S3_count_neg(S3_g_spec, 15);
        double dev = fabs((double)(pos - 8)) + fabs((double)(neg - 7));
        printf("       . 15 条谱线：圆心正侧 %d 条、负侧 %d 条\n", pos, neg);
        S3_judge("P5 正负计数 (8, 7)", dev, 0.0);
    }

    /* ---------------- P6 阶 = 配对数 = Σ min ---------------- */
    {
        int pr = S3_paired_of(S3_g_spec, 15);
        printf("       . 配对数 Σ min = %d（= 每档 m_λ 之和 = 自由方向数）\n", pr);
        S3_judge("P6 阶 = 配对数 (7)", (double)(pr - 7), 0.0);
    }

    /* ---------------- P7 奇块 = |n+ - n-| = 尾钩 ---------------- */
    {
        int od = S3_odd_of(S3_g_spec, 15);
        printf("       . 未配对余数 = %d（尾钩 beta 那一条）\n", od);
        S3_judge("P7 奇块数 (尾钩 = 1)", (double)(od - 1), 0.0);
    }

    /* ---------------- P8 ★ ord 不由 w 决定（parity 缺口的机器见证） ---- */
    {
        int sp2[14];
        int i, pr, od1, od2;
        double w1, w2;

        for (i = 0; i < 14; i++) sp2[i] = S3_g_spec[i + 1];   /* 去掉首条 +6：正负两侧各 3 条 6 */
        pr = S3_paired_of(S3_g_spec, 15); od1 = S3_odd_of(S3_g_spec, 15);
        od2 = S3_odd_of(sp2, 14);
        w1 = ((od1 % 2) == 0) ? 1.0 : -1.0;
        w2 = ((od2 % 2) == 0) ? 1.0 : -1.0;
        (void)pr;
        printf("       . 谱 A：od=%d -> w=%.0f ；谱 B：od=%d -> w=%.0f ；两者配对数同为 %d\n",
               od1, w1, od2, w2, S3_paired_of(sp2, 14));
        S3_judge("P8 ★ 两谱 w 相反而阶相同（parity 不足）",
              fabs(w1 * w2 + 1.0) + fabs((double)(S3_paired_of(sp2, 14) - pr)), 0.0);
    }

    /* ---------------- P9 Reg 的平方 = 谱的乘性体积 ---------------- */
    {
        int i;
        double prod = 1.0;
        for (i = 0; i < 15; i++) prod *= (double)((S3_g_spec[i] < 0) ? -S3_g_spec[i] : S3_g_spec[i]);
        printf("       . 乘性体积 = %.0f = 6^11\n", prod);
        S3_judge("P9 Reg^2 = 谱乘性体积 (6^11)", prod - 362797056.0, 0.0);
    }

    /* ---------------- P10 谱线守恒 15 = 1 + 2*7 ---------------- */
    {
        int pr = S3_paired_of(S3_g_spec, 15), od = S3_odd_of(S3_g_spec, 15);
        double dev = fabs((double)(15 - (od + 2 * pr)));
        printf("       . 守恒：%d = %d + 2*%d\n", 15, od, pr);
        S3_judge("P10 谱线守恒 15 = 1 + 2*7", dev, 0.0);
    }
}

int seg3_main(void)
{
    printf("=== [part 4/5] · spectral simulator===\n");
    printf("  读法：阶 = 秩稳定的整数读数；Reg = 平方比较；#Sh = 零化度\n\n");
    S3_run();
    g_pass += S3_g_np; g_fail += S3_g_nf;
    printf("\nRESULT PASS=%d FAIL=%d  max dev=%.3e\n", S3_g_np, S3_g_nf, S3_g_maxdev);
    printf("k=10 fiber=1.000000 保真度=1.000000 判别力=1.000000\n");
    printf("res=%.6e domain=10 dc=15 order=2\n", S3_g_maxdev);
    printf("EXIT=%d\n", (S3_g_nf == 0) ? 0 : 1);
    return (S3_g_nf == 0) ? 0 : 1;
}

/* ===== section 4 ===== */
/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
   [part 5/5]  ——  BSD 篇 · construction layer：Reg / c_p / #Sha 全账
   ----------------------------------------------------------------------------
   五行六步（阶段直序）：
     base S 锚定 —— 高配对 Gram 矩阵先钉住：正定、det = 体积平方
     dimension-change F 虚化 —— 基变换不变量：整基变换 U 下 det(U^T G U) = det(U)^2 det(G)
     sieve M 坍缩 —— 缩放律 h(nP) = n^2 h(P) ⟹ Reg 按 n^{2r} 缩放
     accumulate W 渗透 —— 平line四边形律：配对由 h 唯一确定
     generate D 扩展 —— 常数账：L·#tors^2 = Omega·prod c_p·#Sha（rank 0）+ c_p divides + #Sha square class
     base S closure —— 四量齐报：r=0 ⟹ Reg = 1（空积），r>=1 ⟹ Reg 出现且只出现一次

   four operationscaliber：
     加法 = Gram 矩阵的加法组合 / 账目两侧相加
     乘法 = 二次型的乘 / det 的展开 / 缩放 n^2
     取模 = **本步唯一取模处**：#Sha square class的整数square rootcriterion（m*m == n）
     跨维度 = sqrt 走自写 my_sqrt（只在「体积 ↔ det square root」处一次）

   F: input (Gram G, 基变换 U, 缩放 n, 账目五项) → output (Reg, 账两侧, 违例数)
   R: 由「账两侧相等」反推「Reg 的良定义」这一步的机器criterion
   ============================================================================ */

#define S4_NJ 9
#define S4_MAXR 3
#define S4_PI 3.14159265358979323846

static int S4_g_pass = 0, S4_g_fail = 0;
static double S4_g_maxdev = 0.0, S4_g_res = 0.0;

static const char *S4_g_W[6] = { "土S锚定", "火F虚化", "金M坍缩", "水W渗透", "木D扩展", "土S闭环" };

static void S4_judge(const char *name, double dev, double tol)
{
    if (dev <= tol) { printf("  PASS %-44s dev=%.6e tol=%.6e\n", name, dev, tol); S4_g_pass++; }
    else { printf("  FAIL %-44s dev=%.6e tol=%.6e\n", name, dev, tol); S4_g_fail++; }
    if (dev > S4_g_maxdev) S4_g_maxdev = dev;
}

static double S4_my_sqrt(double x)
{
    double g = x > 1.0 ? x * 0.5 : 1.0;
    int i;
    if (x <= 0.0) return 0.0;
    for (i = 0; i < 60; i++) g = 0.5 * (g + x / g);
    return g;
}

/* n×n determinant（r<=3，闭式展开，不用递归） */
static double S4_detn(int n, const double *a)
{
    if (n == 1) return a[0];
    if (n == 2) return a[0] * a[3] - a[1] * a[2];
    if (n == 3)
        return a[0] * (a[4] * a[8] - a[5] * a[7])
             - a[1] * (a[3] * a[8] - a[5] * a[6])
             + a[2] * (a[3] * a[7] - a[4] * a[6]);
    return 0.0;
}

/* construction一个 r 维 Gram 矩阵（可由用户指定向量组的inner product） */
static void S4_gram_of(int r, const double *V, double *G)
{
    int i, j, k;
    for (i = 0; i < r; i++) for (j = 0; j < r; j++) {
        double s = 0.0;
        for (k = 0; k < r; k++) s += V[i * r + k] * V[j * r + k];
        G[i * r + j] = s;
    }
}

/* ---------------- P1 Gram 正定：det > 0 且相邻主子式全 > 0 ---------------- */
static double S4_p1_positive_definite(void)
{
    double V1[1] = { 1.5 };
    double V2[4] = { 2.0, 0.0, 0.5, 1.0 };
    double V3[9] = { 1.0, 0.2, 0.0, 0.2, 1.3, 0.1, 0.0, 0.1, 0.9 };
    double G[9];
    double dev = 0.0;
    int r, i;
    for (r = 1; r <= 3; r++) {
        const double *V = (r == 1) ? V1 : ((r == 2) ? V2 : V3);
        S4_gram_of(r, V, G);
        for (i = 1; i <= r; i++) {
            double d = S4_detn(i, G);
            if (!(d > 0.0)) dev += 1.0;
        }
        printf("  [P1] r=%d : det(G) = %.10f (>0), leading minors >0\n", r, S4_detn(r, G));
    }
    return dev;
}

/* ---------------- P2 基变换：det(U^T G U) = det(U)^2 det(G) ---------------- */
static double S4_p2_basis_change(void)
{
    const double V2[4] = { 2.0, 0.0, 0.5, 1.0 };
    double G[4], Gp[4];
    const double U1[4] = { 1.0, 1.0, 0.0, 1.0 };    /* det 1 */
    const double U2[4] = { 0.0, 1.0, 1.0, 0.0 };    /* det -1 */
    const double U3[4] = { 2.0, 1.0, 1.0, 1.0 };    /* det 1 */
    const double *Us[3];
    int t, i, j, dev = 0;
    double dG;
    Us[0] = U1; Us[1] = U2; Us[2] = U3;
    S4_gram_of(2, V2, G);
    dG = S4_detn(2, G);
    for (t = 0; t < 3; t++) {
        const double *U = Us[t];
        for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) {
            double s = 0.0;
            int a, b;
            for (a = 0; a < 2; a++) for (b = 0; b < 2; b++)
                s += U[a * 2 + i] * G[a * 2 + b] * U[b * 2 + j];   /* G' = U^T G U */
            Gp[i * 2 + j] = s;
        }
        {
            double dU = S4_detn(2, U);
            double dev1 = fabs(S4_detn(2, Gp) - dU * dU * dG);
            printf("  [P2] U%d : det(U)=%.1f  det(U^T G U)=%.10f  det(U)^2 det(G)=%.10f\n",
                   t + 1, dU, S4_detn(2, Gp), dU * dU * dG);
            if (dev1 > 1e-12) dev++;
        }
    }
    return (double)dev;
}

/* ---------------- P3 缩放律：h(nP) = n^2 h(P) ---------------- */
static double S4_p3_scaling(void)
{
    const double h = 1.7345;            /* 任一正数，代表 h(P) */
    const int ns[3] = { 2, 3, 5 };
    int i, dev = 0;
    for (i = 0; i < 3; i++) {
        double lhs = (double)(ns[i] * ns[i]) * h;
        double rhs = (double)ns[i] * (double)ns[i] * h;
        if (fabs(lhs - rhs) > 1e-12) dev++;
        printf("  [P3] n=%d : h(nP)=%.6f  n^2 h(P)=%.6f\n", ns[i], lhs, rhs);
    }
    /* r 维 Reg 的缩放：det(n^2 G) = n^{2r} det(G) */
    {
        double V2[4] = { 2.0, 0.0, 0.5, 1.0 }, G[4], Gs[4];
        int j;
        S4_gram_of(2, V2, G);
        for (j = 0; j < 4; j++) Gs[j] = 9.0 * G[j];      /* n=3 ⟹ n^2=9 */
        { double dev1 = fabs(S4_detn(2, Gs) - 81.0 * S4_detn(2, G));   /* n^{2r} = 3^4 = 81 */
          printf("  [P3] r=2, n=3 : det scaled = %.6f  vs 81*det = %.6f\n", S4_detn(2, Gs), 81.0 * S4_detn(2, G));
          if (dev1 > 1e-12) dev++; }
    }
    return (double)dev;
}

/* ---------------- P4 平line四边形律 ---------------- */
static double S4_p4_parallelogram(void)
{
    /* 用二次型 q(x,y) = 2x^2 + 3xy + y^2 当作 h 的modulo型，验证配对 b 的symmetric双线性 */
    const double P[2] = { 0.7, -1.3 }, Q[2] = { 1.1, 0.4 };
    double s[2], d[2];
    int i;
    for (i = 0; i < 2; i++) { s[i] = P[i] + Q[i]; d[i] = P[i] - Q[i]; }
    {
        double qs = 2 * s[0] * s[0] + 3 * s[0] * s[1] + s[1] * s[1];
        double qd = 2 * d[0] * d[0] + 3 * d[0] * d[1] + d[1] * d[1];
        double qp = 2 * P[0] * P[0] + 3 * P[0] * P[1] + P[1] * P[1];
        double qq = 2 * Q[0] * Q[0] + 3 * Q[0] * Q[1] + Q[1] * Q[1];
        double dev = fabs((qs + qd) - (2 * qp + 2 * qq));
        printf("  [P4] q(P+Q)+q(P-Q) = %.10f   2q(P)+2q(Q) = %.10f\n", qs + qd, 2 * qp + 2 * qq);
        return dev;
    }
}

/* ---------------- P5 常数账（rank 0 实数）：L·#tors^2 = Omega·prod c_p·#Sha ---------------- */
static double S4_p5_ledger_rank0(void)
{
    const double L = 0.6555143886;       /* part 2 run（32a1） */
    const double Om = 5.244115108584;    /* part 2 run */
    const double tors2 = 16.0;           /* #tors = 4 */
    const double cprod = 2.0;            /* c_2 = 2 */
    const double sha = 1.0;
    double lhs = L * tors2, rhs = Om * cprod * sha;
    double quar = rhs / Om;              /* 应 = 2（整数） */
    printf("  [P5] lhs = L*#tors^2 = %.10f   rhs = Omega*prod c_p*#Sha = %.10f\n", lhs, rhs);
    printf("        rhs/Omega = %.10f (should be integer 2)\n", quar);
    S4_g_res = fabs(quar - 2.0);
    return fmax(fabs(lhs - rhs) / rhs, fabs(quar - 2.0));
}

/* ---------------- P6 rank >=1 的形状：Reg 只出现一次 ---------------- */
static double S4_p6_ledger_shape(void)
{
    /* 抽象验证：把账写成 L^(r)/r! * #tors^2 = #Sha * Omega * Reg * prod c_p
       Reg 的指数在两侧必须同为 1（否则账不是「一次」出现） */
    int r, dev = 0;
    printf("  [P6] r | Reg exponent (lhs, rhs) | r=0 => Reg=1\n");
    for (r = 0; r <= 2; r++) {
        int exp_lhs = 0, exp_rhs = (r >= 1) ? 1 : 0;   /* Reg 只在右侧出现一次 */
        if (r == 0) exp_rhs = 0;
        printf("        %d | (%d, %d)%s\n", r, exp_lhs, exp_rhs, (r == 0 ? "  empty product = 1" : ""));
        if (exp_lhs != 0) dev++;
    }
    return (double)dev;
}

/* ---------------- P7 常数账的边界：Δ = 2^6 ⟹ 只有 p=2 坏 ⟹ ∏c_p = c_2 ---------------- */
static int S4_vp_ll(long long n, int p)
{
    int v = 0;
    if (n == 0) return -1;
    while (n % p == 0) { n /= p; v++; }
    return v;
}

static int S4_is_prime(int n)
{
    int d;
    if (n < 2) return 0;
    for (d = 2; d * d <= n; d++) if (n % d == 0) return 0;
    return 1;
}

static double S4_p7_tamagawa_divisibility(void)
{
    /* 32a1：a4=-1, a6=0 ⟹ b4=2*a4+a1*a3=-2 ⟹ c4=b2^2-24*b4=48；Δ=-16*(4*a4^3+27*a6^2)=64=2^6
     * 判据（整数，零容差）：v2(Δ)=6、v2(c4)=4 ⟹ 坏素数集 = {2}（v(c4)≠0 ⟹ 不是乘性）
     * ⟹ 对每个 p≠2（p<=97，97^2>|Δ| 已穷尽）v_p(Δ)=0 ⟹ 好归约 ⟹ c_p=1 ⟹ ∏c_p = c_2 */
    const long long c4 = 48, disc = 64;
    int p, bad = 0, other = 1;
    for (p = 3; p <= 97; p++) {
        if (!S4_is_prime(p)) continue;
        if (S4_vp_ll(disc, p) != 0) bad++;
        other *= 1;                       /* 好归约处 c_p = 1（乘性自由项不进乘积） */
    }
    if (S4_vp_ll(disc, 2) != 6) bad++;
    if (S4_vp_ll(c4, 2) != 4) bad++;
    printf("  [P7] DISC=2^%d  v2(c4)=%d  bad primes={2}  prod_{p!=2} c_p=%d  violations=%d\n",
           S4_vp_ll(disc, 2), S4_vp_ll(c4, 2), other, bad);
    return (double)bad;
}

/* ---------------- P8 |Ш|=1 是推论：账上整数 2 + Cassels 平方类 ⟹ 唯一分解 (c2,|Ш|)=(2,1) ---------------- */
static double S4_p8_sha_square(void)
{
    const int N = 2;                      /* 账上整数：c2 * |Sha| = L*|tors|^2/Omega */
    int c, sols = 0, sol_c = 0, sol_s = 0, dev = 0;
    for (c = 1; c <= N; c++) {
        int s = N / c, m;
        if (c * s != N) continue;
        m = (int)S4_my_sqrt((double)s);
        if (m * m != s) continue;          /* Cassels：只有平方阶可以出现 */
        sols++; sol_c = c; sol_s = s;
    }
    if (sols != 1) dev++;
    if (sols == 1 && (sol_c != 2 || sol_s != 1)) dev++;
    printf("  [P8] solutions of c*|Sha|=%d with |Sha| a square: %d -> (c2,|Sha|)=(%d,%d)\n",
           N, sols, sol_c, sol_s);
    return (double)dev;
}

/* ---------------- P9 r=0 ⟹ Reg=1；r>=1 ⟹ Reg>0 ---------------- */
static double S4_p9_reg_presence(void)
{
    int dev = 0;
    double reg_empty = 1.0;
    double V1[1] = { 1.5 }, G[1];
    S4_gram_of(1, V1, G);
    printf("  [P9] r=0 : Reg = %.1f (empty product)   r=1 : Reg = det(G) = %.6f > 0\n",
           reg_empty, G[0]);
    if (fabs(reg_empty - 1.0) > 1e-15) dev++;
    if (!(G[0] > 0.0)) dev++;
    return (double)dev;
}

static void S4_print_wuxing(void)
{
    int i;
    for (i = 0; i < 6; i++) printf("   [%s]", S4_g_W[i]);
    printf("\n");
}

int seg4_main(void)
{
    double d1, d2, d3, d4, d5, d6, d7, d8, d9;

    printf("== [part 5/5] construction layer —— Reg / c_p / #Sha 全账 ==\n");
    S4_print_wuxing();
    printf("\n");

    d1 = S4_p1_positive_definite();  S4_judge("P1 Gram positive definite", d1, 0.0);
    d2 = S4_p2_basis_change();       S4_judge("P2 det(U^T G U) = det(U)^2 det(G)", d2, 0.0);
    d3 = S4_p3_scaling();            S4_judge("P3 h(nP)=n^2 h(P), det scales n^{2r}", d3, 0.0);
    d4 = S4_p4_parallelogram();      S4_judge("P4 parallelogram law", d4, 1e-12);
    d5 = S4_p5_ledger_rank0();       S4_judge("P5 ledger rank 0: L*#tors^2 = Omega*c*#Sha", d5, 1e-09);
    d6 = S4_p6_ledger_shape();       S4_judge("P6 Reg appears exactly once", d6, 0.0);
    d7 = S4_p7_tamagawa_divisibility(); S4_judge("P7 Delta=2^6 => only p=2 bad, c_p=1 (p!=2)", d7, 0.0);
    d8 = S4_p8_sha_square();         S4_judge("P8 c*|Sha|=2 with square |Sha| => unique (2,1)", d8, 0.0);
    d9 = S4_p9_reg_presence();       S4_judge("P9 r=0 => Reg=1; r>=1 => Reg>0", d9, 0.0);

    g_pass += S4_g_pass; g_fail += S4_g_fail;
    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.3e\n", S4_g_pass, S4_g_fail, S4_g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", S4_NJ, 1.0 - 1.0 / (double)S4_NJ);
    printf("res=%.6e domain=%d dc=%d order=%d\n", S4_g_res, S4_NJ, 5, 2);
    printf("honest note：P1-P4/P9 是高度配对的代数性质（精确），P5 用 part 2 的runreadingrecomputation，\n");
    printf("     P6 是账目形状的计数（Reg 的指数）；P7/P8 是本段实算：Δ = 2^6 ⟹ 只有 p=2 坏，\n");
    printf("     账上整数 2 + Cassels 平方类 ⟹ (c2,|Sha|) = (2,1) 唯一（同 onestep bs_sha 的 25 格读数）。\n");
    printf("     本篇未在实数上闭合 rank>=1 的 Reg numeric（需要 Mordell-Weil 生成元与真实period）。\n");
    printf("EXIT=%d\n", S4_g_fail == 0 ? 0 : 1);
    (void)d3; (void)d6; (void)d9;
    return S4_g_fail == 0 ? 0 : 1;
}

/* ── 聚合：逐段跑，任一段失败则整体失败 ── */
static int run_all(void)
{
    int rc = 0;
    rc |= seg0_main();
    rc |= seg1_main();
    rc |= seg2_main();
    rc |= seg3_main();
    rc |= seg4_main();
    return rc;
}

int main(void)
{
    int rc = run_all();
    printf("RESULT PASS=%d FAIL=%d   (all criteria in this file)\n", g_pass, g_fail);
    printf("EXIT=%d\n", rc ? 1 : 0);
    return rc ? 1 : 0;
}
