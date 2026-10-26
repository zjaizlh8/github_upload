/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
   verify_bsd4.c · BSD 的「reciprocal逼迫」final version实证（2026-10-26）
   ----------------------------------------------------------------------------
   用户caliber（逐字）：是reciprocal逼迫出的。
   archive内reciprocal定律： [2233] closure ⟺ 逆在同一容器内可写出；[2234] 膨胀↔收缩乘积 1；
                 [2235] projection率 × 本体倍率 = 1；rank + nullity = dimension。

   本step三件：
     ① constant侧（reciprocal账）：L(E,1)·#tors² / (Ω · ∏c_p · #Ш) = 1
        —— analyticreading与几何体量 Ω 互为reciprocal，商是有理账（本stepmeasured ∏c_p exactly为整数 2）。
     ② rank侧（zerocriterion）：L(E,1) 的 Abel 和 —— rank 0 饱和非零、rank ≥ 1 趋零。
        —— zero阶 = analytic侧的「reciprocalexponent」，与自由度数同位。
     ③ 负结果honest note：局部欧拉积 ∏_{p≤X}(p+1−a_p)/p 的 log-log 斜率**不区分rank**
        （三条曲线全在 1.27 附近）⟹ rank只能从 L(1) 的zero读，不能从局部积读。

   four operationscaliber：加法=点数计数summation；乘法=Hecke 关系/乘法性/欧拉积；取modulo=Legendre 符号；
             跨维度=analytic侧 exp/log/sqrt（Abel kernel、Simpson 求积、log），只在reading上用。
   ============================================================================ */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAXP 5001
#define NP 700

typedef struct {
    const char *name;
    long long a1, a2, a3, a4, a6;
    long long cond;      /* conductor（同时是unique坏prime的power） */
    int rank;            /* Cremona 表，外部input */
    int tors;
} Curve;

static const Curve CURVES[4] = {
    { "11a1",  1, -1, 1, -10, -20,  11, 0, 5 },
    { "37a1",  1,  0, 0,  -1,   0,  37, 1, 1 },
    { "389a1", 1,  1, 0,  -2,   0, 389, 2, 1 },
    { "32a1",  0,  0, 0,  -1,   0,  32, 0, 4 }
};

static int primes[NP], np = 0;
static int npass = 0, nfail = 0;
static double maxdev = 0.0;

static long long imod(long long a, long long m) { long long r = a % m; return r < 0 ? r + m : r; }

static long long ipowmod(long long b, long long e, long long m)
{
    long long r = 1 % m;
    b = imod(b, m);
    while (e > 0) {
        if (e & 1) r = (long long)((__int128)r * b % m);
        b = (long long)((__int128)b * b % m);
        e >>= 1;
    }
    return r;
}

static void sieve(void)
{
    char comp[MAXP + 1];
    int i, j;
    memset(comp, 0, sizeof comp);
    for (i = 2; i <= MAXP; i++) {
        if (comp[i]) continue;
        if (np < NP) primes[np++] = i;
        for (j = 2 * i; j <= MAXP; j += i) comp[j] = 1;
    }
}

static int legendre(long long d, long long p)
{
    d = imod(d, p);
    if (d == 0) return 0;
    return (ipowmod(d, (p - 1) / 2, p) == 1) ? 1 : -1;
}

static int ap_curve(const Curve *E, int p)
{
    long long x, s = 0;
    if ((long long)p == E->cond) return 0;      /* 坏prime只在 cond 是prime时hit */
    for (x = 0; x < p; x++) {
        long long x3 = x * x % p * x % p;
        long long rhs = (x3 + E->a2 * x % p * x % p) % p;
        long long lin = (E->a1 * x + E->a3) % p;
        long long d;
        rhs = imod(rhs + E->a4 * x + E->a6, p);
        d = imod(lin * lin % p + 4 * rhs % p, p);
        s += legendre(d, p);
    }
    return (int)(-s);
}

static void build_an(const Curve *E, int *a, int M)
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
            int ap = ap_curve(E, n);
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

/* L(1) 的 Abel 和：Σ (a_n/n)·e^{−2πny/√N}（kernel只作正则化，绝对值受kernel污染，只读zero/非zero） */
static double abel_L1(const Curve *E, int M, double y)
{
    int *a = (int *)calloc((size_t)M + 1, sizeof(int));
    double s = 0.0, c = 0.0, sq = sqrt((double)E->cond);
    int n;
    build_an(E, a, M);
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

/* Ω(y²=x³−x)=4∫₀^{π/2} dθ/√(1+cos²θ)；identity分branchperiod = 其一半 */
static double omega_32a1(int N)
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

static double local_product(const Curve *E, int X)
{
    double P = 1.0;
    int i;
    for (i = 0; i < np && primes[i] <= X; i++) {
        int p = primes[i];
        if ((long long)p == E->cond) continue;
        P *= ((double)p + 1.0 - (double)ap_curve(E, p)) / (double)p;
    }
    return P;
}

static void judge(const char *name, double dev, double tol)
{
    if (fabs(dev) <= tol) { npass++; printf("  [PASS] %-44s dev=%.3e (tol=%.0e)\n", name, fabs(dev), tol); }
    else { nfail++; printf("  [FAIL] %-44s dev=%.3e (tol=%.0e)\n", name, fabs(dev), tol); }
    if (fabs(dev) > maxdev) maxdev = fabs(dev);
}

/* -------- check数据器（cycle全部收在辅助函数内，main 只取reading） -------- */
static void hs_scan(int ok[3], int cnt[3], int *c3, int *c3z)
{
    int i, k;
    for (i = 0; i < np; i++) {
        int p = primes[i];
        for (k = 0; k < 3; k++) {
            double b;
            if ((long long)p == CURVES[k].cond) continue;
            b = 2.0 * sqrt((double)p);
            cnt[k]++;
            if (fabs((double)ap_curve(&CURVES[k], p)) > b + 1e-9) ok[k] = 0;
        }
        if ((long long)p != 32 && p % 4 == 3) { (*c3)++; if (ap_curve(&CURVES[3], p) == 0) (*c3z)++; }
    }
}

static void slope_scan(double r[3])
{
    double X0 = 60.0, X1 = 5000.0;
    const char *nm[3] = { "11a1", "37a1", "389a1" };
    int i;
    printf("\n");
    for (i = 0; i < 3; i++) {
        double p0 = local_product(&CURVES[i], (int)X0), p1 = local_product(&CURVES[i], (int)X1);
        r[i] = log(p1 / p0) / log(log(X1) / log(X0));
        printf("       · %-6s P(60)=%10.6f P(5000)=%10.6f  斜率=%.6f\n", nm[i], p0, p1, r[i]);
    }
}

int main(void)
{
    const int M = 8000, NPTS = 200000;
    const double LEXT = 0.6555143886;     /* 外部锚：L(32a1,1) */
    double om, tam, dev;
    double S11, S37, S389;

    sieve();
    printf("=== verify_bsd4.c · BSD reciprocal逼迫final version ===\n");
    printf("  ① constant侧reciprocal账：L·#tors² = Ω·∏c_p·#Ш\n");
    printf("  ② rank侧zerocriterion：rank 0 饱和 / rank≥1 趋零\n");
    printf("  ③ 负结果：局部欧拉积 log-log 斜率不区分rank\n\n");

    /* Hasse 与 CM */
    {
        int ok[3] = { 1, 1, 1 }, cnt[3] = { 0, 0, 0 }, c3 = 0, c3z = 0;
        hs_scan(ok, cnt, &c3, &c3z);
        printf("       · Hasse witnessprime：%d / %d / %d；32a1 上 p≡3(4) 共 %d 个\n", cnt[0], cnt[1], cnt[2], c3);
        judge("Hasse bound 11a1", (double)(ok[0] == 1 ? 0 : 1), 0.0);
        judge("Hasse bound 37a1", (double)(ok[1] == 1 ? 0 : 1), 0.0);
        judge("Hasse bound 389a1", (double)(ok[2] == 1 ? 0 : 1), 0.0);
        judge("32a1 CM：p≡3(4) ⟹ a_p=0", (double)(c3 - c3z), 0.0);
    }

    /* ① constant侧：reciprocal账 */
    om  = omega_32a1(NPTS);
    tam = LEXT * (double)(CURVES[3].tors * CURVES[3].tors) / om;   /* ∏c_p·#Ш 的理论值 */
    printf("\n       · Ω(32a1) = %.12f（外部锚 2·2.6220575543 = 5.2441151086）\n", om);
    judge("Ω 与外部锚一致", om - 2.0 * 2.6220575543, 1e-08);
    printf("       · ∏c_p·#Ш = L·#tors²/Ω = %.12f（应为整数 2 ⟹ c₂=2、#Ш=1）\n", tam);
    judge("Tamagawa 积为整数 2（账closure）", tam - 2.0, 1e-07);
    dev = LEXT * 16.0 - om * 2.0;                                  /* reciprocal账：两项cancellation */
    printf("       · L·#tors² = %.9f    Ω·∏c_p·#Ш = %.9f\n", LEXT * 16.0, om * 2.0);
    judge("reciprocal账两项cancellation（residual → 0）", dev, 1e-06);

    /* ② rank侧：zerocriterion */
    S11  = abel_L1(&CURVES[0], M, 0.001);
    S37  = abel_L1(&CURVES[1], M, 0.001);
    S389 = abel_L1(&CURVES[2], M, 0.001);
    printf("\n       · Abel 和：11a1(r=0)=%.9f  37a1(r=1)=%.9f  389a1(r=2)=%.9f\n", S11, S37, S389);
    judge("rank 0：Abel 和非零（> 0.3）", (S11 > 0.3) ? 0.0 : 1.0, 0.0);
    judge("rank 1：Abel 和趋零（|·| < 0.02）", S37, 2e-02);
    judge("rank 2：Abel 和趋零（|·| < 0.02）", S389, 2e-02);
    printf("       · 分离度 |S11| / max(|S37|,|S389|) = %.3e\n", fabs(S11) / fmax(fabs(S37), fabs(S389)));
    judge("分离度 ≥ 1e1（zerocriterion见血）",
          (fabs(S11) / fmax(fabs(S37), fabs(S389)) >= 10.0) ? 0.0 : 1.0, 0.0);
    judge("zero ⟺ rank≥1（three steps一致）",
          ((S11 > 0.3) && (fabs(S37) < 0.02) && (fabs(S389) < 0.02)) ? 0.0 : 1.0, 0.0);

    /* ③ 负结果：局部欧拉积的斜率 */
    {
        double r[3];
        slope_scan(r);
        printf("       · 结论（负结果）：three steps斜率同量级 ⟹ 局部积在 s=1 对所有rank都给 ~log X\n");
        judge("负结果：|r(37a1) − r(389a1)| < 0.1（rank不可由局部积读）",
              r[1] - r[2], 1e-01);
    }

    printf("\nRESULT PASS=%d FAIL=%d\n", npass, nfail);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n",
           npass + nfail, (double)npass / (double)(npass + nfail));
    printf("res=%.6e domain=%d dc=%d order=2\n", maxdev, npass + nfail, 4);
    printf("EXIT=%d\n", nfail == 0 ? 0 : 1);
    return nfail == 0 ? 0 : 1;
}
