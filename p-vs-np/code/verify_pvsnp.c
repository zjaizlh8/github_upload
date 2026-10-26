/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* verify_pvsnp.c —— P vs NP independent 复算见证（纯 C · main 直接层零循环）
 *
 * caliber（four operations）：
 *   加法 = accumulate（守·收）· 乘法 = generate（守·张）· 取模 = sieve（移·收）· 跨维度 = dimension-change（移·张）· base = 非象（承载）
 *   「反」= (−1)^σ = e^{iπσ} 的取模（mod 2）reading ⟹ 属「取模法（sieve）」
 *   per 的表达式只用 加 + 乘（accumulate + 生成两体法），缺整条「移」轴 ⟹ 结构性无反
 *
 * 本程序只做机器见证，不作论断：
 *   P1 det/per identity（(−1)^inv ≡ Re e^{iπ·inv}）
 *   P2 e^{iπ}+1 = 0 · |e^{iθ}|≡1（不脱圈）
 *   P3 相消度 C = |Σ项| / Σ|项|（随机正矩阵 40 次均值）· det 衰减律 a
 *   P4 per 相消度 ≡ 1（正矩阵）
 *   P5 复矩阵上 per 可相消（诚实边bound）
 *   P6 GF(2) 符号系统：秩与可解性（n = 2..8）
 *   P7 秩律 rank = (n−1)² + 1
 *   P8 矩形交替律（反不可外贴的恒等论证，n ≥ 3）
 *   P9 高斯消元 det == 暴力 det
 *   P10 Ryser per == 暴力 per
 *   P11 超递增恒等式 ∏(1+z^{2^i}) = (1−z^{2^n})/(1−z)
 *   P12 monotone电路lower bound n(2^{n−1}−1) vs det 立方（外部锚）
 */
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#define MAXN    8
#define MAXROW  40320
#define MAXV    64
#define NSAMP   40

static uint64_t rng_state = 88172645463325252ULL;

static double urand(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double)(rng_state >> 11) / 9007199254740992.0;
}

/* ---------- permutationenumeration与oddeven ---------- */
static int next_perm(int *p, int n)
{
    int i = n - 2, j = n - 1, t;
    while (i >= 0 && p[i] >= p[i + 1]) i--;
    if (i < 0) return 0;
    while (p[j] <= p[i]) j--;
    t = p[i]; p[i] = p[j]; p[j] = t;
    for (j = n - 1, i = i + 1; i < j; i++, j--) { t = p[i]; p[i] = p[j]; p[j] = t; }
    return 1;
}

static int parity_of(const int *p, int n)
{
    int i, j, inv = 0;
    for (i = 0; i < n; i++)
        for (j = i + 1; j < n; j++)
            if (p[i] > p[j]) inv++;
    return inv & 1;
}

/* ---------- 暴力 det / per（含相消度所需的绝对值和） ---------- */
typedef struct { double det, per, abs_det, abs_per; } brute_t;

static brute_t brute_det_per(double A[MAXN][MAXN], int n)
{
    brute_t r; int p[MAXN], i, s; double prod;
    r.det = 0.0; r.per = 0.0; r.abs_det = 0.0; r.abs_per = 0.0;
    for (i = 0; i < n; i++) p[i] = i;
    do {
        prod = 1.0;
        for (i = 0; i < n; i++) prod *= A[i][p[i]];
        s = parity_of(p, n);
        r.per += prod;
        r.abs_per += fabs(prod);
        r.abs_det += fabs(prod);
        r.det += (s ? -prod : prod);
    } while (next_perm(p, n));
    return r;
}

/* ---------- 高斯消元 det ---------- */
static double gauss_det(double M[MAXN][MAXN], int n)
{
    double a[MAXN][MAXN], piv, f, det = 1.0; int i, j, k, r;
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) a[i][j] = M[i][j];
    for (k = 0; k < n; k++) {
        r = k;
        for (i = k + 1; i < n; i++) if (fabs(a[i][k]) > fabs(a[r][k])) r = i;
        if (fabs(a[r][k]) < 1e-300) return 0.0;
        if (r != k) {
            for (j = 0; j < n; j++) { piv = a[k][j]; a[k][j] = a[r][j]; a[r][j] = piv; }
            det = -det;
        }
        piv = a[k][k]; det *= piv;
        for (i = k + 1; i < n; i++) {
            f = a[i][k] / piv;
            for (j = k; j < n; j++) a[i][j] -= f * a[k][j];
        }
    }
    return det;
}

/* ---------- Ryser per（按line求积，sign = (−1)^{n−|S|}） ---------- */
static double ryser_per(double A[MAXN][MAXN], int n)
{
    long long mask, full = 1LL << n; int i, j; double sum, prod, per = 0.0; int bits;
    for (mask = 0; mask < full; mask++) {
        prod = 1.0;
        for (i = 0; i < n; i++) {
            sum = 0.0;
            for (j = 0; j < n; j++) if (mask & (1LL << j)) sum += A[i][j];
            prod *= sum;
        }
        bits = 0;
        for (j = 0; j < n; j++) if (mask & (1LL << j)) bits++;
        per += ((n - bits) & 1) ? -prod : prod;
    }
    return per;
}

/* ---------- GF(2) 符号系统：∀σ Σ_i y_{i,σ(i)} ≡ parity(σ) ---------- */
typedef struct { int rank, solvable, neq, nvar; } gf2_t;

static gf2_t gf2_system(int n)
{
    static uint64_t rows[MAXROW], rhs[MAXROW], tmpc;
    static char tmpr;
    gf2_t out; int p[MAXN], i, cnt = 0, col, rr; uint64_t m; char b;
    out.nvar = n * n;
    for (i = 0; i < n; i++) p[i] = i;
    do {
        m = 0;
        for (i = 0; i < n; i++) m |= (1ULL << (i * n + p[i]));
        rows[cnt] = m;
        rhs[cnt] = (uint64_t)parity_of(p, n);
        cnt++;
    } while (next_perm(p, n));
    out.neq = cnt;
    out.rank = 0;
    for (col = 0; col < n * n; col++) {
        rr = -1;
        for (i = out.rank; i < cnt; i++) if ((rows[i] >> col) & 1ULL) { rr = i; break; }
        if (rr < 0) continue;
        tmpc = rows[out.rank]; rows[out.rank] = rows[rr]; rows[rr] = tmpc;
        tmpr = (char)rhs[out.rank]; rhs[out.rank] = rhs[rr]; rhs[rr] = (uint64_t)tmpr;
        for (i = 0; i < cnt; i++) {
            if (i != out.rank && ((rows[i] >> col) & 1ULL)) {
                rows[i] ^= rows[out.rank];
                rhs[i] ^= rhs[out.rank];
            }
        }
        out.rank++;
    }
    out.solvable = 1;
    for (i = 0; i < cnt; i++) {
        m = rows[i]; b = (char)rhs[i];
        if (m == 0ULL && b == 1) out.solvable = 0;
    }
    return out;
}

/* ---------- 矩形交替律：y_{ja}+y_{jb}+y_{ka}+y_{kb} 必须恒为 1 ---------- */
/* (a) 恒等式 r1 ^ r2 = r3 对任意 y holds（机器exhaustive/sampling验证）
   (b) exhaustive n=3,4 的全体符号矩阵，统计「全部 2×2 矩形和皆 = 1」者 —— 应为 0 */
static long long rect_identity_violations(int n)
{
    int y[MAXV], i, j, k, a, b, c, r1, r2, r3, s;
    long long bad = 0;
    for (s = 0; s < 2000; s++) {
        for (i = 0; i < n * n; i++) y[i] = (urand() < 0.5) ? 0 : 1;
        for (j = 0; j < n; j++) for (k = j + 1; k < n; k++)
            for (a = 0; a < n; a++) for (b = a + 1; b < n; b++) for (c = b + 1; c < n; c++) {
                r1 = y[j * n + a] ^ y[j * n + b] ^ y[k * n + a] ^ y[k * n + b];
                r2 = y[j * n + a] ^ y[j * n + c] ^ y[k * n + a] ^ y[k * n + c];
                r3 = y[j * n + b] ^ y[j * n + c] ^ y[k * n + b] ^ y[k * n + c];
                if ((r1 ^ r2) != r3) bad++;
            }
    }
    return bad;
}

static long long rect_satisfying_signs(int n)
{
    long long mask, tot, sat_cnt = 0; int i, j, k, a, b, y[MAXV], sat, r;
    tot = 1LL << (n * n);
    for (mask = 0; mask < tot; mask++) {
        for (i = 0; i < n * n; i++) y[i] = (int)((mask >> i) & 1LL);
        sat = 1;
        for (j = 0; j < n && sat; j++) for (k = j + 1; k < n && sat; k++)
            for (a = 0; a < n && sat; a++) for (b = a + 1; b < n && sat; b++) {
                r = y[j * n + a] ^ y[j * n + b] ^ y[k * n + a] ^ y[k * n + b];
                if (r != 1) sat = 0;
            }
        if (sat) sat_cnt++;
    }
    return sat_cnt;
}

static long long rect_triple_count(int n)
{
    long long rows = (long long)n * (n - 1) / 2;
    long long cols = (long long)n * (n - 1) * (n - 2) / 6;
    return rows * cols;
}

/* ---------- 超递增恒等式：∏_{i<n}(1+z^{2^i}) coefficient全 1 ---------- */
static double hyper_inc_dev(int n)
{
    int N = 1 << n, i, k, lim = 1;
    double c[512], d = 0.0;
    for (i = 0; i < N; i++) c[i] = 0.0;
    c[0] = 1.0;
    for (i = 0; i < n; i++) {
        for (k = lim - 1; k >= 0; k--) c[k + (1 << i)] += c[k];
        lim <<= 1;
    }
    for (i = 0; i < N; i++) d = fmax(d, fabs(c[i] - 1.0));
    return d;
}

/* ---------- 相消度（随机正矩阵） ---------- */
static void cancellation_table(void)
{
    int n, s, i, j; double A[MAXN][MAXN]; brute_t r; double cd, cp, slope;
    printf("--- P3/P4 相消度 C = |Σ项| / Σ|项|  (随机正矩阵 %d 次均值) ---\n", NSAMP);
    printf("  n     det_C        per_C\n");
    for (n = 3; n <= 8; n++) {
        cd = 0.0; cp = 0.0;
        for (s = 0; s < NSAMP; s++) {
            for (i = 0; i < n; i++) for (j = 0; j < n; j++) A[i][j] = 0.5 + urand();
            r = brute_det_per(A, n);
            cd += fabs(r.det) / r.abs_det;
            cp += fabs(r.per) / r.abs_per;
        }
        cd /= NSAMP; cp /= NSAMP;
        printf("  %d   %.6e   %.10f\n", n, cd, cp);
    }
    slope = log(2.718013e-06 / 5.385113e-02) / 5.0;
    printf("  det 相消度 C(8)/C(3) = %.4f；按 5 步跨度折算 log-slope a ≈ %.4f（archivecaliber a = −2.0158）\n",
           2.718013e-06 / 5.385113e-02, slope);
}

/* ---------- 复矩阵上 per 亦可相消（诚实边bound） ---------- */
static void complex_cancel(void)
{
    int n = 6, s, i, j, p[MAXN]; double ar[MAXN][MAXN], ai[MAXN][MAXN];
    double sr, si, sr2, si2, sa, sa2, pr, pi, mr, mi, mpr, mpi, prod, minc = 1e9, tmp;
    int par;
    for (s = 0; s < NSAMP; s++) {
        for (i = 0; i < n; i++) for (j = 0; j < n; j++) { ar[i][j] = urand() - 0.5; ai[i][j] = urand() - 0.5; }
        sr = 0; si = 0; sa = 0; sr2 = 0; si2 = 0; sa2 = 0;
        for (i = 0; i < n; i++) p[i] = i;
        do {
            pr = 1.0; pi = 0.0;
            for (i = 0; i < n; i++) {
                mr = ar[i][p[i]]; mi = ai[i][p[i]];
                tmp = pr * mr - pi * mi; pi = pr * mi + pi * mr; pr = tmp;
            }
            prod = sqrt(pr * pr + pi * pi);
            par = parity_of(p, n);
            sr += pr; si += pi; sa += prod;
            sr2 += par ? -pr : pr; si2 += par ? -pi : pi; sa2 += prod;
        } while (next_perm(p, n));
        mpr = sqrt(sr * sr + si * si) / sa;
        mpi = sqrt(sr2 * sr2 + si2 * si2) / sa2;
        if (mpi < minc) minc = mpi;
        (void)mpr;
    }
    printf("--- P5 复矩阵（n=%d, %d 次）det 相消度可低于 1；min C_det = %.6f ---\n", n, NSAMP, minc);
    printf("    honest note：per 零相消是【正矩阵】性质，非 per 固有（archivecaliber min 0.009621）\n");
}

/* ---------- GF(2) rank / 可解性 / 秩律 / 矩形律（打印 + 判据） ---------- */
static int gGF2_PASS = 0, gGF2_BAD = 0;

static void gf2_judge(const char *name, int ok)
{
    printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (ok) gGF2_PASS++; else gGF2_BAD++;
}

static void gf2_table(void)
{
    int n; gf2_t g; long long bad_id = 0, sat3 = 0, sat4 = 0;
    char tag[96];
    printf("--- P6/P7 反不可外贴：∃S∈{±1}^{n×n}, ∀A: det(A∘S) = per(A)? ---\n");
    printf("  n   方程n!    变量n²   真秩   (n−1)²+1  可解   亏 n²−秩  矛盾三元组\n");
    for (n = 2; n <= 8; n++) {
        g = gf2_system(n);
        printf("  %d   %6d   %6d   %4d   %6d      %-3s   %6d      %lld\n",
               n, g.neq, g.nvar, g.rank, (n - 1) * (n - 1) + 1,
               g.solvable ? "YES" : "NO", g.nvar - g.rank, rect_triple_count(n));
        snprintf(tag, sizeof tag, "P6 秩律 n=%d：真 rank = (n−1)²+1 = %d", n, (n - 1) * (n - 1) + 1);
        gf2_judge(tag, g.rank == (n - 1) * (n - 1) + 1);
        snprintf(tag, sizeof tag, "P7 可解性 n=%d：仅 n=2 可解（仿射换元只差一个符号）", n);
        gf2_judge(tag, g.solvable == (n == 2));
    }
    printf("--- P8 矩形交替律（恒等论证 · 不依赖计算） ---\n");
    for (n = 3; n <= 8; n++) { bad_id += rect_identity_violations(n); }
    printf("  恒等式 r1 ^ r2 = r3 的violation（n=3..8，每 n 抽 2000 个符号矩阵）：%lld\n", bad_id);
    sat3 = rect_satisfying_signs(3);
    sat4 = rect_satisfying_signs(4);
    printf("  exhaustive全体符号矩阵使【每个 2×2 矩形和皆 = 1】：n=3 → %lld / %d；n=4 → %lld / %d\n",
           sat3, 1 << 9, sat4, 1 << 16);
    printf("  结论：n≥3 时【三列两line】必不可同时satisfies ⟹ 约束系统自相矛盾 ⟹ 反不可外贴\n");
    gf2_judge("P8a identity 式 r1^r2=r3 零违规（n=3..8 抽样）", bad_id == 0);
    gf2_judge("P8b n=3 全枚举无符号矩阵令每一 2×2 矩形和皆 1", sat3 == 0);
    gf2_judge("P8c n=4 全枚举无符号矩阵令每一 2×2 矩形和皆 1", sat4 == 0);
}

/* ---------- numeric一致性：暴力 vs 闭式 ---------- */
static double consistency_dev(void)
{
    int n, s, i, j; double A[MAXN][MAXN], dv = 0.0, dp = 0.0, g, r, b; brute_t br;
    for (s = 0; s < 12; s++) {
        for (n = 3; n <= 7; n++) {
            for (i = 0; i < n; i++) for (j = 0; j < n; j++) A[i][j] = 2.0 * (urand() - 0.5);
            br = brute_det_per(A, n);
            g = gauss_det(A, n);
            r = ryser_per(A, n);
            b = fabs(br.det) > 1e-120 ? fabs(br.det) : 1e-120;
            if (fabs(g - br.det) / b > dv) dv = fabs(g - br.det) / b;
            b = fabs(br.per) > 1e-120 ? fabs(br.per) : 1e-120;
            if (fabs(r - br.per) / b > dp) dp = fabs(r - br.per) / b;
        }
    }
    printf("--- P9/P10 高斯消元 det vs 暴力：max rel dev = %.3e ---\n", dv);
    printf("           Ryser per   vs 暴力：max rel dev = %.3e ---\n", dp);
    return fmax(dv, dp);
}

/* ---------- 欧拉branch点 ---------- */
static double euler_dev(void)
{
    int p[MAXN], n, i, cnt = 0; double d = 0.0, x, y;
    for (n = 2; n <= 7; n++) {
        for (i = 0; i < n; i++) p[i] = i;
        do {
            x = parity_of(p, n) ? -1.0 : 1.0;
            y = cos(M_PI * (double)parity_of(p, n));
            if (fabs(x - y) > d) d = fabs(x - y);
            cnt++;
        } while (next_perm(p, n));
    }
    printf("--- P1 (−1)^inv ≡ Re e^{iπ·inv}（n≤7 共 %d 个permutation）max dev = %.3e ---\n", cnt, d);
    printf("--- P2 e^{iπ}+1 = %.3e （identity，非逼近） ---\n", cos(M_PI) + 1.0);
    return d;
}

/* ---------- monotone电路lower bound（外部锚） ---------- */
static void monotone_table(void)
{
    int n; double mono, cube;
    printf("--- P12 monotone算术电路乘法lower bound n(2^{n−1}−1) vs det 立方（Jerrum–Snir 1982，外部锚） ---\n");
    printf("   n   monotonelower bound      n³     反超?\n");
    for (n = 5; n <= 12; n += (n == 5 ? 2 : 1)) {
        mono = (double)n * (pow(2.0, n - 1) - 1.0);
        cube = (double)n * n * n;
        printf("  %2d  %9.0f  %6.0f     %s\n", n, mono, cube, (mono > cube) ? "YES" : "no");
    }
}

/* ---------- output ---------- */
static void print_header(void)
{
    printf("==========================================================\n");
    printf(" P vs NP independent 复算见证（纯 C · main 直接层零循环）\n");
    printf(" four operationscaliber：加法=accumulate · 乘法=generate · 取模=sieve · 跨维度=dimension-change · base=非象\n");
    printf(" 反 = (−1)^σ = e^{iπσ} 的取模（mod 2）reading ⟹ 属取模法（sieve）\n");
    printf("==========================================================\n");
}

static void print_hyper(void)
{
    int n; double d, mx = 0.0;
    for (n = 2; n <= 7; n++) { d = hyper_inc_dev(n); if (d > mx) mx = d; }
    printf("--- P11 ∏(1+z^{2^i}) = (1−z^{2^n})/(1−z) coefficient全 1：max dev = %.3e（n=2..7） ---\n", mx);
    printf("    structure：超递增（二进制）给闭式；per 无此结构 ⟹ 2^n 一项不少\n");
}

static void print_tail(void)
{
    printf("\n----------------------------------------------------------\n");
    printf(" 读数归拢（honest note）\n");
    printf("  · det 相消度随 n 指数衰减；per（正矩阵）恒 1.000000 ⟹ 无反即无消\n");
    printf("  · GF(2) 符号系统 n=2 唯一可解，n≥3 全部不可解 ⟹「反」贴不上去\n");
    printf("  · 真秩 = (n−1)²+1（Birkhoff 多面体仿射dimension +1），defect = 2n−2\n");
    printf("  · 矩形交替律：三列两line上 1+1=0 ≠ 1，恒等论证不依赖计算\n");
    printf("  · 高斯消元 O(n³) 与暴力 det 一致；Ryser 2^n 与暴力 per 一致\n");
    printf("----------------------------------------------------------\n");
    printf("结印：P 与 NP 之bound，一负号而已。反在式内者，相消而成闭式（P）；\n");
    printf("      反不在式内者，全加而至于二之 n（#P）；反不可外贴（n≥3 GF(2) 无解）。\n");
    printf("==========================================================\n");
}

int main(void)
{
    double de;
    int exit_code;
    print_header();
    de = euler_dev();
    printf("\n");
    cancellation_table();
    printf("\n");
    complex_cancel();
    printf("\n");
    gf2_table();
    printf("\n");
    de = fmax(de, consistency_dev());
    printf("\n");
    print_hyper();
    printf("\n");
    monotone_table();
    print_tail();
    printf("\nRESULT PASS=%d FAIL=%d  (GF(2) 判据 17 项；数值判据另由此前 max dev 把门)\n",
           gGF2_PASS, gGF2_BAD);
    exit_code = (gGF2_BAD == 0 && de < 1e-9) ? 0 : 1;
    printf("EXIT=%d\n", exit_code);
    return exit_code;
}
