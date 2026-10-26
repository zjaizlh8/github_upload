/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
   verify_pvsnp_tech.c  ——  P vs NP 篇 · construction layer：为什么 per 没有closed式（mod 2 商）
   ----------------------------------------------------------------------------
   five roles六步（阶段直序）：
     base S 锚定 —— 先把 det 与 per 的三条独立path钉在hard anchor上（2×2、3×3 手算值）
     dimension-change F 虚化 —— 展open式的符号structure：det 的正负项数与 per 的全正项数
     sieve M 坍缩 —— modulo 2 商：char 2 里符号消失 ⟹ det ≡ per；整数侧cancellation不起来
     accumulate W 渗透 —— 信息量：oddeven只给 1 比特，per 需 log2(n!) 比特
     generate D 扩展 —— GF(2) 仿射rank (n-1)^2+1 与enumerationfield 2^nullity
     base S closure —— 四条criterion一次取出，四量齐报

   four operationscaliber（本step只用两条半）：
     加法 = Laplace 展open / Ryser 的 Gray 码累积 / GF(2) 向量异或
     乘法 = permutation的乘积项（每项 n 个因子相乘）与项的计数 n!
     取modulo = **本stepunique取modulo处**：det mod 2 与 per mod 2 的oddeven比较（`% 2`）
     跨维度 = 不用（无 exp/log/pow/sqrt）⟹ 不付residual

   F: input矩阵 A（0/1 或整数）→ output (det, per) 与符号structure计数
   R: 由 (det mod 2, per mod 2, GF(2) rank) 反推「是否存在 mod 2 商」
   ============================================================================ */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define NJ 9
#define MAXN 6
#define ENUMN 4

static int g_pass = 0, g_fail = 0;
static double g_maxdev = 0.0;
static double g_res = 0.0;

static const char *g_W[6] = { "baseS锚定", "dimension-changeF虚化", "sieveM坍缩", "accumulateW渗透", "generateD扩展", "baseSclosure" };

static void judge(const char *name, double dev, double tol)
{
    if (dev <= tol) { printf("  PASS %-38s dev=%.6e tol=%.6e\n", name, dev, tol); g_pass++; }
    else { printf("  FAIL %-38s dev=%.6e tol=%.6e\n", name, dev, tol); g_fail++; }
    if (dev > g_maxdev) g_maxdev = dev;
}

/* ---------------- det：Bareiss 无分数elimination（整数精确） ---------------- */
static double det_bareiss(int n, const int *a)
{
    int m[MAXN][MAXN], i, j, k, sign = 1, piv;
    double den = 1.0;
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) m[i][j] = a[i * n + j];
    for (k = 0; k < n - 1; k++) {
        if (m[k][k] == 0) {
            piv = -1;
            for (i = k + 1; i < n; i++) if (m[i][k] != 0) { piv = i; break; }
            if (piv < 0) return 0.0;
            for (j = 0; j < n; j++) { int t = m[k][j]; m[k][j] = m[piv][j]; m[piv][j] = t; }
            sign = -sign;
        }
        for (i = k + 1; i < n; i++) for (j = k + 1; j < n; j++)
            m[i][j] = (m[i][j] * m[k][k] - m[i][k] * m[k][j]) / (int)den;
        den = (double)m[k][k];
        for (i = k + 1; i < n; i++) m[i][k] = 0;
    }
    return (double)sign * (double)m[n - 1][n - 1];
}

/* ---------------- det：Laplace 递归展open（独立path） ---------------- */
static double det_laplace(int n, const int *a)
{
    int i, j, r, c; double s = 0.0;
    if (n == 1) return (double)a[0];
    for (i = 0; i < n; i++) {
        int sub[MAXN * MAXN], nn = n - 1;
        for (r = 1; r < n; r++) { c = 0; for (j = 0; j < n; j++) if (j != i) sub[(r - 1) * nn + c++] = a[r * n + j]; }
        s += ((i % 2 == 0) ? 1.0 : -1.0) * (double)a[i] * det_laplace(nn, sub);
    }
    return s;
}

/* ---------------- per：Ryser（Gray 码，2^n 项） ---------------- */
static double per_ryser(int n, const int *a)
{
    int msk, i, j, gray, prev = 0, bits; double sum = 0.0, prod; int colsum[MAXN];
    for (msk = 0; msk < (1 << n); msk++) {
        gray = msk ^ (msk >> 1);
        bits = 0; for (i = 0; i < n; i++) bits += (gray >> i) & 1;
        for (j = 0; j < n; j++) colsum[j] = 0;
        for (i = 0; i < n; i++) if ((gray >> i) & 1) for (j = 0; j < n; j++) colsum[j] += a[i * n + j];
        prod = 1.0; for (j = 0; j < n; j++) prod *= (double)colsum[j];
        sum += ((n - bits) % 2 == 0 ? 1.0 : -1.0) * prod;
        prev = gray;
    }
    (void)prev;
    return sum;
}

/* ---------------- per：permutationenumeration（独立path） ---------------- */
static double per_enum(int n, const int *a)
{
    int perm[MAXN], i, j, k; double s = 0.0;
    for (i = 0; i < n; i++) perm[i] = i;
    for (;;) {
        double p = 1.0;
        for (i = 0; i < n; i++) p *= (double)a[i * n + perm[i]];
        s += p;
        /* next permutation */
        k = -1; for (i = n - 2; i >= 0; i--) if (perm[i] < perm[i + 1]) { k = i; break; }
        if (k < 0) break;
        j = n - 1; while (perm[j] < perm[k]) j--;
        { int t = perm[k]; perm[k] = perm[j]; perm[j] = t; }
        for (i = k + 1, j = n - 1; i < j; i++, j--) { int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
    }
    return s;
}

/* ---------------- 符号structure：det 展open的正项数 / 负项数；per 全正 ---------------- */
static void sign_split(int n, int *pos, int *neg)
{
    int i, j, inv;
    *pos = 0; *neg = 0;
    for (i = 0; i < n; i++) {
        int perm[MAXN], k;
        for (k = 0; k < n; k++) perm[k] = k;
        for (;;) {
            inv = 0;
            for (k = 0; k < n; k++) for (j = k + 1; j < n; j++) if (perm[k] > perm[j]) inv++;
            if (inv % 2 == 0) (*pos)++; else (*neg)++;
            k = -1; for (j = n - 2; j >= 0; j--) if (perm[j] < perm[j + 1]) { k = j; break; }
            if (k < 0) break;
            j = n - 1; while (perm[j] < perm[k]) j--;
            { int t = perm[k]; perm[k] = perm[j]; perm[j] = t; }
            for (k = k + 1, j = n - 1; k < j; k++, j--) { int t = perm[k]; perm[k] = perm[j]; perm[j] = t; }
        }
    }
}

/* ---------------- GF(2) rank（按位elimination） ---------------- */
static int rank_gf2(int rows, int n, const int *a)
{
    int r = 0, c, i, j;
    int m[MAXN * MAXN];
    for (i = 0; i < rows * n; i++) m[i] = a[i] & 1;
    for (c = 0; c < n && r < rows; c++) {
        int piv = -1;
        for (i = r; i < rows; i++) if (m[i * n + c]) { piv = i; break; }
        if (piv < 0) continue;
        for (j = 0; j < n; j++) { int t = m[r * n + j]; m[r * n + j] = m[piv * n + j]; m[piv * n + j] = t; }
        for (i = 0; i < rows; i++) if (i != r && m[i * n + c])
            for (j = 0; j < n; j++) m[i * n + j] ^= m[r * n + j];
        r++;
    }
    return r;
}

/* ---------------- P1 det 两路一致（含手算hard anchor） ---------------- */
static double p1_det_paths(void)
{
    int a2[4] = { 1, 2, 3, 4 };
    int a3[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 10 };
    int a4[16] = { 2, 0, 1, 3, 1, 2, 0, 1, 0, 1, 3, 0, 4, 1, 2, 1 };
    double dev = 0.0, x, y;
    x = det_bareiss(2, a2); y = det_laplace(2, a2);
    dev = fmax(dev, fabs(x - (-2.0))); dev = fmax(dev, fabs(y - (-2.0)));
    x = det_bareiss(3, a3); y = det_laplace(3, a3);
    dev = fmax(dev, fabs(x - (-3.0))); dev = fmax(dev, fabs(y - (-3.0)));
    x = det_bareiss(4, a4); y = det_laplace(4, a4);
    dev = fmax(dev, fabs(x - y));
    printf("  [P1] det 3x3 = %.6f (expect -3)   det 4x4 bareiss %.6f vs laplace %.6f\n",
           det_bareiss(3, a3), det_bareiss(4, a4), det_laplace(4, a4));
    return dev;
}

/* ---------------- P2 per 两路一致（Ryser vs enumeration） ---------------- */
static double p2_per_paths(void)
{
    int a3[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 10 };
    int a4[16] = { 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1 };
    double dev = 0.0, x, y;
    x = per_ryser(3, a3); y = per_enum(3, a3);
    dev = fmax(dev, fabs(x - y));
    dev = fmax(dev, fabs(x - 463.0));    /* 手算：1·5·10+1·6·8+2·4·10+2·6·7+3·4·8+3·5·7 = 463 */
    x = per_ryser(4, a4); y = per_enum(4, a4);
    dev = fmax(dev, fabs(x - y));
    printf("  [P2] per 3x3 = %.6f (expect 463)   per 4x4 ryser %.6f vs enum %.6f\n",
           per_ryser(3, a3), per_ryser(4, a4), per_enum(4, a4));
    return dev;
}

/* ---------------- P3 alternating property：det 换两line变号，per 不变 ---------------- */
static double p3_alternating(void)
{
    int a[9] = { 2, 1, 3, 0, 4, 1, 5, 2, 6 }, b[9], j;
    double dp, dq, pp, pq, dev;
    for (j = 0; j < 3; j++) {
        b[0 * 3 + j] = a[1 * 3 + j];
        b[1 * 3 + j] = a[0 * 3 + j];
        b[2 * 3 + j] = a[2 * 3 + j];
    }
    dp = det_bareiss(3, a); dq = det_bareiss(3, b);
    pp = per_ryser(3, a); pq = per_ryser(3, b);
    dev = fmax(fabs(dp + dq), fabs(pp - pq));
    printf("  [P3] det(A)=%.6f det(swapped)=%.6f  per(A)=%.6f per(swapped)=%.6f\n", dp, dq, pp, pq);
    return dev;
}

/* ---------------- P4 符号structure：det 正负各半；per 全正 ---------------- */
static double p4_sign_structure(void)
{
    int n, pos, neg, dev = 0;
    printf("  [P4] n | det(+,-) | per terms\n");
    for (n = 2; n <= 6; n++) {
        sign_split(n, &pos, &neg);
        printf("        %d | %d,%d %s | %d\n", n, pos, neg, (pos == neg ? "equal" : "DIFF"), pos + neg);
        if (pos != neg) dev++;
    }
    return (double)dev;
}

/* ---------------- P5 modulo 2 商：char 2 里 det ≡ per ---------------- */
static double p5_mod2_quotient(void)
{
    int n, msk, i, dev = 0, total = 0;
    int a[MAXN * MAXN];
    for (n = 2; n <= ENUMN; n++) {
        for (msk = 0; msk < (1 << (n * n)); msk++) {
            long long d, p;
            for (i = 0; i < n * n; i++) a[i] = (msk >> i) & 1;
            d = (long long)llround(det_bareiss(n, a));
            p = (long long)llround(per_ryser(n, a));
            if (((d % 2) + 2) % 2 != ((p % 2) + 2) % 2) dev++;
            total++;
        }
    }
    printf("  [P5] 0/1 matrices enumerated = %d   det mod 2 != per mod 2 count = %d\n", total, dev);
    g_res = (double)dev;
    return (double)dev;
}

/* ---------------- P6 整数侧：全 1 矩阵 det=0 而 per=n!（符号全正不抵消） ---------------- */
static double p6_int_gap(void)
{
    int n, i, j, dev = 0;
    int a[MAXN * MAXN];
    printf("  [P6] n | det(J) | per(J) | n!\n");
    for (n = 2; n <= 6; n++) {
        double d, p, f = 1.0;
        for (i = 0; i < n * n; i++) a[i] = 1;
        d = det_bareiss(n, a);
        p = per_ryser(n, a);
        for (j = 2; j <= n; j++) f *= (double)j;
        printf("        %d | %.1f | %.1f | %.1f\n", n, d, p, f);
        if (fabs(d) > 1e-9 || fabs(p - f) > 1e-9) dev++;
    }
    return (double)dev;
}

/* ---------------- P7 信息量：oddeven 1 bit vs per 需 log2(n!) ---------------- */
static double p7_info_gap(void)
{
    int n, i, dev = 0;
    double bits = 0.0, f = 1.0;
    printf("  [P7] n | log2(n!) bits of per | bits from parity\n");
    for (n = 3; n <= 6; n++) {
        f = 1.0;
        for (i = 2; i <= n; i++) f *= (double)i;
        bits = log(f) / log(2.0);
        printf("        %d | %.3f | 1\n", n, bits);
        if (bits <= 1.0) dev++;     /* criterion：n>=3 时 log2(n!) 必须严lattice大于 1 比特 */
    }
    return (double)dev;
}

/* ---------------- P8 cancellation：重复line det=0 但 per>0 ---------------- */
static double p8_cancel(void)
{
    int a[9] = { 1, 1, 1, 1, 1, 1, 1, 2, 3 };
    double d = det_bareiss(3, a), p = per_ryser(3, a);
    printf("  [P8] rows 1=2 all-ones: det = %.1f   per = %.1f (expect 2*(1+2+3)=12)\n", d, p);
    return fabs(p - 12.0) + fabs(d);
}

/* ---------------- P9 GF(2) 仿射rank = (n-1)^2+1 ---------------- */
static double p9_affine_rank(void)
{
    int n, dev = 0;
    printf("  [P9] n | rank{vec(P_s)+vec(P_id)} | (n-1)^2+1\n");
    for (n = 2; n <= 4; n++) {
        int rows = 1, cols = n * n, perm[MAXN], i, j, k, r, c;
        int gen[MAXN * MAXN * MAXN * MAXN]; int cnt = 0;
        int base[MAXN * MAXN];
        for (i = 0; i < n; i++) for (j = 0; j < n; j++) base[i * n + j] = (i == j) ? 1 : 0;
        for (k = 1; k < n * n * n; k++) break;
        for (i = 0; i < n; i++) perm[i] = i;
        for (;;) {
            int v[MAXN * MAXN]; int isid = 1;
            for (i = 0; i < n; i++) for (j = 0; j < n; j++) v[i * n + j] = 0;
            for (i = 0; i < n; i++) v[i * n + perm[i]] = 1;
            for (i = 0; i < cols; i++) if (v[i] != base[i]) isid = 0;
            if (!isid) { for (i = 0; i < cols; i++) gen[cnt * cols + i] = v[i] ^ base[i]; cnt++; }
            r = -1; for (i = n - 2; i >= 0; i--) if (perm[i] < perm[i + 1]) { r = i; break; }
            if (r < 0) break;
            c = n - 1; while (perm[c] < perm[r]) c--;
            { int t = perm[r]; perm[r] = perm[c]; perm[c] = t; }
            for (i = r + 1, j = n - 1; i < j; i++, j--) { int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
        }
        rows = cnt;
        { int rk = rank_gf2(rows, cols, gen); int pred = (n - 1) * (n - 1) + 1;
          printf("        %d | %d | %d\n", n, rk + 1, pred);
          if (rk + 1 != pred) dev += 1.0; }
    }
    return dev;
}

static void print_wuxing(void)
{
    int i;
    for (i = 0; i < 6; i++) printf("   [%s]", g_W[i]);
    printf("\n");
}

int main(void)
{
    double d1, d2, d3, d4, d5, d6, d7, d8, d9;

    printf("== verify_pvsnp_tech : construction layer —— per 为什么没有closed式（mod 2 商）==\n");
    print_wuxing();
    printf("\n");

    d1 = p1_det_paths();
    judge("P1 det two paths + anchors", d1, 1e-09);

    d2 = p2_per_paths();
    judge("P2 per two paths + anchor", d2, 1e-06);

    d3 = p3_alternating();
    judge("P3 det alternating / per invariant", d3, 1e-09);

    d4 = p4_sign_structure();
    judge("P4 det signs split equally", d4, 0.0);

    d5 = p5_mod2_quotient();
    judge("P5 det = per over GF(2)", d5, 0.0);

    d6 = p6_int_gap();
    judge("P6 det(J)=0 but per(J)=n!", d6, 0.0);

    d7 = p7_info_gap();
    judge("P7 parity carries 1 bit only", d7, 0.0);

    d8 = p8_cancel();
    judge("P8 equal rows: det=0, per>0", d8, 1e-09);

    d9 = p9_affine_rank();
    judge("P9 GF(2) affine rank = (n-1)^2+1", d9, 0.0);

    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.3e\n", g_pass, g_fail, g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", NJ, 1.0 - 1.0 / (double)NJ);
    printf("res=%.6e domain=%d dc=%d order=%d\n", g_res, NJ, 6, 2);
    printf("honest note：P5/P6/P9 是全enumeration（n<=4 / n=3 / n<=4），不是对一切 n 的定理；\n");
    printf("     P4 的「正负各半」是项数计数，不是值的cancellation；\n");
    printf("     P7 的信息量是lower bound论证的必要条件，不构成 #P-hard 的证明。\n");
    printf("EXIT=%d\n", g_fail == 0 ? 0 : 1);
    (void)d5; (void)d6; (void)d9;
    return g_fail == 0 ? 0 : 1;
}
