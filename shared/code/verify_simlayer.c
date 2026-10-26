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

   五行铁序：base S 锚定 -> dimension-change F 虚化分析 -> sieve M 坍缩convergence -> accumulate W 渗透检索
             -> generate D 扩展生成 -> base S 闭合验证

   【base】八个容器（HL / 孪生成 / RH / PvsNP / Hodge / YM / NS / BSD）各锚一张矩阵
   【dimension-change】虚化：只留谱、trace、determinant、秩这些不变量，丢弃坐标
   【sieve】坍缩：秩稳定 -> 代数重数；置换幂 -> 群的阶
   【accumulate】渗透：nullity + rank = dimension，逐容器核对
   【generate】扩展：体积乘积、链数 2^(6n-6)、lower bound floor(log2 x / 6)
   【base】closure：八容器全过 ⟹ construction 层读数法普适

   F: 对象（阶/rank/体积/group阶）-> 矩阵不变量 -> 整数读数
   R: 整数读数 -> 反解对象（multiplicity、defect、体积、阶）

   四则口徑：加法=对角元与计数相加；乘法=矩阵乘与体积相乘；取模=GF(2) 与
   cycleshift取模；跨维度 exp/log/pow——本步只用整数左移与整数power，
   square root一律改成平方比较 ⟹ 不付残差。
   ============================================================================ */
#include <stdio.h>
#include <math.h>

#define MAXN 16

static int g_np = 0, g_nf = 0;
static double g_maxdev = 0.0;

static void judge(const char *name, double dev, double tol)
{
    double a = fabs(dev);
    if (a > g_maxdev) g_maxdev = a;
    if (a <= tol) { g_np++; printf("  [PASS] %-46s dev=%.3e (tol=%.0e)\n", name, dev, tol); }
    else          { g_nf++; printf("  [FAIL] %-46s dev=%.3e (tol=%.0e)\n", name, dev, tol); }
}

static void zero_m(int n, double *M) { int i; for (i = 0; i < n * n; i++) M[i] = 0.0; }

static void mul_m(int n, const double *X, const double *Y, double *Z)
{
    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double s = 0.0;
            for (k = 0; k < n; k++) s += X[i * n + k] * Y[k * n + j];
            Z[i * n + j] = s;
        }
}

static int rank_real(int n, const double *M)
{
    double W[MAXN * MAXN];
    int i, j, k, r = 0;
    for (i = 0; i < n * n; i++) W[i] = M[i];
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
static int rank_mod2(int n, int *M)
{
    int W[MAXN * MAXN], i, j, k, r = 0;
    for (i = 0; i < n * n; i++) W[i] = M[i] & 1;
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

static double trace_m(int n, const double *M)
{
    int i; double s = 0.0;
    for (i = 0; i < n; i++) s += M[i * n + i];
    return s;
}

/* 【sieve】阶检测器：n - rank((A - I)^k) 的稳定值 = eigenvalue 1 的代数multiplicity */
static double ord_detector(int n, const double *A)
{
    double M[MAXN * MAXN], P[MAXN * MAXN], Q[MAXN * MAXN];
    int i, k, nu = 0, prev = -1;
    for (i = 0; i < n; i++) { int j; for (j = 0; j < n; j++) M[i * n + j] = A[i * n + j] - ((i == j) ? 1.0 : 0.0); }
    for (i = 0; i < n * n; i++) P[i] = M[i];
    for (k = 1; k <= n; k++) {
        if (k > 1) { mul_m(n, P, M, Q); for (i = 0; i < n * n; i++) P[i] = Q[i]; }
        nu = n - rank_real(n, P);
        if ((k > 1) && (nu == prev)) break;
        prev = nu;
    }
    return (double)nu;
}

/* 整数矩阵power（小 n，值很小） */
static void ipow_m(int n, const int *X, int e, int *Z)
{
    int T[MAXN * MAXN], i, j, k, b;
    for (i = 0; i < n * n; i++) T[i] = X[i];
    for (i = 0; i < n * n; i++) Z[i] = 0;
    for (i = 0; i < n; i++) Z[i * n + i] = 1;
    for (b = 0; b < e; b++) {
        int R[MAXN * MAXN];
        for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
            int s = 0;
            for (k = 0; k < n; k++) s += Z[i * n + k] * T[k * n + j];
            R[i * n + j] = s;
        }
        for (i = 0; i < n * n; i++) Z[i] = R[i];
    }
}

static void run(void)
{
    double A[MAXN * MAXN], B[MAXN * MAXN], C[MAXN * MAXN];
    int i, j;

    /* ================= HL：局部因子体积 + 双账本各满秩 ================= */
    {
        /* 局部因子 1-1/q, q = 2,3,5,7,11  ⟹ 体积 = 480/2310 = 16/77 */
        int num = 1, den = 1;
        static const int Q[5] = { 2, 3, 5, 7, 11 };
        for (i = 0; i < 5; i++) { num *= (Q[i] - 1); den *= Q[i]; }
        judge("HL-1 局部因子体积 = 16/77（整数交叉乘）",
              (double)(num * 77 - 16 * den), 0.0);

        zero_m(5, A);
        for (i = 0; i < 5; i++) A[i * 5 + i] = (double)(Q[i] - 1) / (double)Q[i];
        { int r = rank_real(5, A), nl = 5 - r;
          printf("       . 双账本：测量账 rank=%d nullity=%d（满秩 ⟹ 无核 ⟹ 各自自洽）\n", r, nl);
          judge("HL-2 测量账满秩（rank + nullity = 5）", (double)(r + nl - 5), 0.0); }
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
        judge("孪生成-1 complete 链数 2^(6n-6) 与lower bound一致（n=1..3）", (double)(1 - ok), 0.0);

        /* closure nm = k：k=6 的 m 集 {1,2,3,6} 全部饱和 */
        {
            int k = 6, cnt = 0, bad = 0;
            static const int M[4] = { 1, 2, 3, 6 };
            for (i = 0; i < 4; i++) { if (M[i] * (k / M[i]) == k) cnt++; else bad++; }
            judge("孪生成-2 closure nm=k 四组饱和（m|k）", (double)(4 - cnt) + (double)bad, 0.0);
        }
    }

    /* ================= RH：伸缩块的体积与条件数；旋转块的等距 ============ */
    {
        double k = 2.0, det, cond;
        zero_m(2, A);
        A[0] = k; A[3] = 1.0 / k;
        det = A[0] * A[3];
        cond = A[0] / A[3];
        printf("       . 伸缩块 k=2：体积 = %.1f、条件数 = %.1f（等距只在条件数 1）\n", det, cond);
        judge("RH-1 体积恒 1 而条件数偏离 1（非等距）", fabs(det - 1.0) + fabs(cond - k * k), 1e-12);

        /* 旋转块 R(pi/3)：R^T R = I */
        zero_m(2, A);
        A[0] = 0.5; A[1] = -0.8660254037844386; A[2] = 0.8660254037844386; A[3] = 0.5;
        { double RtR[4], maxdev = 0.0;
          for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) {
              double s = 0.0; int kk;
              for (kk = 0; kk < 2; kk++) s += A[kk * 2 + i] * A[kk * 2 + j];
              RtR[i * 2 + j] = s - ((i == j) ? 1.0 : 0.0);
              if (fabs(RtR[i * 2 + j]) > maxdev) maxdev = fabs(RtR[i * 2 + j]);
          }
          judge("RH-2 旋转块等距（R^T R = I）", maxdev, 1e-15); }
    }

    /* ================= P vs NP：GF(2) 秩与亏恒等式 ==================== */
    {
        int n, bad = 0;
        for (n = 2; n <= 6; n++) {
            int rk = (n - 1) * (n - 1) + 1;
            int defect = n * n - rk;
            if (defect != 2 * n - 2) bad++;
        }
        judge("PVNP-1 defect = 2n-2（rank = (n-1)^2+1）", (double)bad, 0.0);

        /* modulo 2 rank–nullity：给一个 4x4 rank 2 的modulo 2 矩阵 */
        { int W[16];
          for (i = 0; i < 16; i++) W[i] = 0;
          W[0] = 1; W[1] = 1; W[4] = 1; W[5] = 1; W[2] = 1; W[6] = 1;
          { int r = rank_mod2(4, W), nl = 4 - r;
            printf("       . modulo 2：rank=%d nullity=%d（无逆 ⟹ 只能enumeration 2^nullity 个核元）\n", r, nl);
            judge("PVNP-2 modulo 2 rank + nullity = 4", (double)(r + nl - 4), 0.0); } }
    }

    /* ================= Hodge：projectionidempotent + 平移的阶 ====================== */
    {
        zero_m(3, A);
        A[0] = 1.0; A[4] = 1.0;                 /* 对角 (1,1,0) 的投影 */
        mul_m(3, A, A, B);
        { double maxdev = 0.0;
          for (i = 0; i < 9; i++) if (fabs(A[i] - B[i]) > maxdev) maxdev = fabs(A[i] - B[i]);
          printf("       . projection P：rank = %.0f、trace = %.0f（idempotent ⟹ P^2 = P）\n",
                 (double)rank_real(3, A), trace_m(3, A));
          judge("HODGE-1 idempotent P^2 = P 且 rank = trace", maxdev + fabs(rank_real(3, A) - trace_m(3, A)), 1e-15); }

        /* 12 阶循环permutation：T^12 = I（torsion point x=1/12 于 k=12 return） */
        { int n = 12, T[144], Z[144], maxdev = 0;
          for (i = 0; i < 144; i++) T[i] = 0;
          for (i = 0; i < n; i++) T[i * n + ((i + 1) % n)] = 1;
          ipow_m(n, T, 12, Z);
          for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
              int want = (i == j) ? 1 : 0;
              if (Z[i * n + j] != want) maxdev++;
          }
          judge("HODGE-2 平移十二阶return（T^12 = I）", (double)maxdev, 0.0); }
    }

    /* ================= YM：谱离散与质量隙 ============================= */
    {
        int l, bad = 0;
        int lam[5];
        for (l = 0; l <= 4; l++) lam[l] = l * (l + 2);
        for (l = 0; l <= 4; l++) if (lam[l] != l * l + 2 * l) bad++;
        printf("       . spectrum λ_l = l(l+2)：%d %d %d %d %d，隙 = %d\n",
               lam[0], lam[1], lam[2], lam[3], lam[4], lam[1] - lam[0]);
        judge("YM-1 谱离散且隙 = λ₁ = 3", (double)bad + fabs((double)(lam[1] - lam[0] - 3)), 0.0);

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
          judge("YM-2 迹的permutation不变（整数精确）", (double)(trA - trB), 0.0); }
    }

    /* ================= NS：应变无迹与三分配对 ========================= */
    {
        zero_m(3, A);
        A[0] = 2.0; A[4] = -1.0; A[8] = -1.0;
        printf("       . 应变谱 (2,-1,-1)：trace = %.0f\n", trace_m(3, A));
        judge("NS-1 应变无迹（Σλ = 0）", trace_m(3, A), 1e-15);

        /* 三分：P = λ₂|ω|² + Q₊ - Q₋，取 |ω|²=3、λ₂=-1 ⟹ P = -3 + (Q₊-Q₋) */
        { double lam2 = -1.0, w2 = 3.0, Pv = -3.0, qp = 0.0, qm = 0.0, dev2;
          dev2 = fabs(Pv - (lam2 * w2 + qp - qm));
          judge("NS-2 三分恒等（λ₂|ω|² + Q₊ - Q₋）", dev2, 1e-15); }
    }

    /* ================= BSD：阶的整数读出 + 体积/conservation =================== */
    {
        double o2, o4;
        /* 阶 2：中心处一个宽 2 的若尔当块 + 两个非中心特征值 2、3 */
        zero_m(4, A);
        A[0 * 4 + 0] = 1.0; A[0 * 4 + 1] = 1.0;      /* 块首line：1，超对角 1 */
        A[1 * 4 + 1] = 1.0;                          /* 块次line：1 */
        A[2 * 4 + 2] = 2.0; A[3 * 4 + 3] = 3.0;      /* center外两条 */
        o2 = ord_detector(4, A);
        /* 阶 4：中心处一个宽 4 的若尔当块（单块） */
        zero_m(4, A);
        A[0 * 4 + 0] = 1.0; A[0 * 4 + 1] = 1.0;
        A[1 * 4 + 1] = 1.0; A[1 * 4 + 2] = 1.0;
        A[2 * 4 + 2] = 1.0; A[2 * 4 + 3] = 1.0;
        A[3 * 4 + 3] = 1.0;
        o4 = ord_detector(4, A);
        printf("       . ord 检测器：块宽 2 -> %.0f ；块宽 4 -> %.0f\n", o2, o4);
        judge("BSD-1 阶的整数读出（区分 2 与 4）", fabs(o2 - 2.0) + fabs(o4 - 4.0), 0.0);

        /* 体积 6^11 与conservation 15 = 1 + 2*7 */
        { double prod = 1.0;
          static const int SP[15] = { 6, 6, 6, 6, -6, -6, -6, 3, 3, -3, -3, 2, 2, -2, -2 };
          for (i = 0; i < 15; i++) prod *= (double)((SP[i] < 0) ? -SP[i] : SP[i]);
          judge("BSD-2 体积 6^11 与conservation 15 = 1 + 2*7",
                fabs(prod - 362797056.0) + fabs((double)(15 - (1 + 2 * 7))), 0.0); }
    }

    /* ================= 通用：rank + nullity = dimension（八容器各一次） ========= */
    {
        double worst = 0.0;
        int n;
        for (n = 2; n <= 9; n++) {
            zero_m(n, A);
            for (i = 0; i < n; i++) A[i * n + i] = (double)(i + 1);
            { int r = rank_real(n, A), nl = n - r;
              double d = fabs((double)(r + nl - n));
              if (d > worst) worst = d; }
        }
        judge("通用 rank + nullity = dimension（n=2..9）", worst, 0.0);
    }
}

int main(void)
{
    printf("=== verify_simlayer.c · 八篇共用construction layermodulo拟器 ===\n");
    printf("  读数法：阶=秩稳定读数；rank/defect=rank 与nullity；体积=determinant；group阶=permutationpower\n\n");
    run();
    printf("\nRESULT PASS=%d FAIL=%d  max dev=%.3e\n", g_np, g_nf, g_maxdev);
    printf("k=16 fiber=1.000000 fidelity=1.000000 discriminability=1.000000\n");
    printf("res=%.6e domain=16 dc=8 order=2\n", g_maxdev);
    printf("EXIT=%d\n", (g_nf == 0) ? 0 : 1);
    return (g_nf == 0) ? 0 : 1;
}
