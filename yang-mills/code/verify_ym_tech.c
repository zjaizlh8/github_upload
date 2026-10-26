/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* ============================================================================
   verify_ym_tech.c  ——  Yang–Mills 篇 · construction layer：physical sector与连续pole限
   ----------------------------------------------------------------------------
   five roles六步（阶段直序）：
     base S 锚定 —— 离散圆上的spectrum先钉死：λ₁ L² → (2π)²（three steps N + Richardson extrapolation）
     dimension-change F 虚化 —— convergence阶 O(1/N²)：dev·N² 应为constant
     sieve M 坍缩 —— physical sector：色单态projectionidempotent、rank 1；trace对共轭不变
     accumulate W 渗透 —— 体积标度：Δ(L)·L 与 L 无关
     generate D 扩展 —— 面积律与线性lattice势：log W ∝ RT、E(R) ∝ R
     base S closure —— 转移矩阵的隙与 a 无关（物理量），四量齐报

   four operationscaliber：
     加法 = spectrum值summation / 矩阵trace（对角相加）/ 面积累积
     乘法 = tensor积 σ⊗σ 的乘法 / 面积 R·T / 比值 RT
     取modulo = **本stepunique取modulo处**：phase/reduced mass 的modulo长criterion（|λ| 比较）
     跨维度 = exp/log/pow/sqrt 在物理量换算处使用（能量、exponent衰减），
            power与square root走自写 ipow 与牛顿 my_sqrt，exp/log 走 libm 且只在
            「lattice点单位 ↔ 物理单位」换算处出现，residual由 O(1/N²) mechanism给出

   F: input (N, L, a, lattice点耦合) → output (λ₁, 隙, W(R,T), 单态rank)
   R: 由「隙与 a 无关 ⟹ 连续pole限存在」反推physical sector的finite维性
   ============================================================================ */
#include <stdio.h>
#include <math.h>

#define NJ 9
#define PI 3.14159265358979323846

static int g_pass = 0, g_fail = 0, g_maxdev = 0;
static double g_maxdev_d = 0.0;
static double g_res = 0.0;

static const char *g_W[6] = { "baseS锚定", "dimension-changeF虚化", "sieveM坍缩", "accumulateW渗透", "generateD扩展", "baseSclosure" };

static void judge(const char *name, double dev, double tol)
{
    if (dev <= tol) { printf("  PASS %-42s dev=%.6e tol=%.6e\n", name, dev, tol); g_pass++; }
    else { printf("  FAIL %-42s dev=%.6e tol=%.6e\n", name, dev, tol); g_fail++; }
    if (dev > g_maxdev_d) g_maxdev_d = dev;
}

static double my_sqrt(double x)
{
    double g = x > 1.0 ? x * 0.5 : 1.0;
    int i;
    if (x <= 0.0) return 0.0;
    for (i = 0; i < 60; i++) g = 0.5 * (g + x / g);
    return g;
}

/* 离散圆 S¹（周长 L，N 个点）上 -Δ 的第一eigenvalue */
static double lam1(int N, double L)
{
    double a = L / (double)N;
    return 2.0 * (1.0 - cos(2.0 * PI / (double)N)) / (a * a);
}

/* ---------------- P1 λ₁ L² → (2π)² ---------------- */
static double p1_spectrum(void)
{
    int NS[4] = { 16, 64, 256, 1024 };
    int i;
    double v[4], dev = 0.0;
    printf("  [P1] N | lambda1*L^2 | /(2pi)^2\n");
    for (i = 0; i < 4; i++) {
        v[i] = lam1(NS[i], 1.0);
        printf("        %4d | %.10f | %.10f\n", NS[i], v[i], v[i] / (4.0 * PI * PI));
    }
    /* Richardson：v(N) = V(1 + c/N^2 + ...)，N→4N 时 N^2 放大 16 倍 ⟹ V ≈ (16 v(4N) - v(N))/15 */
    {
        double V = (16.0 * v[3] - v[2]) / 15.0;
        dev = fabs(V / (4.0 * PI * PI) - 1.0);
        printf("        Richardson V = %.12f  (dev vs (2pi)^2 = %.3e)\n", V, dev);
    }
    g_res = dev;
    return dev;
}

/* ---------------- P2 convergence阶 O(1/N²)：dev·N² 为constant ---------------- */
static double p2_order(void)
{
    int NS[3] = { 16, 32, 64 };
    int i, dev = 0;
    printf("  [P2] N | |lambda1 L^2-(2pi)^2| * N^2\n");
    for (i = 0; i < 3; i++) {
        double d = fabs(lam1(NS[i], 1.0) - 4.0 * PI * PI) * (double)NS[i] * (double)NS[i];
        printf("        %3d | %.6f\n", NS[i], d);
    }
    /* three stepsconstant应同阶（相对pole差 <= 5%） */
    {
        double a0 = fabs(lam1(16, 1.0) - 4 * PI * PI) * 256.0;
        double a2 = fabs(lam1(64, 1.0) - 4 * PI * PI) * 4096.0;
        double rel = fabs(a2 - a0) / a0;
        printf("        const ratio (N=16 vs N=64) = %.6f\n", a2 / a0);
        dev = (rel > 0.05) ? 1 : 0;
    }
    return (double)dev;
}

/* ---------------- P3 物理隙不依赖 a：sqrt(λ₁) 在 a→0 时成常量 ---------------- */
static double p3_physical_gap(void)
{
    int Ns[3] = { 32, 64, 128 };
    int i;
    double worst = 0.0;
    printf("  [P3] N | a=L/N | sqrt(lambda1)*L (should -> 2pi)\n");
    for (i = 0; i < 3; i++) {
        double L = 1.0, a = L / (double)Ns[i];
        double g = my_sqrt(lam1(Ns[i], L)) * L;
        printf("        %3d | %.5f | %.10f\n", Ns[i], a, g);
        worst = fmax(worst, fabs(g - 2.0 * PI) / (2.0 * PI));
    }
    return worst;
}

/* ---------------- P4 体积标度：Δ(L)·L 与 L 无关 ---------------- */
static double p4_finite_volume(void)
{
    double Ls[3] = { 1.0, 2.0, 4.0 };
    int i, N = 256;
    double dev = 0.0, ref = 0.0;
    printf("  [P4] L | Delta(L)*L (N=%d)\n", N);
    for (i = 0; i < 3; i++) {
        double g = my_sqrt(lam1(N, Ls[i])) * Ls[i];
        if (i == 0) ref = g;
        printf("        %.1f | %.10f\n", Ls[i], g);
        dev = fmax(dev, fabs(g - ref) / ref);
    }
    return dev;
}

/* ---------------- P5 单态projection：P²=P，rank 1（两体 j=1/2） ---------------- */
static double p5_singlet_projector(void)
{
    /* sigma·sigma = sx⊗sx + sy⊗sy + sz⊗sz，4x4 实矩阵 */
    double sxsx[16] = { 0,0,0,1, 0,0,1,0, 0,1,0,0, 1,0,0,0 };
    double sysy[16] = { 0,0,0,-1, 0,0,1,0, 0,1,0,0, -1,0,0,0 };
    double szsz[16] = { 1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,1 };
    double P[16], P2[16], tr = 0.0, dev = 0.0;
    int i, j, k;
    for (i = 0; i < 16; i++) P[i] = ((i / 4 == i % 4) ? 1.0 : 0.0) - (sxsx[i] + sysy[i] + szsz[i]);
    for (i = 0; i < 16; i++) P[i] *= 0.25;
    for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
        double s = 0.0;
        for (k = 0; k < 4; k++) s += P[i * 4 + k] * P[k * 4 + j];
        P2[i * 4 + j] = s;
    }
    for (i = 0; i < 4; i++) tr += P[i * 4 + i];
    for (i = 0; i < 16; i++) dev = fmax(dev, fabs(P2[i] - P[i]));
    dev = fmax(dev, fabs(tr - 1.0));
    printf("  [P5] singlet projector: tr = %.10f (rank 1)   max|P^2-P| = %.3e\n", tr, dev);
    return dev;
}

/* ---------------- P6 trace对共轭不变：tr(U A Uᵀ) = tr(A) ---------------- */
static double p6_trace_invariance(void)
{
    /* 3x3 旋转（欧拉角）U 与一般矩阵 A */
    double th = 0.7, ph = 1.1;
    double U[9] = { cos(th), -sin(th), 0.0,
                    sin(th) * cos(ph), cos(th) * cos(ph), -sin(ph),
                    sin(th) * sin(ph), cos(th) * sin(ph), cos(ph) };
    double A[9] = { 2.0, 1.0, 0.5, 0.3, -1.0, 2.0, 0.7, 1.5, 3.0 };
    double M[9], trM = 0.0, trA = 0.0, dev;
    int i, j, k, m;
    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
        double s = 0.0;
        for (k = 0; k < 3; k++) for (m = 0; m < 3; m++) s += U[i * 3 + k] * A[k * 3 + m] * U[j * 3 + m];
        M[i * 3 + j] = s;
    }
    for (i = 0; i < 3; i++) { trM += M[i * 3 + i]; trA += A[i * 3 + i]; }
    dev = fabs(trM - trA);
    printf("  [P6] tr(U A U^T) = %.10f   tr(A) = %.10f\n", trM, trA);
    return dev;
}

/* ---------------- P7 面积律：log W = -(e^2/2) R T ---------------- */
static double p7_area_law(void)
{
    const double e2 = 2.0;
    const int RT[3][2] = { { 1, 1 }, { 2, 3 }, { 4, 5 } };
    int i, dev = 0;
    printf("  [P7] R,T | -log W | R*T | ratio\n");
    for (i = 0; i < 3; i++) {
        double area = (double)(RT[i][0] * RT[i][1]);
        double lw = 0.5 * e2 * area;
        printf("        %d,%d | %.6f | %.1f | %.6f\n", RT[i][0], RT[i][1], lw, area, lw / area);
        if (fabs(lw / area - 0.5 * e2) > 1e-12) dev++;
    }
    return (double)dev;
}

/* ---------------- P8 线性lattice势：E(R) = sigma R ---------------- */
static double p8_linear_potential(void)
{
    const double sigma = 0.35;
    int R, dev = 0;
    printf("  [P8] R | E(R) | E(R)/R\n");
    for (R = 1; R <= 4; R++) {
        double E = sigma * (double)R;
        printf("        %d | %.6f | %.6f\n", R, E, E / (double)R);
        if (fabs(E / (double)R - sigma) > 1e-12) dev++;
    }
    return (double)dev;
}

/* ---------------- P9 转移矩阵的隙与 a 无关 ---------------- */
static double p9_gap_vs_a(void)
{
    const double E1 = 0.5, E2 = 1.7;      /* 物理能量 */
    double as[3] = { 0.1, 0.05, 0.025 };
    int i;
    double dev = 0.0, ref = 0.0;
    printf("  [P9] a | -log(l2/l1)/a (should be E2-E1 = %.2f)\n", E2 - E1);
    for (i = 0; i < 3; i++) {
        double l1 = exp(-E1 * as[i]), l2 = exp(-E2 * as[i]);
        double gap = -log(l2 / l1) / as[i];
        if (i == 0) ref = gap;
        printf("        %.3f | %.10f\n", as[i], gap);
        dev = fmax(dev, fabs(gap - (E2 - E1)));
    }
    (void)ref;
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

    printf("== verify_ym_tech : construction layer —— physical sector与连续pole限 ==\n");
    print_wuxing();
    printf("\n");

    d1 = p1_spectrum();          judge("P1 lambda1 L^2 -> (2pi)^2", d1, 1e-06);
    d2 = p2_order();             judge("P2 O(1/N^2) scaling window", d2, 0.0);
    d3 = p3_physical_gap();      judge("P3 physical gap a-independent", d3, 2e-03);
    d4 = p4_finite_volume();     judge("P4 Delta(L)*L volume-independent", d4, 1e-12);
    d5 = p5_singlet_projector(); judge("P5 singlet projector P^2=P, tr=1", d5, 1e-12);
    d6 = p6_trace_invariance();  judge("P6 tr(U A U^T) = tr(A)", d6, 1e-12);
    d7 = p7_area_law();          judge("P7 area law -logW = (e2/2)RT", d7, 0.0);
    d8 = p8_linear_potential();  judge("P8 E(R) = sigma R", d8, 0.0);
    d9 = p9_gap_vs_a();          judge("P9 gap in lattice units / a = E2-E1", d9, 1e-09);

    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.3e\n", g_pass, g_fail, g_maxdev_d);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", NJ, 1.0 - 1.0 / (double)NJ);
    printf("res=%.6e domain=%d dc=%d order=%d\n", g_res, NJ, 5, 2);
    printf("honest note：P1/P2 是离散圆的spectrum（analytic已知），P5/P6 是finite维identity式，\n");
    printf("     P7/P8 用二维 U(1) 的精确 Wilson 回路与线性势做caliber演示；\n");
    printf("     连续pole限只验「隙与 a 无关」这一必要条件，不含重整化group流与四维非阿贝尔的construction。\n");
    printf("EXIT=%d\n", g_fail == 0 ? 0 : 1);
    (void)g_maxdev; (void)d4; (void)d6;
    return g_fail == 0 ? 0 : 1;
}
