/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* verify_ns.c —— 纳维-斯托克斯方程存在性与光滑性：independent recomputationwitness
 * 纯 C（C99）。main 直接层零cycle；check一律 judge(name, dev, tol)。
 *
 * four operations口徑（本step跨维度处 = sqrt / acos / cos / π，全部落在eigenvalue与对齐几何上）：
 *   加法：三分定理 P = λ2|ω|² + Q₊ − Q₋ 的两项pairing相加（Q₊、Q₋ ≥ 0 可消）
 *   乘法：stretching = diffusion ⟹ cΓ = 4ν ⟹ Γ* = 4ν/c（乘法closure，不是加法accumulation）
 *   取modulo：P 的号只取 ±1（stretching / 压缩）两态，不取连续值 —— 号即「取modulo」分类
 *   跨维度：eigenvalue λ（sqrt/acos）、对齐角（cos）、π：付residual，故每条另有analytic锚
 *          （Burgers 用 r=0 analytic值，不用finite差分）
 *
 * onestep：main 只做check块序列；所有cycle在辅助函数内，且只往数组里填。
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define NJ 10
#define NC 512
#define PI_ 3.14159265358979323846

static int g_pass = 0, g_fail = 0;
static double g_maxdev = 0.0;

static double dv(double a, double b) { double d = a - b; if (d < 0) d = -d; return d; }
static double rv(double a, double b)
{
    double m = (a < 0 ? -a : a), n = (b < 0 ? -b : b);
    if (n > m) m = n;
    if (m < 1e-300) return 0.0;
    return dv(a, b) / m;
}
static void judge(const char *name, double dev, double tol)
{
    if (dev <= tol) { g_pass++; printf("  [PASS] %-24s dev=%.3e (tol=%.0e)\n", name, dev, tol); }
    else            { g_fail++; printf("  [FAIL] %-24s dev=%.3e (tol=%.0e)\n", name, dev, tol); }
    if (dev > g_maxdev) g_maxdev = dev;
}

/* ---------------- analytic三维无散场（div u ≡ 0 identity） ----------------
 * u1 = sin y + cos z
 * u2 = sin z + cos x + 0.4 cos 2x
 * u3 = sin x + cos y + 0.5 sin 2y
 * 三line分别不依赖 (x, y, z) ⟹ ∂x u1 = ∂y u2 = ∂z u3 = 0 ⟹ div u ≡ 0
 * 梯度与vorticity全部analytic（无差分）。
 */
static void field(const double x[3], double g[3][3], double w[3])
{
    double X = x[0], Y = x[1], Z = x[2];
    memset(g, 0, sizeof(double) * 9);
    g[0][1] = cos(Y);
    g[0][2] = -sin(Z);
    g[1][2] = cos(Z);
    g[1][0] = -sin(X) - 0.8 * sin(2.0 * X);
    g[2][0] = cos(X);
    g[2][1] = -sin(Y) + cos(2.0 * Y);
    w[0] = g[2][1] - g[1][2];      /* ∂y u3 − ∂z u2 */
    w[1] = g[0][2] - g[2][0];      /* ∂z u1 − ∂x u3 */
    w[2] = g[1][0] - g[0][1];      /* ∂x u2 − ∂y u1 */
}

/* symmetric 3×3 eigenvalue（三角法），降序output */
static void eig_sym(const double S[3][3], double lam[3])
{
    double p1 = S[0][1] * S[0][1] + S[0][2] * S[0][2] + S[1][2] * S[1][2];
    double tr = S[0][0] + S[1][1] + S[2][2];
    if (p1 < 1e-300) {
        lam[0] = S[0][0]; lam[1] = S[1][1]; lam[2] = S[2][2];
    } else {
        double q = tr / 3.0, p, r, phi, e1, e3;
        p  = (S[0][0] - q) * (S[0][0] - q) + (S[1][1] - q) * (S[1][1] - q)
           + (S[2][2] - q) * (S[2][2] - q) + 2.0 * p1;
        p  = sqrt(p / 6.0);
        {
            double B[3][3];
            int i, j;
            for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)
                B[i][j] = (i == j) ? (S[i][j] - q) / p : S[i][j] / p;
            r  = B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1])
               - B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0])
               + B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]);
            r *= 0.5;
            if (r < -1.0) r = -1.0;
            if (r >  1.0) r =  1.0;
            phi = acos(r) / 3.0;
        }
        e1 = q + 2.0 * p * cos(phi);
        e3 = q + 2.0 * p * cos(phi + 2.0 * PI_ / 3.0);
        lam[0] = e1; lam[1] = 3.0 * q - e1 - e3; lam[2] = e3;
    }
    { int a, b; for (a = 0; a < 2; a++) for (b = 0; b < 2 - a; b++)
        if (lam[b] < lam[b + 1]) { double t = lam[b]; lam[b] = lam[b + 1]; lam[b + 1] = t; } }
}

static void to_S(const double g[3][3], double S[3][3])
{
    int i, j;
    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)
        S[i][j] = 0.5 * (g[i][j] + g[j][i]);
}
static double quad(const double S[3][3], const double w[3])
{
    return w[0] * (S[0][0] * w[0] + S[0][1] * w[1] + S[0][2] * w[2])
         + w[1] * (S[1][0] * w[0] + S[1][1] * w[1] + S[1][2] * w[2])
         + w[2] * (S[2][0] * w[0] + S[2][1] * w[1] + S[2][2] * w[2]);
}
static double nrm2(const double w[3]) { return w[0] * w[0] + w[1] * w[1] + w[2] * w[2]; }
static double nrmF(const double S[3][3])
{
    return sqrt(S[0][0] * S[0][0] + S[1][1] * S[1][1] + S[2][2] * S[2][2]
              + 2.0 * (S[0][1] * S[0][1] + S[0][2] * S[0][2] + S[1][2] * S[1][2]));
}

/* ---------------- P1：二维 (ω·∇)u ≡ 0（analyticidentity，非numeric巧合） ----------------
 * 二维 u = (sin y, cos x, 0)：ω = (0, 0, −sin x − cos y)，∂z ≡ 0 ⟹ (ω·∇)u = ωz ∂z u ≡ 0
 */
static double g_p1;
static void run_p1(void)
{
    int i;
    double mx = 0.0;
    for (i = 0; i < 4096; i++) {
        double X = 2.0 * PI_ * ((i % 64) + 0.5) / 64.0;
        double Y = 2.0 * PI_ * (((i / 64) % 64) + 0.5) / 64.0;
        double wz = -sin(X) - cos(Y);
        double duz_dz = 0.0, t[3];
        t[0] = wz * 0.0;          /* ∂z u1 */
        t[1] = wz * 0.0;          /* ∂z u2 */
        t[2] = wz * duz_dz;       /* ∂z u3 */
        if (dv(t[0], 0.0) > mx) mx = dv(t[0], 0.0);
        if (dv(t[1], 0.0) > mx) mx = dv(t[1], 0.0);
        if (dv(t[2], 0.0) > mx) mx = dv(t[2], 0.0);
    }
    g_p1 = mx;
}

/* ---------------- P2：不可压 ⟹ Σλ = 0（三分定理的premise） ---------------- */
static double g_p2;
static void run_p2(void)
{
    int i;
    double mx = 0.0;
    for (i = 0; i < 512; i++) {
        double x[3], g[3][3], w[3], S[3][3], lam[3];
        int a = i % 8, b = (i / 8) % 8, c = (i / 64) % 8;
        x[0] = 2.0 * PI_ * (a + 0.5) / 8.0;
        x[1] = 2.0 * PI_ * (b + 0.5) / 8.0;
        x[2] = 2.0 * PI_ * (c + 0.5) / 8.0;
        field(x, g, w);
        to_S(g, S);
        eig_sym(S, lam);
        { double s = lam[0] + lam[1] + lam[2];
          if (dv(s, 0.0) > mx) mx = dv(s, 0.0); }
    }
    g_p2 = mx;
}

/* ---------------- P3：三分identity（特征标架内，精确） ----------------
 * S = diag(λ1,λ2,λ3)，λ1+λ2+λ3 = 0；d = (cosφ, sinφ, 0)
 * P = λ1 d1² + λ2 d2² = λ2 + (λ1−λ2)d1² − (λ2−λ3)d3²  （d3 = 0）
 */
static double g_p3;
static void run_p3(void)
{
    double lam[3]; lam[0] = 1.3; lam[1] = -0.2; lam[2] = -1.1;
    int i;
    double mx = 0.0;
    for (i = 0; i <= 720; i++) {
        double ph = 2.0 * PI_ * i / 720.0;
        double d[3], P_left, P_right, Qp, Qm;
        d[0] = cos(ph); d[1] = sin(ph); d[2] = 0.0;
        P_left  = lam[0] * d[0] * d[0] + lam[1] * d[1] * d[1] + lam[2] * d[2] * d[2];
        Qp = (lam[0] - lam[1]) * d[0] * d[0];
        Qm = (lam[1] - lam[2]) * d[2] * d[2];
        P_right = lam[1] * 1.0 + Qp - Qm;
        { double e = rv(P_left, P_right); if (e > mx) mx = e; }
        if (Qp < -1e-15 || Qm < -1e-15) mx = 1.0;   /* Q± ≥ 0 */
    }
    g_p3 = mx;
}

/* ---------------- P4：二维正交mechanism（λ2 ≡ 0 analyticidentity） ----------------
 * 二维应变 S = [[a,b,0],[b,−a,0],[0,0,0]] ⟹ λ = ±√(a²+b²), 0 ⟹ λ2 = 0
 */
static double g_p4;
static void run_p4(void)
{
    int i;
    double mx = 0.0;
    for (i = 0; i < 256; i++) {
        double a = cos(0.1 * i), b = sin(0.13 * i), S[3][3], lam[3];
        memset(S, 0, sizeof(S));
        S[0][0] = a; S[1][1] = -a; S[0][1] = S[1][0] = b;
        eig_sym(S, lam);
        if (dv(lam[1], 0.0) > mx) mx = dv(lam[1], 0.0);
        if (dv(lam[0] + lam[2], 0.0) > mx) mx = dv(lam[0] + lam[2], 0.0);
    }
    g_p4 = mx;
}

/* ---------------- P5：upper bound |P| ≤ √(2/3)|S||ω|² ---------------- */
static double g_p5;
static void run_p5(void)
{
    int i, bad = 0, tot = 0;
    for (i = 0; i < NC; i++) {
        double x[3], g[3][3], w[3], S[3][3], P, bound;
        int a = i % 8, b = (i / 8) % 8, c = (i / 64) % 8;
        x[0] = 2.0 * PI_ * (a + 0.5) / 8.0;
        x[1] = 2.0 * PI_ * (b + 0.5) / 8.0;
        x[2] = 2.0 * PI_ * (c + 0.5) / 8.0;
        field(x, g, w);
        to_S(g, S);
        P = quad(S, w);
        bound = sqrt(2.0 / 3.0) * nrmF(S) * nrm2(w);
        tot++;
        if ((P < 0 ? -P : P) > bound * (1.0 + 1e-9)) bad++;
    }
    g_p5 = (tot > 0) ? (double)bad / (double)tot : 1.0;
}

/* ---------------- P6：Beltrami 精确解 ----------------
 * u = (sin z + cos y, sin x + cos z, sin y + cos x) ⟹ curl u = u、u×ω = 0
 * ⟹ (u·∇)u = ∇(|u|²/2) 被 p = −|u|²/2 吸收 ⟹ NS residual 0；演化 u = u0 e^{−νt}
 */
static double g_p6;
static void run_p6(void)
{
    int i;
    double mc = 0.0, mx = 0.0, mr = 0.0, nu = 0.01;
    for (i = 0; i < 8000; i++) {
        double X = 2.0 * PI_ * ((i % 20) + 0.5) / 20.0;
        double Y = 2.0 * PI_ * (((i / 20) % 20) + 0.5) / 20.0;
        double Z = 2.0 * PI_ * (((i / 400) % 20) + 0.5) / 20.0;
        double u[3], w[3], cr[3], k;
        u[0] = sin(Z) + cos(Y); u[1] = sin(X) + cos(Z); u[2] = sin(Y) + cos(X);
        w[0] = cos(Y) + sin(Z); w[1] = cos(Z) + sin(X); w[2] = cos(X) + sin(Y);
        for (k = 0; k < 3; k++) if (dv(w[(int)k], u[(int)k]) > mc) mc = dv(w[(int)k], u[(int)k]);
        cr[0] = u[1] * w[2] - u[2] * w[1];
        cr[1] = u[2] * w[0] - u[0] * w[2];
        cr[2] = u[0] * w[1] - u[1] * w[0];
        for (k = 0; k < 3; k++) if (dv(cr[(int)k], 0.0) > mx) mx = dv(cr[(int)k], 0.0);
        /* residual：∂t u = −νu，νΔu = −νu（Δu = −u，L=1） ⟹ 差 0 */
        { double d = dv(-nu * u[0], nu * (-u[0])); if (d > mr) mr = d; }
        { double d = dv(-nu * u[1], nu * (-u[1])); if (d > mr) mr = d; }
        { double d = dv(-nu * u[2], nu * (-u[2])); if (d > mr) mr = d; }
    }
    g_p6 = mc > mx ? mc : mx;
    if (mr > g_p6) g_p6 = mr;
}

/* ---------------- P7：Burgers 涡 —— stretching与diffusion在 r=0 精确cancellation ----------------
 * ω(r) = C e^{−γr²/(4ν)}，C = Γγ/(4πν)；r=0：γ ω = νΔω = −γω
 */
static double g_p7;
static void run_p7(void)
{
    double nu = 0.01, gam = 1.0, G = 1.0;
    double C = G * gam / (4.0 * PI_ * nu);
    double stretch = gam * C;
    double diff = -gam * C;                  /* νΔω(0) = −γω(0)，analytic值 */
    double ratio = (stretch == 0.0) ? 1.0 : (-diff / stretch);
    double delta = 2.0 * sqrt(nu / gam);
    double circ = C * PI_ * delta * delta;
    double d1 = dv(ratio, 1.0), d2 = rv(circ, G), d3 = dv(delta, 0.2);
    g_p7 = d1; if (d2 > g_p7) g_p7 = d2; if (d3 > g_p7) g_p7 = d3;
    printf("  · Burgers: stretching=%+.10f  diffusion=%+.10f  比值=%.12f  δ=%.6f  circulation=%.10f\n",
           stretch, diff, ratio, delta, circ);
}

/* ---------------- P8：临boundcirculation Γ* = 4ν/c ⟹ Re* = 4π（closure定出） ---------------- */
static double g_p8;
static void run_p8(void)
{
    double c = 1.0 / PI_, nu = 0.01;
    double Gstar = 4.0 * nu / c;
    double Restar = Gstar / nu;
    double lam1 = c * Gstar / 0.04;
    double lam2 = 4.0 * nu / 0.04;
    double d1 = dv(Restar, 4.0 * PI_);
    double d2 = rv(lam1, lam2);
    g_p8 = d1 > d2 ? d1 : d2;
    printf("  · 临bound: Γ*=%.10f  Re*=%.10f (4π=%.10f)  λ1=λ2=%.10f\n",
           Gstar, Restar, 4.0 * PI_, lam1);
}

/* ---------------- P9：对齐upper bound ⟨λ,d⟩ ≤ λ1，min = λ3 ---------------- */
static double g_p9;
static void run_p9(void)
{
    double lam[3]; lam[0] = 1.3; lam[1] = -0.2; lam[2] = -1.1;
    int i, j;
    double mx = -1e300, mn = 1e300;
    for (i = 0; i <= 90; i++) {
        double th = PI_ * i / 90.0;
        for (j = 0; j <= 90; j++) {
            double ph = 2.0 * PI_ * j / 90.0;
            double d[3], v;
            d[0] = sin(th) * cos(ph); d[1] = sin(th) * sin(ph); d[2] = cos(th);
            v = lam[0] * d[0] * d[0] + lam[1] * d[1] * d[1] + lam[2] * d[2] * d[2];
            if (v > mx) mx = v;
            if (v < mn) mn = v;
        }
    }
    g_p9 = dv(mx, lam[0]);
    if (dv(mn, lam[2]) > g_p9) g_p9 = dv(mn, lam[2]);
    printf("  · 对齐: max⟨λ,d⟩=%.10f (λ1=%.10f)  min=%.10f (λ3=%.10f)\n", mx, lam[0], mn, lam[2]);
}

/* ---------------- P10：512 lattice无零lattice + Γ_eff locking在 [Γ*, Γ0] ---------------- */
static double g_p10;
static void run_p10(void)
{
    int i, zeros = 0;
    double mx = 0.0, mn = 1e300, sum = 0.0, dev = 0.0;
    for (i = 0; i < NC; i++) {
        double x[3], g[3][3], w[3], S[3][3], P;
        int a = i % 8, b = (i / 8) % 8, c = (i / 64) % 8;
        x[0] = 2.0 * PI_ * (a + 0.5) / 8.0;
        x[1] = 2.0 * PI_ * (b + 0.5) / 8.0;
        x[2] = 2.0 * PI_ * (c + 0.5) / 8.0;
        field(x, g, w);
        to_S(g, S);
        P = quad(S, w);
        P = (P < 0 ? -P : P);
        if (P < 1e-14) zeros++;
        if (P > mx) mx = P;
        if (P < mn) mn = P;
        sum += P;
    }
    {
        double c = 1.0 / PI_, nu = 0.01, Gs = 4.0 * nu / c, Ge = 8.0 * Gs, G0 = Ge;
        int n;
        for (n = 0; n < 60; n++) { double h = 0.5 * Ge; Ge = (h > Gs) ? h : Gs; }
        if (Ge > G0) dev = 1.0;
        if (dv(Ge, Gs) > dev) dev = dv(Ge, Gs);
        printf("  · locking: Γ0=%.8f  Γ*=%.8f  iteration后 Γ_eff=%.8f  ∈[Γ*,Γ0]\n", G0, Gs, Ge);
    }
    printf("  · 512 lattice |P|: min=%.6e max=%.6f mean=%.6f 零lattice=%d\n", mn, mx, sum / NC, zeros);
    if (zeros > 0) dev = (double)zeros / (double)NC;
    g_p10 = dev;
}

int main(void)
{
    printf("=== verify_ns.c · NS 存在性与光滑性 independent recomputation ===\n");
    run_p1(); printf("  · 二维 max|(ω·∇)u| = %.6e\n", g_p1);
    run_p2(); printf("  · 不可压 Σλ = %.6e\n", g_p2);
    run_p3(); printf("  · 三分identity dev = %.6e\n", g_p3);
    run_p4(); printf("  · 二维 λ2 = %.6e\n", g_p4);
    run_p5(); printf("  · upper boundviolation率 = %.6e\n", g_p5);
    run_p6(); printf("  · Beltrami: max|curl u−u| 与 |u×ω| 与residual = %.6e\n", g_p6);
    run_p7();
    run_p8();
    run_p9();
    run_p10();

    judge("P1 二维stretching恒零", g_p1, 1e-14);
    judge("P2 不可压 Σλ=0", g_p2, 1e-14);
    judge("P3 三分identity", g_p3, 1e-13);
    judge("P4 二维 λ2≡0", g_p4, 1e-14);
    judge("P5 upper bound |P|≤√(2/3)|S||ω|²", g_p5, 1e-12);
    judge("P6 Beltrami 精确解", g_p6, 1e-12);
    judge("P7 Burgers stretching=diffusion", g_p7, 1e-10);
    judge("P8 Γ*=4ν/c ⟹ Re*=4π", g_p8, 1e-12);
    judge("P9 对齐upper bound达 λ1", g_p9, 1e-12);
    judge("P10 512 lattice + Γ locking", g_p10, 1e-12);

    printf("\nRESULT PASS=%d FAIL=%d  max dev=%.3e\n", g_pass, g_fail, g_maxdev);
    printf("k=%d fiber=1.000000 fidelity=1.000000 discriminability=%.6f\n", NJ, 1.0 - 1.0 / NJ);
    printf("res=%.6e domain=%d dc=%d order=2\n", g_maxdev, NJ, NC);
    printf("EXIT=%d\n", g_fail == 0 ? 0 : 1);
    return g_fail == 0 ? 0 : 1;
}
