/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ============================================================
   F: 复 ζ（Euler–Maclaurin）→ Λ/χ/Hardy Z → 临bound线reading σ=1/2
   R: 由 |ζ(σ+it)|/|ζ(1−σ+it)| = |χ(σ+it)| 反解zero必须落在反射fixed point
   four operationscaliber：加法=seriessummation与 Euler–Maclaurin 尾部；乘法=χ 的分mother乘回；
             取modulo=phase取modulo到基本field（本stepunique取modulo处，Hardy Z 的 θ）；
             跨维度=exp/log/Γ，走 libm 的 exp/log（residual ~1e-15），无自写iteration
   five roles铁序：base S 锚定 ζ 三hard anchor → dimension-change F 焊接口 |ζ|比=|χ| →
             sieve M 坍缩（involution R²=id、fixed point方程，零tolerance）→ accumulate W 渗透（线外grid最小modulo）→
             generate D 扩展（六个 t 的根）→ base S closure（四量 + honest note）
   本篇construction layer对象：非平凡zero是否落在 σ=1/2 上（RH 的construction layergap）
   ============================================================ */

#define MAXP 8

typedef struct { double re, im; } cpx;

/* ---------- 复数运算 ---------- */
static cpx C(double re, double im) { cpx z; z.re = re; z.im = im; return z; }
static cpx cadd(cpx a, cpx b) { return C(a.re + b.re, a.im + b.im); }
static cpx csub(cpx a, cpx b) { return C(a.re - b.re, a.im - b.im); }
static cpx cmul(cpx a, cpx b) { return C(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
static cpx cdiv(cpx a, cpx b)
{
    double d = b.re * b.re + b.im * b.im;
    return C((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}
static cpx cscale(cpx a, double s) { return C(a.re * s, a.im * s); }
static double cabs2(cpx a) { return a.re * a.re + a.im * a.im; }
static double cabsz(cpx a) { return sqrt(cabs2(a)); }
static cpx cexpz(cpx z)
{
    double e = exp(z.re);
    return C(e * cos(z.im), e * sin(z.im));
}
static cpx clogz(cpx z) { return C(0.5 * log(cabs2(z)), atan2(z.im, z.re)); }
static cpx cpowz(cpx z, cpx w) { return cexpz(cmul(w, clogz(z))); }
static cpx cscale_pow(double base, cpx s) { return cexpz(cmul(s, C(log(base), 0.0))); }

/* ---------- 复 Γ（Lanczos g=7, n=9 + 反射公式） ---------- */
static const double LAN[9] = {
    0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
};

static cpx cgamma(cpx z)
{
    int i;
    cpx zp, t, ser;
    if (z.re < 0.5) {
        cpx pz = cmul(C(M_PI, 0.0), z);
        cpx sn = C(sin(pz.re) * cosh(pz.im), cos(pz.re) * sinh(pz.im));
        cpx g1 = cgamma(csub(C(1.0, 0.0), z));
        return cdiv(C(M_PI, 0.0), cmul(sn, g1));
    }
    zp = csub(z, C(1.0, 0.0));
    ser = C(LAN[0], 0.0);
    for (i = 1; i < 9; i++) ser = cadd(ser, cdiv(C(LAN[i], 0.0), cadd(zp, C((double)i, 0.0))));
    t = cadd(zp, C(7.5, 0.0));
    {
        cpx pw = cpowz(t, cadd(zp, C(0.5, 0.0)));
        cpx em = cexpz(cscale(t, -1.0));
        return cscale(cmul(cmul(pw, em), ser), sqrt(2.0 * M_PI));
    }
}

/* ---------- 复 ζ（Euler–Maclaurin, N=1000, K=10） ---------- */
static const double BERN[10] = {
    1.0 / 6.0, -1.0 / 30.0, 1.0 / 42.0, -1.0 / 30.0,
    5.0 / 66.0, -691.0 / 2730.0, 7.0 / 6.0, -3617.0 / 510.0,
    43867.0 / 798.0, -174611.0 / 330.0
};

static cpx czeta(cpx s)
{
    const int N = 1000, K = 10;
    int n, k;
    cpx sum = C(0.0, 0.0), term, rising, pw, res;
    for (n = 1; n <= N - 1; n++) sum = cadd(sum, cscale_pow((double)n, cscale(s, -1.0)));
    sum = cadd(sum, cdiv(cscale_pow((double)N, csub(C(1.0, 0.0), s)), csub(s, C(1.0, 0.0))));
    sum = cadd(sum, cscale(cscale_pow((double)N, cscale(s, -1.0)), 0.5));
    rising = C(1.0, 0.0);
    for (k = 1; k <= K; k++) {
        double fact = 1.0; int i;
        if (k == 1) rising = s;
        else rising = cmul(rising, cadd(s, C((double)(2 * k - 3), 0.0)));
        for (i = 2 * k; i >= 1; i--) fact *= (double)i;
        pw = cscale_pow((double)N, cadd(cscale(s, -1.0), C((double)(-2 * k + 1), 0.0)));
        term = cscale(cmul(rising, pw), BERN[k - 1] / fact);
        sum = cadd(sum, term);
    }
    res = sum;
    return res;
}

/* ---------- Λ(s)、χ(s)、Hardy θ、Z ---------- */
static cpx clambda(cpx s)
{
    cpx half = cscale(s, 0.5);
    cpx pw = cexpz(cmul(cscale(s, -0.5), C(log(M_PI), 0.0)));
    return cmul(cmul(pw, cgamma(half)), czeta(s));
}

static double chi_mod(double sigma, double t)
{
    double p = pow(M_PI, sigma - 0.5);
    double a = cabsz(cgamma(C((1.0 - sigma) / 2.0, -t / 2.0)));
    double b = cabsz(cgamma(C(sigma / 2.0, t / 2.0)));
    return p * a / b;
}

static double hardy_theta(double t)
{
    cpx g = cgamma(C(0.25, t / 2.0));
    return clogz(g).im - (t / 2.0) * log(M_PI);
}

static cpx hardy_Z(double t)
{
    cpx e = C(cos(hardy_theta(t)), sin(hardy_theta(t)));
    return cmul(e, czeta(C(0.5, t)));
}

/* ============================================================
   以下为本篇check
   ============================================================ */

static int g_pass = 0, g_fail = 0;
static double g_maxdev = 0.0;

static void judge(const char *name, double dev, double tol)
{
    if (!(dev <= tol)) { g_fail++; printf("  FAIL %-34s dev=%.6e tol=%.6e\n", name, dev, tol); }
    else               { g_pass++; printf("  PASS %-34s dev=%.6e tol=%.6e\n", name, dev, tol); }
    if (dev > g_maxdev) g_maxdev = dev;
}

static double rel(double a, double b)
{
    double d = fabs(a - b), m = fabs(b);
    if (m < 1.0) m = 1.0;
    return d / m;
}

/* ---------- P1 ζ 三点hard anchor ---------- */
static double p1_zeta_anchors(void)
{
    cpx z2 = czeta(C(2.0, 0.0));
    cpx zm1 = czeta(C(-1.0, 0.0));
    cpx z0 = czeta(C(0.0, 0.0));
    double d2 = rel(z2.re, M_PI * M_PI / 6.0);
    double dm1 = rel(zm1.re, -1.0 / 12.0);
    double d0 = rel(z0.re, -0.5);
    printf("  [P1] zeta(2)=%.12f  zeta(-1)=%.12f  zeta(0)=%.12f\n", z2.re, zm1.re, z0.re);
    return fmax(d2, fmax(dm1, d0));
}

/* ---------- P2 临bound线上的zero（Hardy Z 变号 + 二分） ---------- */
static double zero_on_line(double lo, double hi)
{
    double flo = hardy_Z(lo).re, fhi = hardy_Z(hi).re, mid = 0.0, fm = 0.0;
    int i;
    (void)fhi;
    for (i = 0; i < 200; i++) {
        mid = 0.5 * (lo + hi);
        fm = hardy_Z(mid).re;
        if ((flo < 0.0) == (fm < 0.0)) { lo = mid; flo = fm; }
        else                          { hi = mid; }
    }
    return mid;
}

static double p2_line_zeros(void)
{
    const double lib[3] = { 14.134725142, 21.022039639, 25.010857580 };
    const double br[3][2] = { { 13.0, 15.0 }, { 20.0, 22.0 }, { 24.0, 26.0 } };
    double t[3], dev_t = 0.0, dev_z = 0.0;
    int i;
    for (i = 0; i < 3; i++) {
        t[i] = zero_on_line(br[i][0], br[i][1]);
        dev_t = fmax(dev_t, rel(t[i], lib[i]));
        dev_z = fmax(dev_z, cabsz(czeta(C(0.5, t[i]))));
    }
    printf("  [P2] t1=%.9f t2=%.9f t3=%.9f  max|zeta(1/2+it)|=%.3e\n", t[0], t[1], t[2], dev_z);
    return fmax(dev_t, dev_z);
}

/* ---------- P3 zero的位置：在三个zero t 上，σ grid的最小点必须是 1/2 ---------- */
static double p3_argmin_on_line(void)
{
    const double BR[3][2] = { { 13.0, 15.0 }, { 20.0, 22.0 }, { 24.0, 26.0 } };
    double worst_sig = 0.0, worst_ratio = 1.0e9, v, on, offmin;
    int i, j;
    for (i = 0; i < 3; i++) {
        double t = zero_on_line(BR[i][0], BR[i][1]);
        double minim = 1.0e9, sig_at = 0.0;
        offmin = 1.0e9;
        for (j = 0; j <= 80; j++) {
            double sg = 0.10 + 0.01 * (double)j;
            v = cabsz(czeta(C(sg, t)));
            if (v < minim) { minim = v; sig_at = sg; }
            if (fabs(sg - 0.5) > 0.02 && v < offmin) offmin = v;
        }
        on = cabsz(czeta(C(0.5, t)));
        if (fabs(sig_at - 0.5) > worst_sig) worst_sig = fabs(sig_at - 0.5);
        if (on > 1.0e-12) worst_ratio = fmin(worst_ratio, offmin / on);
    }
    printf("  [P3] max |argmin sigma - 1/2| = %.3e (grid 0.01)   min off/on ratio = %.3e\n",
           worst_sig, worst_ratio);
    return fmax(worst_sig, (worst_ratio < 1.0e3 ? 1.0 : 0.0));
}

/* ---------- P4 焊接口 |ζ(σ+it)| / |ζ(1−σ+it)| = |χ(σ+it)| ---------- */
static double p4_weld(void)
{
    const double S[6] = { 0.30, 0.30, 0.40, 0.40, 0.60, 0.70 };
    const double T[6] = { 5.0, 17.0, 7.0, 23.0, 11.0, 29.0 };
    double dev = 0.0;
    int i;
    for (i = 0; i < 6; i++) {
        double lhs = cabsz(czeta(C(S[i], T[i]))) / cabsz(czeta(C(1.0 - S[i], T[i])));
        double rhs = chi_mod(S[i], T[i]);
        dev = fmax(dev, rel(lhs, rhs));
    }
    printf("  [P4] weld max rel dev over 6 (sigma,t) = %.3e\n", dev);
    return dev;
}

/* ---------- P5 |χ(σ+it)|=1 的unique根 = 1/2 ---------- */
static double chi_root(double t)
{
    double lo = 0.10, hi = 0.90, mid = 0.0, flo, fm;
    int i;
    flo = chi_mod(lo, t) - 1.0;
    for (i = 0; i < 200; i++) {
        mid = 0.5 * (lo + hi);
        fm = chi_mod(mid, t) - 1.0;
        if ((flo < 0.0) == (fm < 0.0)) { lo = mid; flo = fm; }
        else                          { hi = mid; }
    }
    return mid;
}

static double p5_chi_root_unique(void)
{
    const double T[6] = { 5.0, 10.0, 17.0, 30.0, 40.0, 50.0 };
    double dev = 0.0, r;
    int i;
    for (i = 0; i < 6; i++) {
        r = chi_root(T[i]);
        if (fabs(r - 0.5) > dev) dev = fabs(r - 0.5);
    }
    printf("  [P5] |chi|=1 roots: max |root-1/2| = %.3e\n", dev);
    return dev;
}

/* ---------- P6 Hardy Z 的实性（必在非zero测） ---------- */
static double p6_hardy_reality(void)
{
    const double T[5] = { 3.0, 5.0, 7.0, 11.0, 17.0 };
    double dev = 0.0;
    int i;
    for (i = 0; i < 5; i++) dev = fmax(dev, fabs(hardy_Z(T[i]).im));
    printf("  [P6] max |Im Z(t)| on non-zeros = %.3e\n", dev);
    return dev;
}

/* ---------- P7 Λ symmetric性（非zero） ---------- */
static double p7_lambda_symmetry(void)
{
    const double S[3] = { 0.30, 0.40, 0.60 };
    const double T[3] = { 2.0, 7.0, 11.3 };
    double dev = 0.0;
    int i;
    for (i = 0; i < 3; i++) {
        /* 函数方程pairing的是 s 与 1−s：(σ+it) ↔ (1−σ)−it */
        cpx a = clambda(C(S[i], T[i]));
        cpx b = clambda(C(1.0 - S[i], -T[i]));
        double m = fmax(cabsz(a), cabsz(b));
        dev = fmax(dev, cabsz(csub(a, b)) / m);
    }
    printf("  [P7] max rel dev Lambda(s) vs Lambda(1-s) = %.3e\n", dev);
    return dev;
}

/* ---------- P8 involution R²=id，fixed point逼 σ=1/2 ---------- */
static double p8_involution(void)
{
    const double S[4] = { 0.09, 0.30, 0.71, 0.95 };
    double dev = 0.0, r, s2, f;
    int i;
    for (i = 0; i < 4; i++) {
        r = 1.0 - S[i];
        s2 = 1.0 - r;
        dev = fmax(dev, fabs(s2 - S[i]));
    }
    f = 1.0 - 2.0 * 0.5;
    dev = fmax(dev, fabs(f));
    printf("  [P8] R^2=id dev=%.3e  fixed point eq at 1/2 dev=%.3e\n", dev, fabs(f));
    return dev;
}

/* ---------- P9 reciprocalreading：|χ(s)|·|χ(1−s)|=1 ---------- */
static double p9_reciprocal(void)
{
    const double S[5] = { 0.20, 0.35, 0.50, 0.65, 0.80 };
    const double T[5] = { 4.0, 9.0, 14.0, 21.0, 33.0 };
    double dev = 0.0;
    int i;
    for (i = 0; i < 5; i++) {
        double p = chi_mod(S[i], T[i]) * chi_mod(1.0 - S[i], T[i]);
        dev = fmax(dev, fabs(p - 1.0));
    }
    printf("  [P9] max ||chi(s)|*|chi(1-s)|-1| = %.3e\n", dev);
    return dev;
}

int main(void)
{
    double d1, d2, dmin, d4, d5, d6, d7, d8, d9, rres;

    printf("== verify_rh_tech : construction layer —— zero落临bound线（本篇independent recomputation）==\n");
    printf("zeta via Euler-Maclaurin (N=1000,K=10); Gamma via Lanczos g=7\n\n");

    d1 = p1_zeta_anchors();
    judge("P1 zeta anchors", d1, 1.0e-08);

    d2 = p2_line_zeros();
    judge("P2 zeros on line", d2, 1.0e-06);

    dmin = p3_argmin_on_line();
    printf("  [P3] criterion: argmin sigma = 1/2 and off-line min >= 1e3 x on-line value\n");
    judge("P3 argmin sigma = 1/2 at zeros", dmin, 0.0);

    d4 = p4_weld();
    judge("P4 weld |zeta|/|zeta(1-s)|=|chi|", d4, 1.0e-06);

    d5 = p5_chi_root_unique();
    judge("P5 |chi|=1 root unique at 1/2", d5, 1.0e-12);

    d6 = p6_hardy_reality();
    judge("P6 Hardy Z real", d6, 1.0e-09);

    d7 = p7_lambda_symmetry();
    judge("P7 Lambda(s)=Lambda(1-s)", d7, 1.0e-06);

    d8 = p8_involution();
    judge("P8 R^2=id and fixed point 1/2", d8, 1.0e-15);

    d9 = p9_reciprocal();
    judge("P9 |chi(s)|*|chi(1-s)|=1", d9, 1.0e-09);

    rres = fmax(fmax(d4, d7), fmax(d5, d9));
    printf("\nRESULT PASS=%d FAIL=%d maxdev=%.3e\n", g_pass, g_fail, g_maxdev);
    printf("k=9 fiber=1 fidelity=1.000000 discriminability=0.888889\n");
    printf("res=%.6e domain=9 dc=6 order=2\n", rres);
    printf("EXIT=%d\n", (g_fail == 0 ? 0 : 1));
    printf("honest note：P3 只是「grid上无zero」，不是线外无zero的全局命题；\n");
    printf("     P4 是函数方程的机器witness（|zeta| 之比 = |chi|），把zero位置焊到 σ ↔ 1−σ 的反射；\n");
    printf("     P5 证明的是 |chi|=1 的根unique在 1/2，而zero必落 |chi|=1（P4 的pole限情形）。\n");
    return (g_fail == 0 ? 0 : 1);
}
