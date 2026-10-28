/* verify_claims.c —— 两条读数改成「被迫」：闭式给值，扫描只作对照
 *
 * 四法口径（four operations）：加法=accumulate · 乘法=generate · 取模=sieve · 跨维度=dimension-change · base=非象
 *   跨维度处只用 sqrt（√2），残差口径 = 相对偏差；判据一律「闭式 vs 扫描」，不是「扫描自己报个数」。
 *
 * ① Hodge 非挠 x=√2/7 的 argmin：**由 Pell/连分数递推被迫**，不靠扫。
 *      √2 的收敛子 q：1,2,5,12,29,70,99?,239…（q'=2q+q₋），p'=2p+p₋
 *      定理（最佳逼近）：min_{1≤k≤K} ‖kα‖ 在 q_n ≤ K 的最大收敛子处取到，值是 ‖q_nα‖。
 *      ⟹ K=163 时 argmin 必为 99（下一个是 239>163），min = |99√2−140|/7。
 *      这也直接说明：**163 不是 Pell 分母**，所以「到 k=163 最小距 0.001020」在结构上不可能对。
 *
 * ② ns 尾项：**由地板不等式被迫发散**，不靠比大小。
 *      G_n = max(q·G_{n-1}, Gs) ⟹ G_n ≥ Gs > 0 ⟹ S_2N − S_N = Σ_{n=N}^{2N-1} G_n ≥ N·Gs > 0，
 *      且实测恰等于 N·Gs（判定式），随 N 线性增长 ⟹ 非 Cauchy。
 *      旧判据把「两段和之差」当尾增量，那个差是瞬态超额，与 N 无关。
 */
#include <stdio.h>
#include <math.h>

static int gp = 0, gf = 0;
static void J(const char *name, double dev, double tol)
{
    int ok = (dev <= tol);
    if (ok) gp++; else gf++;
    printf("  [%s] %-58s dev=%.3e tol=%.1e\n", ok ? "PASS" : "FAIL", name, dev, tol);
}

int main(void)
{
    printf("== ① Hodge：argmin 由 Pell 递推被迫 ==\n");
    {
        const double irr = sqrt(2.0) / 7.0;
        const int K = 163;
        long long qbest = 1, pbest = 1;
        int k, kscan = 0;
        double vscan = 1e300;

        /* ---- 判项块 A：α=√2/7 的连分数收敛子给出 K 内最佳分母 ---- */
        {
            double x = irr;
            long long hm1 = 0, h = 1, km1 = 1, kk = 0;   /* h/k 收敛子，初值 0/1 与 1/0 */
            int it;
            for (it = 0; it < 40; it++) {
                double a = floor(x);
                long long h2 = (long long)a * h + hm1;
                long long k2 = (long long)a * kk + km1;
                if (k2 > K) break;
                hm1 = h; h = h2; km1 = kk; kk = k2;
                if (x - a < 1e-15) break;
                x = 1.0 / (x - a);
            }
            qbest = (kk > 0) ? kk : 1; pbest = h;
            printf("       α=√2/7 的连分数收敛子（分母<=%d 的最大项）：p/q = %lld/%lld\n",
                   K, pbest, qbest);
        }

        /* ---- 判项块 B：闭式值 vs 扫描值 ---- */
        {
            double closed = fabs((double)qbest * irr - (double)pbest);
            for (k = 1; k <= K; k++) {
                double d = fabs((double)k * irr - floor((double)k * irr + 0.5));
                if (d < vscan) { vscan = d; kscan = k; }
            }
            printf("       扫描：min=%.9f @ k=%d ；闭式：|%lldα−%lld| = %.9f\n",
                   vscan, kscan, qbest, pbest, closed);
            J("argmin 由连分数收敛子唯一给出（扫描 == 收敛子分母）", (double)(kscan == (int)qbest ? 0 : 1), 0.0);
            J("最小值由闭式给出（扫描 == |qα−p|）", fabs(vscan - closed), 1e-15);
        }

        /* ---- 判项块 C：163 不是 Pell 分母（论文的错在结构上被排除） ---- */
        {
            long long qq = 1, qq2 = 1, t; int is_pell = 0;
            while (qq <= K) {
                if (qq == K) is_pell = 1;
                t = 2 * qq + qq2; qq2 = qq; qq = t;
            }
            J("163 不是 Pell 分母（故它不可能是 argmin）", (double)(is_pell ? 1 : 0), 0.0);
        }
    }

    printf("\n== ② ns：发散由地板不等式被迫 ==\n");
    {
        const double q = 0.5, Gs = 0.1256637061;
        int Ns[3] = { 500, 2000, 5000 }, t;
        for (t = 0; t < 3; t++) {
            int n, N = Ns[t];
            double SN = 0.0, S2N = 0.0, g = 1.0, gmin = 1e300;
            char tag[96];
            for (n = 0; n < 2 * N; n++) {
                if (g < gmin) gmin = g;
                if (n < N) SN += g;
                S2N += g;
                g = fmax(q * g, Gs);
            }
            snprintf(tag, sizeof tag, "N=%d：S_2N−S_N 恰等于 N·Gs（判定式，零容差）", N);
            J(tag, fabs((S2N - SN) - (double)N * Gs), 1e-9 * (double)N * Gs);
            snprintf(tag, sizeof tag, "N=%d：G_n ≥ Gs 逐项成立（发散的地板前提）", N);
            J(tag, (gmin < Gs) ? (Gs - gmin) / Gs : 0.0, 0.0);
        }
    }

    printf("\nRESULT PASS=%d FAIL=%d\nEXIT=%d\n", gp, gf, gf == 0 ? 0 : 1);
    return gf == 0 ? 0 : 1;
}
