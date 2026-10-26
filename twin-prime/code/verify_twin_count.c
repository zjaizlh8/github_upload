/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* verify_twin_count.c —— 孪generateprime对计数independent recomputation（本篇论文 §6 表的数据witness）
 * caliber：pi2(x) = #{ p : p, p+2 皆素, p+2 <= x }
 * 与archivecalibercross-check：6m 中心者（排除三连体 (3,5),(5,7)）另计一列。
 * 纯 C，逐点sieve（位集），无外部表、无拟合。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LIMIT 1000000000ULL  /* 1e9 */

static unsigned char *sieve; /* 位集：1 位 1 数 */

static int is_prime(unsigned long long n)
{
    if (n < 2ULL) return 0;
    if (n == 2ULL) return 1;
    if ((n & 1ULL) == 0ULL) return 0;
    return (sieve[n >> 3] >> (n & 7ULL)) & 1U ? 1 : 0;
}

/* 传统 Eratosthenes（位集）—— 只在 main 里调一次，cycle封在函数内 */
static void build_sieve(unsigned long long n)
{
    unsigned long long i, j;
    unsigned long long bytes = (n >> 3) + 2ULL;
    sieve = (unsigned char *)malloc((size_t)bytes);
    if (!sieve) { printf("MALLOC_FAIL\n"); exit(1); }
    memset(sieve, 0xFF, (size_t)bytes);
    sieve[0] &= (unsigned char)~0x03U;      /* 0,1 非素 */
    for (i = 2ULL; i * i <= n; i++) {
        if (!is_prime(i)) continue;
        for (j = i * i; j <= n; j += i)
            sieve[j >> 3] &= (unsigned char)~(1U << (j & 7ULL));
    }
}

/* 分段统计：返回 (x 以内孪generate对数, 其中 6m 中心者) */
static void count_twins(unsigned long long x, unsigned long long *all, unsigned long long *six)
{
    unsigned long long p;
    unsigned long long a = 0ULL, s = 0ULL;
    for (p = 3ULL; p + 2ULL <= x; p += 2ULL) {
        if (is_prime(p) && is_prime(p + 2ULL)) {
            a++;
            if (((p + 1ULL) % 6ULL) == 0ULL) s++;   /* 中心 6m ⟺ p+1 ≡ 0 mod 6 */
        }
    }
    *all = a; *six = s;
}

static void run_targets(void)
{
    unsigned long long targets[5] = {100000ULL, 1000000ULL, 10000000ULL,
                                     100000000ULL, 1000000000ULL};
    int i;
    unsigned long long tot = 0ULL, six = 0ULL;
    printf("RESULT 孪generateprime对计数independent recomputation（纯 C 位集sieve）\n");
    printf("caliber pi2(x)=#{p: p,p+2 皆素, p+2<=x}；six=其中中心 6m 者\n");
    printf("%-12s %-14s %-14s %-10s\n", "x", "pi2(x)", "six(6m中心)", "差(三连体)");
    for (i = 0; i < 5; i++) {
        unsigned long long x = targets[i];
        count_twins(x, &tot, &six);
        printf("%-12llu %-14llu %-14llu %-10llu\n", x, tot, six, tot - six);
    }
    printf("LIMIT=%llu\n", LIMIT);
}

int main(void)
{
    build_sieve(LIMIT);
    run_targets();
    free(sieve);
    return 0;
}
