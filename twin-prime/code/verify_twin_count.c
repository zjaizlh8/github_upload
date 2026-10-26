/*
 * NOTE: comments and output labels in this file were translated from Chinese.
 * The original (Chinese) source is the machine witness in the D:\\dpgo archive;
 * the translation keeps line structure, check names and all numeric values unchanged.
 */
/* verify_twin_count.c —— 孪生素数对计数independent recomputation（本篇论文 §6 表的数据见证）
 * caliber：pi2(x) = #{ p : p, p+2 皆素, p+2 <= x }
 * 与archivecalibercross-check：6m 中心者（排除三连体 (3,5),(5,7)）另计一列。
 * 纯 C，逐点筛法（位集），无外部表、无拟合。
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

/* 传统 Eratosthenes（位集）—— 只在 main 里调一次，循环封在函数内 */
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

/* 分段统计：返回 (x 以内孪生成对数, 其中 6m 中心者) */
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

/* 外部锚：五档 pi2(x) 的公认值（不来自本程序，来自独立文献/独立实现的公开计数） */
static const unsigned long long KNOWN_X[5]   = {100000ULL, 1000000ULL, 10000000ULL,
                                                100000000ULL, 1000000000ULL};
static const unsigned long long KNOWN_PI2[5] = {1224ULL, 8169ULL, 58980ULL,
                                                440312ULL, 3424506ULL};

/* 一次扫描把五档结果装进数组（循环只做填充，不判断不打印） */
static void fill_counts(unsigned long long *pi2, unsigned long long *six)
{
    int i;
    unsigned long long a, s;
    for (i = 0; i < 5; i++) {
        count_twins(KNOWN_X[i], &a, &s);
        pi2[i] = a; six[i] = s;
    }
}

static void print_table(const unsigned long long *pi2, const unsigned long long *six)
{
    int i;
    printf("TABLE 孪生素数对计数（独立复算 · 纯 C 位集 sieve）\n");
    printf("口径 pi2(x)=#{p: p,p+2 皆素, p+2<=x}；six=其中中心 6m 者\n");
    printf("%-12s %-14s %-14s %-10s\n", "x", "pi2(x)", "six(6m中心)", "差(非6m中心)");
    for (i = 0; i < 5; i++)
        printf("%-12llu %-14llu %-14llu %-10llu\n",
               KNOWN_X[i], pi2[i], six[i], pi2[i] - six[i]);
    printf("LIMIT=%llu\n", LIMIT);
}

/* 判据（6 项）：五档 pi2 对公认值 + 一条结构律「非 6m 中心者恰 1 个 = (3,5)」 */
static int judge_counts(const unsigned long long *pi2, const unsigned long long *six)
{
    int i, bad = 0, pass = 0;
    for (i = 0; i < 5; i++) {
        int ok = (pi2[i] == KNOWN_PI2[i]);
        printf("  [%s] P%-2d pi2(%llu) = %llu （公认 %llu）\n",
               ok ? "PASS" : "FAIL", i + 1, KNOWN_X[i], pi2[i], KNOWN_PI2[i]);
        if (ok) pass++; else bad++;
    }
    {
        int ok = 1;
        for (i = 0; i < 5; i++)
            if (pi2[i] - six[i] != 1ULL) ok = 0;
        printf("  [%s] P6   非 6m 中心者恒 1（唯 (3,5)）· 五档一致\n", ok ? "PASS" : "FAIL");
        if (ok) pass++; else bad++;
    }
    printf("RESULT PASS=%d FAIL=%d\n", pass, bad);
    return bad;
}

int main(void)
{
    unsigned long long pi2[5], six[5];
    int bad;
    build_sieve(LIMIT);
    fill_counts(pi2, six);           /* 一趟扫描：只填数组 */
    print_table(pi2, six);           /* 一次取出：只打印 */
    bad = judge_counts(pi2, six);    /* 一次取出：只判据 */
    printf("EXIT=%d\n", bad ? 1 : 0);
    free(sieve);
    return bad ? 1 : 0;
}
