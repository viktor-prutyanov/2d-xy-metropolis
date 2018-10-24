#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cerrno>
#include <getopt.h>
#include <iostream>
#include <cmath>
#include <climits>

static long int L = 0;
static unsigned int smpl = 0;
static uint8_t *s = NULL;

#define S(i, j) (s[i * L + j])
#define NR_STATES (1 << UINT8_WIDTH)
#define E_pair(s1, s2) (cos(2 * (s1 * M_PI - s2 * M_PI) / NR_STATES))

inline uint8_t S_p(int i, int j)
{
    int I = (i < L) ? ((i >= 0) ? i : (L + (i % L))) : (i % L);
    int J = (j < L) ? ((j >= 0) ? j : (L + (j % L))) : (j % L);

    return S(I, J);
}

static double T = 1.0;
static unsigned long int nr_iters = 0;
static double E = 0;

unsigned long int get_pos_ulong_opt(char *optarg, const char *name)
{
    char *endptr;
    unsigned long int arg = strtoul(optarg, &endptr, 10);

    errno = 0;

    if (errno == ERANGE) {
        fprintf(stderr, "args error: %s is out of range\n", name);
    } else if (*endptr != '\0') {
        fprintf(stderr, "args error: invalid %s\n", name);
        errno = EINVAL;
    } else if (arg <= 0) {
        fprintf(stderr, "args error: nonpositive %s\n", name);
        errno = EINVAL;
    }

    return arg;
}

double get_pos_double_opt(char *optarg, const char *name)
{
    char *endptr;
    double arg = strtod(optarg, &endptr);

    errno = 0;

    if (errno == ERANGE) {
        fprintf(stderr, "args error: %s is out of range.\n", name);
    } else if (*endptr != '\0') {
        fprintf(stderr, "args error: invalid %s.\n", name);
        errno = EINVAL;
    } else if (arg <= 0) {
        fprintf(stderr, "args error: nonpositive %s.\n", name);
        errno = EINVAL;
    }

    return arg;
}

double get_m(uint8_t *s)
{
    double mx = 0;
    double my = 0;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            mx += cos(2 * M_PI * S(i, j) / NR_STATES);
            my += sin(2 * M_PI * S(i, j) / NR_STATES);
        }
    }

    return hypot(mx, my);
}

double get_E(uint8_t *s)
{
    double E = 0;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            E += cos(2 * (S(i, j) * M_PI - S_p(i + 1, j) * M_PI) / NR_STATES) +
                cos(2 * (S(i, j) * M_PI - S_p(i, j + 1) * M_PI) / NR_STATES);
        }
    }

    return -E;
}

/*
double get_E2(int8_t *s)
{
    double E = get_E(s) * 1.0 / (L * L);
    double E2 = 0;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            double e = pow(E + S(i, j) * (S_p(i + 1, j) + S_p(i - 1, j) + S_p(i, j + 1) + S_p(i, j - 1)) / 2, 2);
            E2 += e;
        }
    }

    return E2 * 1.0 / (L * L);
}
*/

double get_dE(uint8_t *s, uint8_t ds, int i, int j)
{
    uint8_t s0 = S(i, j);
    uint8_t s1 = s0 + ds;

    return (E_pair(s0, S_p(i + 1, j)) + E_pair(s0, S_p(i - 1, j)) +
            E_pair(s0, S_p(i, j + 1)) + E_pair(s0, S_p(i, j - 1))) -
        (E_pair(s1, S_p(i + 1, j)) + E_pair(s1, S_p(i - 1, j)) +
         E_pair(s1, S_p(i, j + 1)) + E_pair(s1, S_p(i, j - 1)));
}

double get_C_v(uint8_t *s, double T)
{
    return 0.0;
//    return get_E2(s) / (T * T);
}

void dump(uint8_t *s, const char *name)
{
    FILE *file;

    if (!name) {
        return;
    }

    file = fopen(name, "w");
    if (!file) {
        return;
    }

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L - 1; j++) {
            fprintf(file, "%u ", S(i, j));
        }
        fprintf(file, "%u\n", S(i, L - 1));
    }

    fclose(file);
}

int main(int argc, char *argv[])
{
    int opt = 0;
    const char usage_str[] = "usage: %s [-L size] [-h help] [-n iterations]"
                             " [-T temperature] [-s sampling] [-o out_file]\n";
    char *out_filename = NULL;

    while ((opt = getopt(argc, argv, "L:n:T:o:s:h")) != -1) {
        switch(opt) {
        case 'L':
            L = get_pos_ulong_opt(optarg, "size of lattice side");
            break;
        case 'n':
            nr_iters = get_pos_ulong_opt(optarg, "number of iterations");
            break;
        case 'T':
            T = get_pos_double_opt(optarg, "temperature");
            break;
        case 'o':
            out_filename = strdup(optarg);
            break;
        case 's':
            smpl = get_pos_ulong_opt(optarg, "sampling period");
            break;
        case 'h':
        default:
            fprintf(stderr, usage_str, argv[0]);
            return -1;
            break;
        }
    }

    if (L < 1) {
        return -1;
    }

    if (!out_filename) {
        out_filename = (char *)malloc(256);
        sprintf(out_filename, "%lf_%ld.txt", T, L);
    }

    s = (uint8_t *)malloc(L * L * sizeof(*s));

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            S(i, j) = rand() % NR_STATES;
        }
    }

    E = get_E(s);
    //printf("%d.0,%lf\n", -1, E);

    for (size_t iter = 0; iter < nr_iters; iter++) {
        int k = rand() % (L * L), i = k / L, j = k % L;
        uint8_t ds = 1;
        double dE = get_dE(s, ds, i, j);

        if ((dE <= 0) || (rand() * 1.0 / RAND_MAX < exp(-dE / T))) {
            S(i, j) += ds;
            E += dE;
        }

        if (smpl && !((iter + 1) % smpl)) {
            printf("%lu.0,%lf\n", iter, E / (L * L));
        }
    }

    //printf("%lf,%lf,%lf,%lf\n", T, E * 1.0 / (L * L), get_m(s) * 1.0 / (L * L), get_C_v(s, T));

    if (!smpl) {
        printf("%lf,%lf,%lf\n", T, E / (L * L), get_m(s) / (L * L));
    }
    dump(s, out_filename);

    free(s);
    free(out_filename);

    return 0;
}
