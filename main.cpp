#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cerrno>
#include <getopt.h>
#include <iostream>
#include <cmath>
#include <climits>
#include <vector>

//#define MKL
#define NR_RND 1024

#ifdef MKL
#include "mkl_vsl.h"
#endif

static long int L = 0;
static unsigned int smpl = 0;
static double *s = NULL;

#define S(i, j) (s[i * L + j])
#define E_pair(s1, s2) (cos(s1 - s2))

inline double S_p(int i, int j)
{
    int I = (i < L) ? ((i >= 0) ? i : (L + (i % L))) : (i % L);
    int J = (j < L) ? ((j >= 0) ? j : (L + (j % L))) : (j % L);

    return S(I, J);
}

static double T = 1.0;
static size_t nr_iters = 0;
static size_t a = 0;
static double E = 0;
static double E2 = 0;

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

double get_mx(double *s)
{
    double mx = 0;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            mx += cos(S(i, j));
        }
    }

    return mx;
}

double get_my(double *s)
{
    double my = 0;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            my += sin(S(i, j));
        }
    }

    return my;
}

double get_E(double *s)
{
    double e = 0;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            e += E_pair(S(i, j), S_p(i + 1, j)) + E_pair(S(i, j), S_p(i, j + 1));
        }
    }

    return -e;
}

double get_E_s(double *s, double s0, int i, int j)
{
    return -(E_pair(s0, S_p(i + 1, j)) + E_pair(s0, S_p(i - 1, j)) + E_pair(s0, S_p(i, j + 1)) + E_pair(s0, S_p(i, j - 1))) / 2;
}

/*
double get_dE2(double *s, double s_prev, int i, int j)
{
    double s0 = S(i, j);

    return pow(get_E_s(s, s0, i, j), 2) - pow(get_E_s(s, s_prev, i, j), 2);
}

double get_E2(double *s)
{
    double E2 = 0;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            double s0 = S(i, j);

            E2 += pow(get_E_s(s, s0, i, j), 2);
        }
    }

    return E2;
}
*/
double get_dE(double *s, double s1, int i, int j)
{
    double s0 = S(i, j);

    return (E_pair(s0, S_p(i + 1, j)) + E_pair(s0, S_p(i - 1, j)) +
            E_pair(s0, S_p(i, j + 1)) + E_pair(s0, S_p(i, j - 1))) -
        (E_pair(s1, S_p(i + 1, j)) + E_pair(s1, S_p(i - 1, j)) +
         E_pair(s1, S_p(i, j + 1)) + E_pair(s1, S_p(i, j - 1)));
}

double get_Cv(double T, double E, double E2)
{
    return (E2 / (L * L) - E * E / (L * L * L * L)) / (T * T);
}

void dump(double *s, const char *name)
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
            fprintf(file, "%.3lf ", S(i, j));
        }
        fprintf(file, "%.3lf\n", S(i, L - 1));
    }

    fclose(file);
}

int main(int argc, char *argv[])
{
    int opt = 0;
    const char usage_str[] = "usage: %s [-L size] [-h help] [-n iterations]"
                             " [-T temperature] [-s sampling] [-o out_file]"
                             " [-a avg_msrmnts]\n";
    char *out_filename = NULL;
#ifdef MKL
    VSLStreamStatePtr stream;
#endif
    double sum_E = 0;
    double sum_m = 0;
    double mx = 0, my = 0;

    while ((opt = getopt(argc, argv, "L:n:T:o:s:a:h")) != -1) {
        switch(opt) {
        case 'L':
            L = get_pos_ulong_opt(optarg, "size of lattice side");
            break;
        case 'n':
            nr_iters = get_pos_ulong_opt(optarg, "number of iterations");
            break;
        case 'a':
            a = get_pos_ulong_opt(optarg, "number of iterations for"
                    " average measurements");
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

    if ((L < 1) || (a > nr_iters)) {
        return -1;
    }

    if (!out_filename) {
        out_filename = (char *)malloc(256);
        sprintf(out_filename, "%lf_%ld.txt", T, L);
    }

    s = (double *)malloc(L * L * sizeof(*s));

#ifdef MKL
    vslNewStream(&stream, VSL_BRNG_RDRAND, 1);
    //vslNewStream(&stream, VSL_BRNG_MT19937, 1);
    double *d_rnd_buf1 = (double *)malloc(NR_RND * sizeof(*d_rnd_buf1));
    double *d_rnd_buf2 = (double *)malloc(NR_RND * sizeof(*d_rnd_buf2));
    int *i_rnd_buf = (int *)malloc(NR_RND * sizeof(*i_rnd_buf));
#endif
    //srand(time(NULL));

#ifdef MKL
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, L * L, s, 0, 2 * M_PI);
#else
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            S(i, j) = (rand() * 1. / RAND_MAX) * 2 * M_PI;
        }
    }
#endif

    std::vector<double> Es, ms;

    E = get_E(s);
    mx = get_mx(s);
    my = get_my(s);

    for (size_t iter = 0; iter < nr_iters; iter++) {
#ifdef MKL
        if (iter % NR_RND == 0) {
            viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, NR_RND, i_rnd_buf, 0, L * L);
            vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, NR_RND, d_rnd_buf1, - M_PI * 0.5, M_PI * 0.5);
            vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, NR_RND, d_rnd_buf2, 0, 1);
        }
        int k = i_rnd_buf[iter % NR_RND];
#else
        int k = rand() % (L * L);
#endif
        int i = k / L, j = k % L;
#ifdef MKL
        double ds = d_rnd_buf1[iter % NR_RND];
#else
        double ds = (rand() * 1. / RAND_MAX - 0.5) * M_PI;
#endif
        double s0 = S(i, j);
        double s1 = S(i, j) + ds;
        double dE;

        if (s1 < 0) {
            s1 += 2 * M_PI;
        } else if (s1 >= 2 * M_PI) {
            s1 -= 2 * M_PI;
        }

        dE = get_dE(s, s1, i, j);

#ifdef MKL
        if ((dE <= 0) || (d_rnd_buf2[iter % NR_RND] < exp(-dE / T))) {
#else
        if ((dE <= 0) || (rand() * 1.0 / RAND_MAX < exp(-dE / T))) {
#endif
            S(i, j) = s1;
            E += dE;
            if (iter > (nr_iters - a)) {
                mx += cos(s1) - cos(s0);
                my += sin(s1) - sin(s0);
            }
        }

        if (iter == (nr_iters - a)) {
            mx = get_mx(s);
            my = get_my(s);
        }

        if (iter > (nr_iters - a)) {
            sum_m += hypot(mx, my);
            sum_E += E;
            Es.push_back(E);
            ms.push_back(hypot(mx, my));
        }

        if (smpl && !((iter + 1) % smpl)) {
            printf("%lu.0,%lf\n", iter, E / (L * L));
        }
    }

    double E_avg = sum_E / a / (L * L);
    double E_var = 0;
    for (auto E : Es) {
        E_var += pow(E_avg - E / (L * L), 2);
    }
    E_var /= a;

    double m_avg = sum_m / a / (L * L);
    double m_var = 0;
    for (auto m : ms) {
        m_var += pow(m_avg - m / (L * L), 2);
    }
    m_var /= a;

    double m4_avg = 0;
    double m2_avg = 0;
    for (auto m : ms) {
        m4_avg += pow(m / (L * L), 4);
        m2_avg += pow(m / (L * L), 2);
    }
    double U = 1 - m4_avg * a / (2 * m2_avg * m2_avg);

    if (!smpl) {
        printf("%lf,%lf,%lf,%lf,%lf,%lf,", T, sum_E / a / (L * L),
                sum_m / a / (L * L), E_var / (T * T), m_var / T, U);
    }
    dump(s, out_filename);

    free(s);
    free(out_filename);

#ifdef MKL
    vslDeleteStream(&stream);
    free(d_rnd_buf1);
    free(d_rnd_buf2);
    free(i_rnd_buf);
#endif

    return 0;
}
