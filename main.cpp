#include <cmath>
#include <cstdio>
#include <fftw3.h>        //БДПФ
#include <iostream>
#include <stdlib.h>

using namespace std;
#define Pi static_cast<long double> (3.14159265358979323846)
unsigned int N = 4;

long double f (long double x, long double y, long double z);
long double u (long double x, long double y, long double z);
void        change_F (long double* G, long double* F);
long double norm (long double* G, long double* U);

int main()
{
    unsigned int m;
    int          n;
    long double* F;
    long double* U;
    long double* G;
    fftwl_plan   p;

    n = static_cast<int> (N - 1);
    m = static_cast<unsigned int> (pow (n, 3));

    G = static_cast<long double*> (fftw_malloc (sizeof (long double) * m));
    F = static_cast<long double*> (fftw_malloc (sizeof (long double) * m));
    U = static_cast<long double*> (fftw_malloc (sizeof (long double) * m));
    
    m = static_cast<unsigned int> (pow (n, 2));

    for (unsigned int i = 1; i < N; i++)
    {
        for (unsigned int j = 1; j < N; j++)
        {
            for (unsigned int k = 1; k < N; k++)
            {
                F[(i - 1) * m + (j - 1) * (N - 1) + k - 1] =
                    f (Pi * i / N, Pi * j / N, Pi * k / N);
                U[(i - 1) * m + (j - 1) * (N - 1) + k - 1] =
                    u (Pi * i / N, Pi * j / N, Pi * k / N);
            }
        }
    }

    p = fftwl_plan_r2r_3d (n, n, n, F, G, FFTW_RODFT00, FFTW_RODFT00,
                           FFTW_RODFT00, FFTW_ESTIMATE);

    fftwl_execute (p);
    change_F (G, F);
    fftwl_execute (p);

    printf ("norm = %Lf\n", norm (G, U));

    fftwl_destroy_plan (p);

    fftwl_free (F);
    fftwl_free (G);
    fftwl_free (U);

    return 0;
}

long double f (long double x, long double y, long double z)
{
    return 3 * sinl (x) * sinl (y) * sinl (z);
}

long double u (long double x, long double y, long double z)
{
    return sinl (x) * sinl (y) * sinl (z);
}

void change_F (long double* G, long double* F)
{
    unsigned int m;
    long double* c = new long double[N - 1];
    m              = static_cast<unsigned int> (powl (N - 1, 2));
    for (unsigned int i = 1; i < N; i++)
        c[i - 1] = powl (sinl (i * Pi / (2 * N)), 2);

    for (unsigned int i = 0; i < N - 1; i++)
    {
        for (unsigned int j = 0; j < N - 1; j++)
        {
            for (unsigned int k = 0; k < N - 1; k++)
                F[i * m + j * (N - 1) + k] =
                    G[i * m + j * (N - 1) + k] * powl (Pi, 2) /
                    (powl (2 * N, 5) * (c[i] + c[j] + c[k]));
        }
    }
}

long double norm (long double* G, long double* U)
{
    long double norm = static_cast<long double> (0);
    for (unsigned int i = 0; i < powl (N - 1, 3); i++)
        norm += powl (U[i] - G[i], 2);
    norm = sqrtl (norm / powl (N + 1, 3));
    return norm;
}