#include <cmath>
#include <cstdio>
#include <fftw3.h>        //БДПФ
#include <iostream>
#include <stdlib.h>

using namespace std;

unsigned int N = 5;

double f (double x, double y, double z);
double u (double x, double y, double z);
void        change_F (double* G, double* F);
double norm (double* G, double* U);

int main()
{
    unsigned int m;
    int          n;
    double* F;
    double* U;
    double* G;
    fftw_plan p;

    n = static_cast<int> (N - 1);
    m = static_cast<unsigned int> (pow (n, 3));

    G = static_cast<double*> (fftw_malloc (sizeof (double) * m));
    F = static_cast<double*> (fftw_malloc (sizeof (double) * m));
    U = static_cast<double*> (fftw_malloc (sizeof (double) * m));

    m = static_cast<unsigned int> (pow (n, 2));

    for (unsigned int i = 1; i < N; i++)
    {
        for (unsigned int j = 1; j < N; j++)
        {
            for (unsigned int k = 1; k < N; k++)
            {
                F[(i - 1) * m + (j - 1) * (N - 1) + k - 1] =
                    f (M_PI * i / N, M_PI * j / N, M_PI * k / N);
                U[(i - 1) * m + (j - 1) * (N - 1) + k - 1] =
                    u (M_PI * i / N, M_PI * j / N, M_PI * k / N);
            }
        }
    }

    p = fftw_plan_r2r_3d (n, n, n, F, G, FFTW_RODFT00, FFTW_RODFT00,
                           FFTW_RODFT00, FFTW_ESTIMATE);

    fftw_execute (p);
    change_F (G, F);
    fftw_execute (p);

    printf ("N = %u,\tnorm = %f\n", N, norm (G, U));

    fftw_destroy_plan (p);

    fftw_free (F);
    fftw_free (G);
    fftw_free (U);

    return 0;
}

double f (double x, double y, double z)
{
    return 3 * sin (x) * sin (y) * sin (z);
}

double u (double x, double y, double z)
{
    return sin (x) * sin (y) * sin (z);
}

void change_F (double* G, double* F)
{
    unsigned int m;
    double* c = new double[N - 1];
    m              = static_cast<unsigned int> (pow (N - 1, 2));
    for (unsigned int i = 1; i < N; i++)
        c[i - 1] = pow (sin (i * M_PI / (2 * N)), 2);

    for (unsigned int i = 0; i < N - 1; i++)
    {
        for (unsigned int j = 0; j < N - 1; j++)
        {
            for (unsigned int k = 0; k < N - 1; k++)
                F[i * m + j * (N - 1) + k] =
                    G[i * m + j * (N - 1) + k] * pow (M_PI, 2) /
                    (pow (2 * N, 5) * (c[i] + c[j] + c[k]));
        }
    }
}

double norm (double* G, double* U)
{
    double norm = static_cast<double> (0);
    for (unsigned int i = 0; i < pow (N - 1, 3); i++)
        norm += pow (U[i] - G[i], 2);
    norm = sqrt (norm / pow (N + 1, 3));
    return norm;
}