#include <cmath>
#include <fftw3.h>        //БДПФ
#include <iostream>

using namespace std;

double f (double x, double y, double z);
double u (double x, double y, double z);
void   change_F (double* G, double* F, unsigned int N);
double norm (double* G, double* U, unsigned int N);

int main()
{
    unsigned int N;        // количество точек апроксимации
    unsigned int m;        // временное хранение, площадь сетки
    int n;                 // размерность вектора

    double* F;        // численное решение
    double* U;        // аналитическое решение
    double* G;        // временное хранение

    fftw_plan p;        // настройки ДБПФ

    N = 5;
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

    fftw_execute (p);          // выполнение прямого ДБПФ
    change_F (G, F, N);        // решение ЛАУ
    fftw_execute (p);          // выполнение обратного ДБПФ

    cout << "N = " << N << ",\tnorm = " << norm (G, U, N) << endl;

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

double u (double x, double y, double z) { return sin (x) * sin (y) * sin (z); }

void change_F (double* G, double* F, unsigned int N)
{
    unsigned int m;
    double*      c;

    m = static_cast<unsigned int> (pow (N - 1, 2));
    c = new double[N - 1];

    for (unsigned int i = 1; i < N; i++)
        c[i - 1] = pow (2 * N * sin (i * M_PI / (2 * N)) / M_PI, 2);

    for (unsigned int i = 0; i < N - 1; i++)
    {
        for (unsigned int j = 0; j < N - 1; j++)
        {
            for (unsigned int k = 0; k < N - 1; k++)
                F[i * m + j * (N - 1) + k] =
                    G[i * m + j * (N - 1) + k] /
                    (pow (2 * N, 3) * (c[i] + c[j] + c[k]));
        }
    }
}

double norm (double* G, double* U, unsigned int N)
{
    double norm;
    norm = static_cast<double> (0);

    for (unsigned int i = 0; i < pow (N - 1, 3); i++)
        norm += pow (U[i] - G[i], 2);

    return sqrt (norm / pow (N + 1, 3));
}