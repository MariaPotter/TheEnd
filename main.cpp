#include <cmath>
#include <fftw3.h>        //БДПФ
#include <iostream>

using namespace std;

double f (double G, double y, double z);
double u (double G, double y, double z);
void   Runthrough_Method (double* in, double* out, unsigned int N);        //Прогонка
double norm (double* G, double* U, unsigned int N);

int main (int argc, char* argv[])
{
    unsigned int N;        // количество точек апроксимации
    unsigned int m;        // временное хранение, площадь сетки
    int n;        // размерность вектора

    double* F;        // численное решение
    double* U;        // аналитическое решение
    double* G;        // временное хранение

    fftw_plan p;        // настройки ДБПФ

    N = 128;
    if (argc > 1) sscanf (argv[1], "%u", &N);

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

    fftw_execute (p);        // выполнение прямого ДБПФ
    Runthrough_Method (G, F, N);        // решение ЛАУ
    fftw_execute (p);        // выполнение обратного ДБПФ

    cout << "N = " << N << ",   norm = " << norm (G, U, N) << endl;

    fftw_destroy_plan (p);

    fftw_free (F);
    fftw_free (G);
    fftw_free (U);

    return 0;
}

double f (double G, double y, double z)
{
    return 3 * sin (G) * sin (y) * sin (z);
}

void change_F (double* in, double* out, unsigned int N);

double u (double G, double y, double z) { return sin (G) * sin (y) * sin (z); }

void Runthrough_Method (double* in, double* out, unsigned int N)
{
    double* y = new double[N];
    double* q = new double[N];
    double* p = new double[N];
    double* l = new double[N];
    double* z = new double[N];

    double* a = new double[N];
    double* b = new double[N];
    double* c = new double[N];
    double* d = new double[N];
    double* e = new double[N];
    double* f = new double[N];

    for (unsigned int i = 0; i < N; i++)
    {
        if (i > 1) e[i] = in[i * N + i - 2];
        else e[i] = 0;

        if (i > 0) c[i] = in[i * N + i - 1];
        else c[i] = 0;

        d[i] = in[i * (1 + N)];

        if (i < N - 1) a[i] = in[i * N + i + 1];
        else a[i] = 0;

        if (i < N - 2) b[i] = in[i * N + i + 2];
        else b[i] = 0;

        y[i] = 0;
        q[i] = 0;
        p[i] = 0;
        l[i] = 0;
        z[i] = 0;
        f[i] = in[i * N + N + 1];
    }
    y[0] = 0;
    q[0] = d[0];
    p[0] = b[0] / q[0];
    ;
    l[0] = a[0] / q[0];
    z[0] = f[0] / q[0];

    y[1] = c[1];
    q[1] = d[1] - l[0] * y[1];
    p[1] = b[1] / q[1];
    l[1] = (a[1] - p[0] * y[1]) / q[1];
    z[1] = (f[1] - z[0] * y[1]) / q[1];

    for (unsigned int i = 2; i < N - 2; i++)
    {
        y[i] = c[i] - l[i - 2] * e[i];
        q[i] = d[i] - p[i - 2] * e[i] - l[i - 1] * y[i];
        p[i] = b[i] / q[i];
        l[i] = (a[i] - p[i - 1] * y[i]) / q[i];
        z[i] = (f[i] - z[i - 2] * e[i] - z[i - 1] * y[i]) / q[i];
    }

    y[N - 2] = c[N - 2] - l[N - 4] * e[N - 2];
    q[N - 2] = d[N - 2] - p[N - 4] * e[N - 2] - l[N - 3] * y[N - 2];
    p[N - 2] = 0;
    l[N - 2] = (a[N - 2] - p[N - 3] * y[N - 2]) / q[N - 2];
    z[N - 2] =
        (f[N - 2] - z[N - 3] * e[N - 2] - z[N - 3] * y[N - 2]) / q[N - 2];

    y[N - 1] = c[N - 1] - l[N - 3] * e[N - 1];
    q[N - 1] = d[N - 1] - p[N - 3] * e[N - 1] - l[N - 2] * y[N - 1];
    p[N - 1] = 0;
    l[N - 1] = 0;
    z[N - 1] =
        (f[N - 1] - z[N - 2] * e[N - 1] - z[N - 2] * y[N - 1]) / q[N - 1];

    out[N - 1] = z[N - 1];
    out[N - 2] = z[N - 2] - l[N - 2] * out[N - 1];

    for (int i = static_cast<int> (N) - 3; i > -1; i--)
    {
        out[i] = z[i] - l[i] * out[i + 1] - p[i] * out[i + 2];
    }
    out[N - 1] = 0;

    change_F (in, out, N);

    delete[] y;
    delete[] q;
    delete[] p;
    delete[] l;
    delete[] z;

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] e;
}

double norm (double* G, double* U, unsigned int N)
{
    double norm;
    norm = static_cast<double> (0);

    for (unsigned int i = 0; i < pow (N - 1, 3); i++)
        norm += pow (U[i] - G[i], 2);

    return sqrt (norm / pow (N + 1, 3));
}

void change_F (double* in, double* out, unsigned int N)
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
                out[i * m + j * (N - 1) + k] =
                    in[i * m + j * (N - 1) + k] /
                    (pow (2 * N, 3) * (c[i] + c[j] + c[k]));
        }
    }
}
