#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#include <cmath>

//Точное решение задачи Коши
double exact_solution (double x)
{
    return 2 / (exp(2 * x) + 1);
}

//y'
double f(double x, double y)
{
    return -2*y + y * y;
}

//Функция вычисляющая факториал
int fact(size_t N)
{
    int res = 1;
    for (size_t i = 1; i <= N; i++)
    {
        res *= i;
    }
    return res;
}

/*
* val_taylor - значения решения в точках x_k, k=-2...N (N+3 значений)
*/
void taylor(double val_taylor[], double x[], size_t N, double y0)
{
    double coeffs[7];
    coeffs[0] = y0;
    coeffs[1] = -1;
    coeffs[2] = 0;
    coeffs[3] = 2;
    coeffs[4] = 0;
    coeffs[5] = -16;
    coeffs[6] = 0;

    for (size_t i = 0; i < N + 3; i++)
    {
        val_taylor[i] = 0;
        for (int j = 0; j < 7; j++)
        {
            val_taylor[i] += coeffs[j]*pow(x[i], j)/fact(j);
        }
    }
}

void adams(double adams_values[], double y[], double x[], int N, double h)
{
    double** matrix = new double* [200];
    for (int i = 0; i < 200; i++)
    {
        matrix[i] = new double[200];
    }

    std::copy(y, y + 5, matrix[0]);
    for (int i = 0; i < 5; i++)
    {
        matrix[1][i] = h * f(x[i], matrix[0][i]);
    }

    for (int i = 2; i < 6; i++)
    {
        for (int j = 0; j < 6 - i; j++)
        {
            matrix[i][j] = matrix[i - 1][j + 1] - matrix[i - 1][j];
        }
    }

    for (int i = 1; i < N - 1; i++)
    {
        matrix[0][4 + i] = matrix[0][4 + i - 1] + matrix[1][4 + i - 1] + matrix[2][4 + i - 2] / 2 + 5 * matrix[3][4 + i - 3] / 12 + 3 * matrix[4][4 + i - 4] / 8 + 251 * matrix[5][4 + i - 5] / 720;
        matrix[1][4 + i] = h * f(x[4 + i], matrix[0][4 + i]);
        matrix[2][4 + i - 1] = matrix[1][4 + i] - matrix[1][4 + i - 1];
        matrix[3][4 + i - 2] = matrix[2][4 + i - 1] - matrix[2][4 + i - 2];
        matrix[4][4 + i - 3] = matrix[3][4 + i - 2] - matrix[3][4 + i - 3];
        matrix[5][4 + i - 4] = matrix[4][4 + i - 3] - matrix[4][4 + i - 4];
    }

    if (N - 2 >= 0) std::copy(matrix[0] + 5, matrix[0] + N +3 , adams_values);

    for (int i = 0; i < 200; i++)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}

void runge_kutta(double y[], double x[], int N, double h, double y0)
{
    double k_1, k_2, k_3, k_4;
    y[0] = y0;
    for (int i = 1; i <= N; i++)
    {
        k_1 = h * f(x[i + 1], y[i - 1]);
        k_2 = h * f(x[i + 1] + 0.5 * h, y[i - 1] + 0.5 * k_1);
        k_3 = h * f(x[i + 1] + 0.5 * h, y[i - 1] + 0.5 * k_2);
        k_4 = h * f(x[i + 1] + h, y[i - 1] + k_3);
        y[i] = y[i - 1] + (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6;
    }
}

void euler_I(double y[], double x[], int N, double h, double y0)
{
    y[0] = y0;
    for (int i = 1; i <= N; i++)
    {
        y[i] = y[i - 1] + h * f(y[i + 1] + 0.5 * h, y[i - 1] + 0.5 * h * f(x[i + 1], y[i - 1]));
    }
}

void euler_II(double y[], double x[], int N, double h, double y0)
{
    y[0] = y0;
    for (int i = 1; i <= N; i++)
    {
        y[i] = y[i - 1] + 0.5 * h * (f(y[i + 1], y[i - 1]) + f(x[i + 2], y[i - 1] + h * f(x[i + 1], y[i - 1])));
    }
}

#endif // FUNCTIONS_H_INCLUDED
