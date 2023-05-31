#include <iostream>
#include "functions.h"

using namespace std;

int main()
{
    setlocale(0, "Russian");
    cout << "Задание 6. Численное решение Задачи Коши для обыкновенного дифференциального уравнения первого порядка." << endl << endl;
    cout << "Дано: y'(x)=-2y(x)+(y(x))^2;  y(0)=1;  h=0.1  N=10" << endl;
    cout << "Точное решение задачи Коши: 2 / (exp(2 * x) + 1)" << endl;
    /**Начальные данные*/
    double x0 = 0;
    double y0 = 1;

    double h, x[1000], y[1000];
    int N;
    cout << "Введите параметры задачи: шаг h = ";
    cin >> h;
    while (h <= 0)
    {
        cout << endl << "Введено недопустимое значение h! Повторите попытку.";
        cin >> h;
    }
    cout << "N = ";
    cin >> N;
    while (N <= 0)
    {
        cout << endl << "Введено недопустимое значение N! Повторите попытку.";
        cin >> N;
    }

    cout << "Точные значения: " << endl;
    for (int k = -2; k <= N; k++)
    {
        x[k+2] = x0 + k * h;
        y[k+2] = exact_solution(x[k+2]);
        cout << "x[" << k << "] = " << x[k+2] << "; y[" << k << "] = " << y[k+2] << ";" << endl;
    }
    cout << "------------------------------------------------" << endl;

    cout << "Метод Тейлора" << endl;
    double val_taylor[1000]; //приближенные значения в точках x_k, k=-2..N, полученные методом разложения в ряд Тейлора
    taylor(val_taylor, x, N, y0);
    for (int k = -2; k <= N; k++)
    {
        cout << "x[" << k << "] = " << x[k+2] << "; y[" << k << "] = " << val_taylor[k+2] << ";" << endl;
        cout << "Абсолютная погрешность: " << fabs(val_taylor[k+2] - y[k+2]) << endl;
    }
    cout << "------------------------------------------------" << endl;

    cout << "Экстраполяционный метод Адамса" << endl;
    double val_adams[1000];
    adams(val_adams, val_taylor, x, N, h);
    for (int k = 0; k < N-2; k++)
    {
        cout << "x[" << k+3 << "] = " << x[k+5] << "; y[" << k+3 << "] = " << val_adams[k] << ";" << endl;
    }
    cout << "Абсолютная погрешность: " << fabs(val_adams[N-3] - y[N]) << endl;
    cout << "------------------------------------------------" << endl;

    cout << "Метод Рунге-Кутта" << endl;
    double val_runge_kutta[1000];
    runge_kutta(val_runge_kutta, x, N, h, y0);
    for (int k = 1; k <= N; k++)
    {
        cout << "x[" << k << "] = " << x[k+2] << "; y[" << k << "] = " << val_runge_kutta[k] << ";" << endl;
    }
    cout << "Абсолютная погрешность: " << fabs(val_runge_kutta[N] - y[N]) << endl;
    cout << "------------------------------------------------" << endl;

    cout << "Метод Эйлера I" << endl;
    double val_euler_I[1000];
    euler_I(val_euler_I, x, N, h, y0);
    for (int k = 1; k <= N; k++)
    {
        cout << "x[" << k << "] = " << x[k+2] << "; y[" << k << "] = " << val_euler_I[k] << ";" << endl;
    }
    cout << "Абсолютная погрешность: " << fabs(val_euler_I[N] - y[N]) << endl;
    cout << "------------------------------------------------" << endl;

    cout << "Метод Эйлера IШ" << endl;
    double val_euler_II[1000];
    euler_II(val_euler_II, x, N, h, y0);
    for (int k = 1; k <= N; k++)
    {
        cout << "x[" << k << "] = " << x[k+2] << "; y[" << k << "] = " << val_euler_II[k] << ";" << endl;
    }
    cout << "Абсолютная погрешность: " << fabs(val_euler_II[N] - y[N]) << endl;
    cout << "------------------------------------------------" << endl;

    return 0;
}
