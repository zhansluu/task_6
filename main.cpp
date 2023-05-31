#include <iostream>
#include "functions.h"

using namespace std;

int main()
{
    setlocale(0, "Russian");
    cout << "������� 6. ��������� ������� ������ ���� ��� ������������� ����������������� ��������� ������� �������." << endl << endl;
    cout << "����: y'(x)=-2y(x)+(y(x))^2;  y(0)=1;  h=0.1  N=10" << endl;
    cout << "������ ������� ������ ����: 2 / (exp(2 * x) + 1)" << endl;
    /**��������� ������*/
    double x0 = 0;
    double y0 = 1;

    double h, x[1000], y[1000];
    int N;
    cout << "������� ��������� ������: ��� h = ";
    cin >> h;
    while (h <= 0)
    {
        cout << endl << "������� ������������ �������� h! ��������� �������.";
        cin >> h;
    }
    cout << "N = ";
    cin >> N;
    while (N <= 0)
    {
        cout << endl << "������� ������������ �������� N! ��������� �������.";
        cin >> N;
    }

    cout << "������ ��������: " << endl;
    for (int k = -2; k <= N; k++)
    {
        x[k+2] = x0 + k * h;
        y[k+2] = exact_solution(x[k+2]);
        cout << "x[" << k << "] = " << x[k+2] << "; y[" << k << "] = " << y[k+2] << ";" << endl;
    }
    cout << "------------------------------------------------" << endl;

    cout << "����� �������" << endl;
    double val_taylor[1000]; //������������ �������� � ������ x_k, k=-2..N, ���������� ������� ���������� � ��� �������
    taylor(val_taylor, x, N, y0);
    for (int k = -2; k <= N; k++)
    {
        cout << "x[" << k << "] = " << x[k+2] << "; y[" << k << "] = " << val_taylor[k+2] << ";" << endl;
        cout << "���������� �����������: " << fabs(val_taylor[k+2] - y[k+2]) << endl;
    }
    cout << "------------------------------------------------" << endl;

    cout << "����������������� ����� ������" << endl;
    double val_adams[1000];
    adams(val_adams, val_taylor, x, N, h);
    for (int k = 0; k < N-2; k++)
    {
        cout << "x[" << k+3 << "] = " << x[k+5] << "; y[" << k+3 << "] = " << val_adams[k] << ";" << endl;
    }
    cout << "���������� �����������: " << fabs(val_adams[N-3] - y[N]) << endl;
    cout << "------------------------------------------------" << endl;

    cout << "����� �����-�����" << endl;
    double val_runge_kutta[1000];
    runge_kutta(val_runge_kutta, x, N, h, y0);
    for (int k = 1; k <= N; k++)
    {
        cout << "x[" << k << "] = " << x[k+2] << "; y[" << k << "] = " << val_runge_kutta[k] << ";" << endl;
    }
    cout << "���������� �����������: " << fabs(val_runge_kutta[N] - y[N]) << endl;
    cout << "------------------------------------------------" << endl;

    cout << "����� ������ I" << endl;
    double val_euler_I[1000];
    euler_I(val_euler_I, x, N, h, y0);
    for (int k = 1; k <= N; k++)
    {
        cout << "x[" << k << "] = " << x[k+2] << "; y[" << k << "] = " << val_euler_I[k] << ";" << endl;
    }
    cout << "���������� �����������: " << fabs(val_euler_I[N] - y[N]) << endl;
    cout << "------------------------------------------------" << endl;

    cout << "����� ������ I�" << endl;
    double val_euler_II[1000];
    euler_II(val_euler_II, x, N, h, y0);
    for (int k = 1; k <= N; k++)
    {
        cout << "x[" << k << "] = " << x[k+2] << "; y[" << k << "] = " << val_euler_II[k] << ";" << endl;
    }
    cout << "���������� �����������: " << fabs(val_euler_II[N] - y[N]) << endl;
    cout << "------------------------------------------------" << endl;

    return 0;
}
