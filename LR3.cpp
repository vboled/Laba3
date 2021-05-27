#include "iostream"

#include <fstream>
#include <typeinfo>
#include <cmath>
#include <malloc.h>


using namespace std;


const double PI = 3.141592653589793;
const double Eps = 1e-5;
const double EPSILON = 1e-10;

const double L = 10.;
const double T = 1;
const double a = 1.;

double f(double x, char s)
{
    switch(s)
    {
        case '1': return sin(PI * x);
        case '2': return x * (1 - x);
        case '3': return (x + 0.5) * (x + 1);
        case '4': {
            if (x >= 10 / 3 && x <= 20 / 3)
                return 1;
            return 0;
        }
        case '5': {
            if (x >= 10 / 3 && x <= 20 / 3)
                return 1;
            return 0;
        }
    }
    return -1;
}

double g(double x, char s)
{
    switch(s)
    {
        case '1': return 0;
        case '2': return 0;
        case '3': return cos(x + 0.5);
        case '4': return 0;
        case '5': return 0;
    }
    return -1;
}

double psix(double t)
{
    return 0;
}

double phi(double t, char s)
{
    switch(s)
    {
        case '1': return 0.0;
        case '2': return 0.0;
        case '3': return 0.5;
        case '4': return 0;
        case '5': return 0;
    }
    return -1;
}

double psi(double t, char s)
{
    switch(s)
    {
        case '1': return 0.0;
        case '2': return 0.0;
        case '3': return 3 - 2 * t;
        case '4': return 0;
        case '5': return 0;
    }
    return -1;
}

double ddf(double x, char s, double h)
{
    switch(s)
    {
        case '1': return -pow(PI, 2) * sin(PI * x);
        case '2': return -2;
        case '3': return 2;
        case '4': return 0;
        case '5': return (f(x + h, s) - 2 * f(x, s) + f(x - h, s)) / h / h;
    }
    return -1;
}

void Krest(int n, int m,char s)
{
    double h = L / (n - 1);
    double tau = T / (m - 1);


    double *Y0 = new double [n];
    double *Y1 = new double [n];
    double *Y2 = new double [n];
    double *Y_spec;

    for (int i = 1; i < n - 1; i++)
        Y0[i] = f(i * h, s);
    Y0[0] = phi(0, s);
    Y0[n - 1] = psi(0, s);

    ofstream f;
    string filename = "LR3_";
    filename.push_back(s);
    filename += ".txt";
    f.open(filename);
    f << L << " ";
    for (size_t i = 0; i < n; i++)
        f << T << " ";
    f << endl;
    f << n << " ";
    for (size_t i = 0; i < n; i++)
        f << m << " ";
    f << endl;
    for (int i = 1; i < n - 1; i++) {
        Y1[i] = Y0[i] + tau * g(i * h, s) + pow(a * tau, 2) / 2 * ddf(i * h, s, h);
    }
    Y1[0] = phi(tau, s);
    Y1[n - 1] = psi(tau, s);

    f << 0 << " ";
    for (size_t i = 0; i < n; i++)
        f << Y0[i] << " ";
    f << endl;
    f << tau << " ";
    for (size_t i = 0; i < n; i++)
        f << Y1[i] << " ";
    f << endl;

    for (int j = 2; j < m; j++)
    {
        for (int i = 1; i < n - 1; i++)
        {
            Y2[i] = pow(a * tau / h, 2) * (Y1[i + 1] + Y1[i - 1]) + 2 * (1 - pow( a * tau / h, 2)) * Y1[i] - Y0[i];
        }
        Y2[0] = phi(j * tau, s);
        Y2[n - 1] = psi(j * tau, s);

        f << tau * j << " ";
        for (size_t i = 0; i < n; i++)
            f << Y2[i] <<" ";
        f << endl;

        Y_spec = Y0;
        Y0 = Y1;
        Y1 = Y2;
        Y2 = Y_spec;
    }

    f.close();
    f.clear();
    delete [] Y0;
    delete [] Y1;
    delete [] Y2;
    delete [] Y_spec;
}

int main()
{
    setlocale(LC_ALL, "Russian");

    cout << "Task:" << endl;
    cout << "1" << endl;
    cout << "2" << endl;
    cout << "3" << endl;
    char s;
    cin >> s;
    int n;
    cout << "steps L:" << endl;
    cin >> n;
    int m;
    cout << "steps T:" << endl;
    cin >> m;
    Krest(n,m,s);
}

