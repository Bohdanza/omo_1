#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <utility>

using namespace std;

const double EPS = 1e-5;
const double F_EPS = 1e-6;

double f(double x)
{
    return (x - 3) * exp(x - 4) - 2 * cos(2 * x) - 3 * x + 8;
}
double df(double x)
{
    return (x - 2) * exp(x - 4) + 4 * sin(2 * x) - 3;
}
double d2f(double x)
{
    return (x - 1) * exp(x - 4) + 8 * cos(2 * x);
}

double phi(double x)
{
    return 4 + log((2 * cos(2 * x) + 3 * x - 8) / (x - 3));
}

void header()
{
    cout << setw(4) << "k" << "  " << setw(18) << "x_k" << "  " << setw(18) << "x_k - x_(k-1)" << "  " << setw(18)
         << "f(x_k)" << '\n';
}

void row(int k, double x, double dx, double fx)
{
    cout << setw(4) << k << "  " << showpos << fixed << setprecision(10) << setw(18) << x << "  " << scientific << setw(18);
    if (k == 0)
        cout << "-";
    else
        cout << dx;
    cout << "  " << setw(18) << fx << noshowpos << '\n';
}

void iterate(const char *name, double x0, const function<double(double)> &step)
{
    cout << name << '\n';
    header();
    row(0, x0, 0, f(x0));
    double x = x0;
    for (int k = 1;; ++k)
    {
        double xn = step(x), fxn = f(xn);
        row(k, xn, xn - x, fxn);
        if (!isfinite(xn) || (fabs(xn - x) < EPS && fabs(fxn) < F_EPS))
            break;
        x = xn;
    }
    cout << '\n';
}

pair<double, double> sampling()
{
    cout << "Initial approximation:\n";
    cout << setw(8) << "x" << "  " << setw(12) << "f(x)" << '\n';
    pair<double, double> br{NAN, NAN};
    double prevX = 0, prevF = 0;
    for (int i = 0; i <= 40; ++i)
    {
        double x = -10 + 0.5 * i, fx = f(x);
        cout << showpos << fixed << setprecision(2) << setw(8) << x << "  " << setprecision(6) << setw(12) << fx << noshowpos
             << '\n';
        if (i > 0 && prevF * fx < 0)
            br = {prevX, x};
        prevX = x;
        prevF = fx;
    }
    cout << '\n';
    return br;
}

double fourierPoint(pair<double, double> br)
{
    return f(br.first) * d2f(br.first) > 0 ? br.first : br.second;
}

int main()
{
    pair<double, double> br = sampling();

    iterate("Simple iteration method", (br.first + br.second) / 2, phi);
    
    double x0 = fourierPoint(br), df0 = df(x0);
    iterate("Modified Newton method", x0, [df0](double x) { return x - f(x) / df0; });
    
    iterate("Steffensen method", x0, [](double x) 
    {
        double fx = f(x);
        return x - fx * fx / (f(x + fx) - fx);
    });

    return 0;
}
