#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

/*
Функция, первый варинт.
*/
double func(double x) {
    return exp(sqrt(x)) + 2.1 * exp(-2.1 * x);
}

/*
Метод сканирования.
f - функция
l - левая границы
r - правая граница
n - количество узлов сетки, не считая узел на правой границе отрезка
*/
double stupid_search(double (*f)(double), double l, double r, unsigned long long n) {
    double d = (r - l) / n;
    double mnx = r;
    double mnf = f(mnx);
    for (double x = l; x < r; x += d) {
        double t = f(x);
        if (t < mnf) {
            mnf = t;
            mnx = x;
        }
    }
    return mnx;
}

/*
Метод дихотомии.
f - функция
l - левая границы
r - правая граница
n - количество итераций
*/
double dichotomy_method(double (*f)(double), double l, double r, unsigned int n) {
    for (unsigned long long i = 0; i < n; i++) {
        double mid = (r + l) / 2;
        double d = (r - l) / 30;
        double a = mid - d;
        double b = mid + d;
        if (f(a) < f(b)) {
            r = b;
        } else {
            l = a;
        }
    }
    return (l + r) / 2;
}

/*
Метод золотого сечения.
f - функция
l - левая границы
r - правая граница
n - количество итераций
*/
double golden_section_search(double (*f)(double), double l, double r, unsigned int n) {
    const double phi = 0.6180339887; // (sqrt(5) - 1) / 2
    double x1 = r - (r - l) * phi;
    double x2 = l + (r - l) * phi;
    double f1 = f(x1);
    double f2 = f(x2);
    for (unsigned long long i = 1; i < n; i++) {
        if (f1 < f2) {
            r = x2;
            x2 = x1;
            x1 = r - (r - l) * phi;
            f2 = f1;
            f1 = f(x1);
        } else {
            l = x1;
            x1 = x2;
            x2 = l + (r - l) * phi;
            f1 = f2;
            f2 = f(x2);
        }
    }
    if (f1 < f2) {
        r = x2;
    } else {
        l = x1;
    }
    return (l + r) / 2;
}

/*
Метод чисел Фибоначчи.
f - функция
l - левая границы
r - правая граница
n - количество итераций
*/
double fibonacci_search(double (*f)(double), double l, double r, unsigned int n) {
    double* fib = new double[n + 2];
    fib[0] = 1;
    fib[1] = 2;
    for (unsigned int i = 2; i - 2 < n; i++) {
        fib[i] = fib[i - 2] + fib[i - 1];
    }
    double x1 = l + (r - l) * fib[n - 1] / fib[n + 1];
    double x2 = l + (r - l) * fib[n] / fib[n + 1];
    double f1 = f(x1);
    double f2 = f(x2);
    for (unsigned int i = n; i > 1; i--) {
        if (f1 < f2) {
            r = x2;
            x2 = x1;
            x1 = l + (r - l) * fib[i - 2] / fib[i];
            f2 = f1;
            f1 = f(x1);
        } else {
            l = x1;
            x1 = x2;
            x2 = l + (r - l) * fib[i - 1] / fib[i];
            f1 = f2;
            f2 = f(x2);

        }
    }
    delete[] fib;
    if (f1 < f2) {
        r = x2;
    } else {
        l = x1;
    }
    return (l + r) / 2;
}


int main() {
    cout << fixed << setprecision(8);
    double ss = stupid_search(func, -10, 10, 10000000);
    cout << "stupid search: min: " << func(ss) << "  arg: " << ss << endl;
    double dm = dichotomy_method(func, -10, 10, 50);
    cout << "dichotomy method: min: " << func(dm) << "  arg: " << dm << endl;
    double gss = golden_section_search(func, -10, 10, 50);
    cout << "golden section search: min: " << func(gss) << "  arg: " << gss << endl;
    double fs = fibonacci_search(func, -10, 10, 50);
    cout << "fibonacci search: min: " << func(fs) << "  arg: " << fs << endl;
    return 0;
}
