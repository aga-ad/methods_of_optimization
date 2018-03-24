#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double sqr(double x) {
    return x*x;
}

// The first variant
// 6*x^2 + 7*x^6 - x*y^2 + 2*x^3*y^3 + 8*y^4 + 5*y^6
double f(const double arg[2]) {
    double x = arg[0];
    double y = arg[1];
    return 6 * pow(x, 2) + 7 * pow(x, 6) - x * pow(y, 2) + 2 * pow(x, 3) * pow(y, 3) + 8 * pow(y, 4) + 5 * pow(y, 6);
}

/*
Вычисление градиента численным методом
point - точка
f - Функция
d - результирующий градиент
*/
template<unsigned int NARG>
void derivative(const double point[NARG], double (*f)(const double arg[NARG]), double d[NARG]) {
    const double eps = 0.0000001;
    double point1[NARG];
    for (int i = 0; i < NARG; i++) {
        point1[i] = point[i];
    }
    for (int i = 0; i< NARG; i++) {
        point1[i] += eps;
        d[i] = (f(point1) - f(point)) / eps;
        point1[i] -= eps;
    }
}

/*
Квадрат нормы вектора
vec - вектор
*/
template<unsigned int NARG>
double norm2(const double vec[NARG]) {
    double res = 0;
    for (int i = 0; i < NARG; i++) {
        res += sqr(vec[i]);
    }
    return res;
}

/*
Градиентный спуск
point - начальная точка
f - Функция
lambda - параметр начального шага
lambda_mult - изменение параметра шага
eps2 - квадрат нормы производной для остановки
max-steps - максимальное количество шагов
*/
template<unsigned int NARG>
void gradient_descent(double point[NARG], double (*f)(const double arg[NARG]), double lambda, double lambda_mult, double eps2, size_t max_steps) {
    double d[NARG];
    double normal;
    for (size_t i = 0; i < max_steps; i++) {
        derivative<NARG>(point, f, d);
        normal = norm2<NARG>(d);
        //cout << i << " " << normal << endl;
        if (normal < eps2) {
            return;
        }
        for (int j = 0; j < NARG; j++) {
            point[j] = point[j] * d[j] * lambda / normal;
        }
        lambda *= lambda_mult;
    }

}

/*
Покоординатный поиск, в качестве метода для оптимизации функции от одной переменной используется примитивный перебор по сетке от l до r с количеством узлов n
point - начальная точка
f -функция
l - левая граница сетки
r - права граница сетки
n - количество узлов сетки
eps - разница в функциях для остановки
max_steps - максимальное количество проходов по всем переменным
*/
template<unsigned int NARG>
void coordinate_search(double point[NARG], double (*f)(const double arg[NARG]), double l, double r, int n, double eps, size_t max_steps) {
    double last = f(point);
    double d, mnx, mnf;
    for (size_t i = 0; i < max_steps; i++) {
        for (int j = 0; j < NARG; j++) {
            d = (r - l) / n;
            mnx = r;
            point[j] = mnx;
            mnf = f(point);
            for (double x = l; x < r; x += d) {
                point[j] = x;
                double t = f(point);
                if (t < mnf) {
                    mnf = t;
                    mnx = x;
                }
            }
            point[j] = mnx;
        }
        if (last - mnf < eps) {
            return;
        }
        last = mnf;
    }
}

int main() {
    cout << fixed << setprecision(8);

    {
        double point[2] = {-2.2, 3.5};
        gradient_descent<2>(point, f, 1, 0.9, 0.000000000001, 10000);
        cout << "point = " << point[0] << " " << point[1] << "    local min = " << f(point) << endl;
    }
    {
        double point[2] = {-2.2, 3.5};
        coordinate_search<2>(point, f, -20, 20, 200, 0.000001, 10000);
        cout << "point = " << point[0] << " " << point[1] << "    local min = " << f(point) << endl;
    }
    return 0;
}
