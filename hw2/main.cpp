#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double sqr(double x) {
    return x*x;
}

// f(x, y) = 3 * x^2 + 8 * x^4 - 12 * x * y + 4 * y^2 - 4 * x^2 * y^3 + y^6
// global minimum is -2.1817 in (0.840911, 1.1511)
double f(const double arg[2]) {
    double x = arg[0];
    double y = arg[1];
    return 3 * pow(x, 2) + 8 * pow(x, 4) - 12 * x * y + 4 * pow(y, 2) - 4 * pow(x, 2) * pow(y, 3) + pow(y, 6);
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
double norm(const double vec[NARG]) {
    double res = 0;
    for (int i = 0; i < NARG; i++) {
        res += sqr(vec[i]);
    }
    return sqrt(res);
}

/*
Градиентный спуск
point - начальная точка
f - Функция
lambda - параметр начального шага
lambda_mult - изменение параметра шага
eps - норма производной для остановки
max-steps - максимальное количество шагов
*/
template<unsigned int NARG>
void gradient_descent(double point[NARG], double (*f)(const double arg[NARG]), double lambda, double lambda_mult, double eps, size_t max_steps) {
    double d[NARG];
    double normal;
    for (size_t i = 0; i < max_steps; i++) {
        derivative<NARG>(point, f, d);
        normal = norm<NARG>(d);
        if (normal < eps) {
            return;
        }
        for (int j = 0; j < NARG; j++) {
            point[j] = point[j] - d[j] * lambda / normal;
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
        double point[2] = {2.2, 3.5};
        gradient_descent<2>(point, f, 100, 0.999, 0.00000000001, 100000);
        cout << "point = " << point[0] << " " << point[1] << "    local min = " << f(point) << endl;
        //finds point = -0.46639822 -0.55766976    local min = -0.66506061 local minimum
    }

    {
        double point[2] = {2.2, 3.5};
        gradient_descent<2>(point, f, 100, 0.9, 0.00000000001, 100000);
        cout << "point = " << point[0] << " " << point[1] << "    local min = " << f(point) << endl;
        // finds point = 0.84091106 1.15109548    local min = -2.18169756 global minimum
    }

    {
        double point[2] = {-12.2, 13.5};
        coordinate_search<2>(point, f, -20, 20, 10000, 0.000001, 10000);
        cout << "point = " << point[0] << " " << point[1] << "    local min = " << f(point) << endl;
        //finds point = 0.84000000 1.15200000    local min = -2.18162320 global minimum

    }
    return 0;
}
