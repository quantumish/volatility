#include <functional>
#include <iostream>
#include <fstream>
#include <random>
#include <cassert>
#include <cmath>
#include "date.h"

struct FuncPair {
    std::function<double(double)> f;
    std::function<double(double)> f_prime;

    FuncPair(std::function<double(double)> function, std::function<double(double)> derivative);
    FuncPair(std::function<double(double)> function, double eps);
};

FuncPair::FuncPair(std::function<double(double)> function, std::function<double(double)> derivative)
    :f(function), f_prime(derivative)
{}

FuncPair::FuncPair(std::function<double(double)> function, double eps)
    :f(function)
{
    f_prime = [this, eps](double x) -> double {
        return (f(x+eps) - f(x))/eps;
    };
}

double newton_raphson(double initial, double threshold, const FuncPair& func)
{
    double y = func.f(initial);
    double x = initial;
    while (fabs(y) > threshold) {
        x = x - func.f(x)/func.f_prime(x);
        y = func.f(x);
    }
    return x;
}

struct constants {
    double C;
    double K;
    double r;
    double t;
    double S;
};

#define EPSILON 0.1
#define LARGE_NUM 10000
double black_scholes (double x, constants c)
{
    auto N = [](double x) -> double {
        return 0.5 * erfc(-x * M_SQRT1_2);
    };
    double d_1 = (log(c.S/c.K) + (c.r + pow(x, 2)/2)*c.t)/sqrt(pow(x,2) * c.t);
    double d_2 = (log(c.S/c.K) + (c.r - pow(x, 2)/2)*c.t)/sqrt(pow(x,2) * c.t);
    return c.S*N(d_1) - c.K*exp(-c.r*c.t)*N(d_2) - c.C;
}

int rdn(int y, int m, int d) { /* Rata Die day one is 0001-01-01 */
    if (m < 3) y--, m += 12;
    return 365*y + y/4 - y/100 + y/400 + (153*m - 457)/5 + d - 306;
}

using namespace date;
int main(int argc, char** argv)
{
    std::ifstream file (argv[1]);
    if (!file) {
        std::cout << "Bad file." << "\n";
        return 1;
    }
    std::string line;
    constants c = {0,0,0,0,0};
    while(getline(file, line)) {
        std::tm e;
        int m_0 = 3, d_0 = 3, y_0 = 2020;
        int m,d,y;
        sscanf(line.c_str(), "%lf %0d/%0d/%0d %lf", &(c.C), &d, &m, &y, &(c.K));
        y+=2000;
        //        std::cout << o.tm_mday << '/' << o.tm_mon << '/' << o.tm_year << " " << e.tm_mday << '/' << e.tm_mon << '/' << e.tm_year << " ";
        c.S = 58.29;
        c.r = 0.05;
        auto d0 = year{y_0}/m_0/d_0;
        auto d1 = year{y}/m/d;
        c.t = (sys_days{d1} - sys_days{d0}).count();
        auto wrapper = [c] (double x) -> double {
            return black_scholes(x, c);
        };
        FuncPair option (wrapper, 0.001);
        double x = newton_raphson(strtod(argv[2], NULL), 0.001, option);
        std::cout << c.S << " " << c.K << " " << c.r << " " << c.t << " " << c.C << "\t\t";
        std::cout << x << "\n";
    }
}
