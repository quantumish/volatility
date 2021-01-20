#include <functional>
#include <iostream>
#include <random>
#include <cassert>
#include <cmath>


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
        std::cout << x << " " << y << '\n';
    }
    return x;
}


#define EPSILON 0.1
#define LARGE_NUM 10000
double black_scholes (double x)
{
    double C = 3.7888;
    double K = 50;
    double r = 0.05;
    double t = 0.5;
    double S = 52;
    auto N = [](double x) -> double {
        return 0.5 * erfc(-x * M_SQRT1_2);
    };
    double d_1 = (log(S/K) +  (r + pow(x, 2)/2)*t)/sqrt(pow(x,2) * t);
    double d_2 = (log(S/K) + (r - pow(x, 2)/2)*t)/sqrt(pow(x,2) * t);
    return S*N(d_1) - K*exp(-r*t)*N(d_2) - C;
}

int main(int argc, char** argv) {
    FuncPair option (black_scholes, 0.001);
    std::cout << black_scholes(strtod(argv[1], NULL)) << "\n";
    double x = newton_raphson(0.18, 0.01, option);
    std::cout << x << "\n";
}
