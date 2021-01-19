#include <functional>
#include <iostream>
#include <random>
#include <cassert>
#include <cmath>

template <typename T>
double newton_raphson(double initial, double threshold, const T& func)
{
    double y = func.f(initial);
    double x = initial;
    while (fabs(y) > threshold) {       
        x = x - func.f(x)/func.f_prime(x);
        y = func.f(x);
        std::cout << x << " " << y << '\n';
        assert(2<1);
    }
    return x;
}

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

double black_scholes (double x)
{
    double C = 17.50;
    double K = 25.00;
    double r = 0.05;
    double t = 3;
    double S = 44.57;
    double rho = 1;
    auto N = [](double x) -> double {
        return 1/(sqrt(2*M_PI)) * exp(-0.5 * pow(x, 2));
    };
    double d_1 = (log(S/K) + (r + pow(rho, 2)/2)*t)/(rho * sqrt(t));
    double d_2 = d_1 - (rho * sqrt(t));
    return S*N(d_1) - K*exp(-r*t)*N(d_2) - C;
}

int main() {
    FuncPair option (black_scholes, 0.001);
    std::cout << black_scholes(0) << "\n";
    // double x = newton_raphson<FuncPair>(0, 0.01, option);
    // std::cout << x << "\n";
}
