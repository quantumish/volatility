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

#define EPSILON 0.1
#define LARGE_NUM 10000
double black_scholes (double x)
{
    double C = 4.5;
    double K = 50;
    double r = 0.05;
    double t = 0.5;
    double S = 52;
    auto N = [](double x) -> double {
        double sum = 0;
        for (double i = -LARGE_NUM; i <= x; i+=EPSILON) {           
            sum += (1/(sqrt(2*M_PI)) * exp(-0.5 * pow(i, 2))) * EPSILON;
        }
        return sum;
    };
    std::cout << x << "\n";
    double d_1 = (log(S/K) +  (r + pow(x, 2)/2)*t)/sqrt(pow(x,2) * t);
    double d_2 = (log(S/K) + (r - pow(x, 2)/2)*t)/sqrt(pow(x,2) * t);
    std::cout << d_1 << " " << d_2 << " " << N(d_1) << " " << N(d_2) << " " << S*N(d_1) << " " << K*exp(-r*t)*N(d_2) << " " << S << " "  << K*exp(-r*t) << "\n";
    return S*N(d_1) - K*exp(-r*t)*N(d_2);
}

int main(int argc, char** argv) {
    FuncPair option (black_scholes, 0.001);
    std::cout << black_scholes(strtod(argv[1], NULL)) << "\n";
    double x = newton_raphson<FuncPair>(0, 0.01, option);
    std::cout << x << "\n";
}
