#include <functional>
#include <iostream>
#include <random>
#include <cassert>
#include <cmath>

template <class ...Args>
struct FuncPair {
    std::function<double(double)> f;
    std::function<double(double)> f_prime;

    FuncPair(std::function<double(double)> function, std::function<double(double)> derivative, Args... args);            
    FuncPair(std::function<double(double)> function, double eps, Args... args);
};

FuncPair::FuncPair(std::function<double(double)> function, std::function<double(double)> derivative, Args... args)
    :f(function, args), f_prime(derivative, args)
{}

FuncPair::FuncPair(std::function<double(double)> function, double eps, Args... args)
    :f(function)
{
    f_prime = [this, eps](double x, Args... args) -> double {
        return (f(x+eps, args) - f(x, args))/eps;
    };
}

template <class ...Args>
double newton_raphson(double initial, double threshold, const FuncPair& func, Args... args)
{
    double y = func.f(initial, args);
    double x = initial;
    while (fabs(y) > threshold) {
        x = x - func.f(x, args)/func.f_prime(x, args);
        y = func.f(x, args);
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
double black_scholes (double x)
{
    double C = 1.47;
    double K = 44;
    double r = 0.05;
    double t = (double)3/365;
    double S = 44.75;
    auto N = [](double x) -> double {
        return 0.5 * erfc(-x * M_SQRT1_2);
    };
    double d_1 = (log(S/K) + (r + pow(x, 2)/2)*t)/sqrt(pow(x,2) * t);
    double d_2 = (log(S/K) + (r - pow(x, 2)/2)*t)/sqrt(pow(x,2) * t);
    return S*N(d_1) - K*exp(-r*t)*N(d_2) - C;
}

int main(int argc, char** argv)
{
    std::ifstream file (argv[1]);
    if (!file) {
        std::cout << "Bad file." << "\n";
        return 1;
    }
    std::string line;
    while(getline(file, line)) {
        
        sscanf()
    }    
    FuncPair option (black_scholes, 0.001);
    double x = newton_raphson(0.6, 0.001, option);
    std::cout << x << "\n";
}
