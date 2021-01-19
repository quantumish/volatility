#include <functional>
#include <iostream>
#include <cmath>

template <typename T>
double newton_raphson(double initial, double threshold, const T& func)
{
    double y = func.f(initial);
    double x = initial;
    while (fabs(y) > threshold) {       
        x = x - func.f(x)/func.f_prime(x);
        y = func.f(x);
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

double quad (double x) {return pow(x,2)+(3*x)+2.25;}

int main() {
    FuncPair quadratic (quad, 0.001);
    double x = newton_raphson<FuncPair>(0, 0.0001, quadratic);
    std::cout << x << "\n";
}
