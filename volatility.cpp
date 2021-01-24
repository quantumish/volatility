#include <functional>
#include <iostream>
#include <fstream>
#include <random>
#include <cassert>
#include <cmath>
#include <ctime>

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

int main(int argc, char** argv)
{
    std::ifstream file (argv[1]);
    if (!file) {
        std::cout << "Bad file." << "\n";
        return 1;
    }
    std::string line;
    constants c;
    while(getline(file, line)) {
        std::tm e;
        std::tm o = {0,0,0,3,31,120};
        sscanf(line.c_str(), "%*s Equity,,%f,,UAL US %0d/%0d/%0d C%f",  &c.C, &e.tm_mday, &e.tm_mon, &e.tm_year, &c.K);        
        e.tm_year+=100;
        std::time_t et = std::mktime(&e);
        std::time_t ot = std::mktime(&o);
        c.S = 58.29;
        c.r = 0.05;
        c.t = std::difftime(et, ot) / (60*60*24*365);
        std::cout << c.S << " " << c.K << " " << c.r << " " << c.t << " " << c.C << "\n";
    }    
    auto wrapper = [c] (double x) -> double {
        return black_scholes(x, c);
    };
    FuncPair option (wrapper, 0.001);
    double x = newton_raphson(2, 0.001, option);
    std::cout << x << "\n";
}
