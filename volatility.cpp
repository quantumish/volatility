#include <functional>
#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
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

int main(int argc, char** argv)
{
    std::ifstream file (argv[1]);
    if (!file) {
        std::cout << "Bad file." << "\n";
        return 1;
    }
    std::string line;
    constants c = {0,0,0,0,0};
    int count = 0;
    double sum = 0;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << std::left << std::setw(10) << "S" << std::left << "\t"
              << std::setw(10) << "K" << std::left << "\t"
              << std::setw(10) << "r" << std::left << "\t"
              << std::setw(10) << "t" << std::left << "\t"
              << std::setw(10) << "C" << "\t\t"
              << "σ" << "\n";
    while(getline(file, line)) {
        sscanf(line.c_str(), "%lf %lf %lf", &(c.C), &(c.t), &(c.K));
        c.S = 58.29;
        c.r = 0.05;
        auto wrapper = [c] (double x) -> double {
            return black_scholes(x, c);
        };
        FuncPair option (wrapper, 0.001);
        double x = newton_raphson(strtod(argv[2], NULL), 0.001, option);
        std::cout << std::setw(10) << c.S << std::left << "\t"
                  << std::setw(10) << c.K << std::left << "\t"
                  << std::setw(10) << c.r << std::left << "\t"
                  << std::setw(10) << c.t << std::left << "\t"
                  << std::setw(10) << c.C << "\t\t";
        std::cout << x << "\n";
        sum+=x;
        count++;
    }
    std::cout << "Mean σ: " << sum/count << "\n";
}
