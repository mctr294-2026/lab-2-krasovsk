#include <cmath>
#include <functional>

#include "roots.hpp"

bool bisection(std::function<double(double)> f, double a, double b, double *root){
    if (f(a)*f(b)<0){
        double c;
        c = (a+b)/2;

        while (std::abs(f(c)) > 1e-6){
            if (f(a)*f(c) < 0){
                b = c;
                c = (a+b)/2;
            } else {
                a = c;
                c = (a+b)/2;
            }
        }
        *root = c;
        return true;
    } else return false;
}

bool regula_falsi(std::function<double(double)> f, double a, double b, double *root){
    if (f(a)*f(b) < 0){
        double c;
        c = a-(f(a)*(b-a))/(f(b)-f(a));
        while (std::abs(f(c)) > 1e-6){
            if (f(a)*f(c) < 0){
                b = c;
            } else if (f(b)*f(c) < 0) {
                a = c;
            }
            c = a-(f(a)*(b-a))/(f(b)-f(a));
        }
        *root = c;
        return true;
    } else return false;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double* root){

    for (int i = 0; i < 1000000; ++i) {
        if (std::abs(g(c)) < 1e-15) return false;

        double d = c - f(c) / g(c);

        if (d < a || d > b) return false;

        if (std::abs(d - c) <= 1e-6) {
            *root = d;
            return true;
        }
        c = d;
    }

    *root = c;
    return true;
}

bool secant(std::function<double(double)> f, double a, double b, double c, double *root)
{
    double denom;
    denom = (f(b) - f(c));
    
    double d = c - f(c) * ((b-c)/denom);

    while (std::abs(d-c) > 1e-6) {
        b = c;
        c = d;
        denom = (f(b) - f(c));
        d = c - f(c) * ((b-c) / denom);
    }
    *root = c;
    return true;
}