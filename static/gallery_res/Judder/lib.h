#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct judder {
    double J;
    double F_a;
    double F_b;
    double CFF(
        const double & L)
    {
        return a * log(L) + b;    
    }
    double α(
        const double & F)
    {
        return 1 / double(F);    
    }
    double β(
        const double & L)
    {
        return log10(L);    
    }
    judder(
        const double & F,
        const double & L,
        const double & S,
        const std::function<double(double, double, double)> & P,
        const double & a,
        const double & b,
        const double & M,
        const double & L_a,
        const double & L_b)
    {
        // F_a = M⋅ CFF(L_a)
        F_a = M * CFF(L_a);
        // F_b = M⋅ CFF(L_b)
        F_b = M * CFF(L_b);
        // J = P(α(F), β(L), S)
        J = P(α(F), β(L), S);
    }
};

struct second {
    double E;
    second(
        const std::vector<double> & O,
        const std::vector<double> & M)
    {
        const long dim_0 = O.size();
        assert( M.size() == dim_0 );
        double sum_0 = 0;
        for(int i=1; i<=O.size(); i++){
            sum_0 += abs(log(O.at(i-1)) - log(M.at(i-1))) / double(log(O.at(i-1)));
        }
        // E = sum_i |log(O_i) - log(M_i)|/log(O_i) 
        E = sum_0;
    }
};

struct third {
    double L_b;
    third(
        const double & a,
        const double & b,
        const double & F_b,
        const double & F_a,
        const double & L_a)
    {
        // L_b = 10^((a F_blog((L_a))+b(F_b-F_a))/(aF_a))
        L_b = pow(10, ((a * F_b * log((L_a)) + b * (F_b - F_a)) / double((a * F_a))));
    }
};

