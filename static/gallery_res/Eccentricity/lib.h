#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct first {

    Eigen::MatrixXd 𝑔(
        const Eigen::VectorXd & x,
        const Eigen::VectorXd & `$x_0$`,
        const double & 𝜃,
        const double & 𝜎,
        const double & `$𝑓_𝑠$`)
    {
        const long n = x.size();
        assert( `$x_0$`.size() == n );

        Eigen::Matrix<double, 1, 2> 𝑔_0;
        𝑔_0 << cos(𝜃), sin(𝜃);
        return exp(-pow((x - `$x_0$`).lpNorm<2>(), 2) / double((2 * pow(𝜎, 2)))) * 2 * M_PI * `$𝑓_𝑠$` * x * 𝑔_0.unaryExpr<double(*)(double)>(&std::cos);    
    }
    double 𝜏(
        const double & `$𝑓_𝑠$`)
    {
    
    

    
        return m(log(`$𝑓_𝑠$`) - log(𝑓_left_curly_bracket_𝑠0_right_curly_bracket), 0);    
    }
    double 𝜁(
        const double & `$𝑓_𝑠$`)
    {
    
    

    
        return exp(𝑝[9-1] * 𝜏(`$𝑓_𝑠$`)) - 1;    
    }
    double Ψ(
        const double & 𝑒,
        const double & `$𝑓_𝑠$`)
    {
    
    

    
        return m(0, 𝑝[0-1] * 𝜏(`$𝑓_𝑠$`) + 𝑝[1-1] * 𝜏(`$𝑓_𝑠$`) + 𝑝[2-1] + (𝑝[3-1] * pow(𝜏(`$𝑓_𝑠$`), 2) + 𝑝[4-1] * 𝜏(`$𝑓_𝑠$`) + 𝑝[5-1]) * 𝜁(`$𝑓_𝑠$`) * 𝑒 + (𝑝[6-1] * pow(𝜏(`$𝑓_𝑠$`), 2) + 𝑝[7-1] * 𝜏(`$𝑓_𝑠$`) + 𝑝[8-1]) * 𝜁(`$𝑓_𝑠$`) * pow(𝑒, 2));    
    }
    double A(
        const double & 𝑒)
    {
    
    

    
        return log(64) * 2.3 / double((0.106 * (𝑒 + 2.3)));    
    }
    double d(
        const double & L)
    {
    
    

    
        return 7.75 - 5.75 * (pow((L * a / double(846)), 0.41) / double((pow((L * a / double(846)), 0.41) + 2)));    
    }
    double 𝑠(
        const double & 𝑒,
        const double & `$𝑓_𝑠$`)
    {
    
    

    
        return 𝜁(`$𝑓_𝑠$`) * (q[0-1] * pow(𝑒, 2) + q[1-1] * 𝑒) + q[2-1];    
    }
    first(
        const std::function<double(double, double)> & m,
        const Eigen::Matrix<double, 10, 1> & 𝑝,
        const double & 𝑓_left_curly_bracket_𝑠0_right_curly_bracket,
        const double & a,
        const Eigen::Matrix<double, 3, 1> & q)
    {
    
    }
};

