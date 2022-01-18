#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct kappa {
    double x;
    double y;
    double r;
    double θ;
    double a;
    double κ(
        const double & θ)
    {
        return (8 * (3 - pow(sin(θ), 2)) * pow(sin(θ), 4)) / double((a * pow((pow(sin(2 * θ), 2) + 4), (3 / double(2)))));    
    }
    double phi(
        const double & θ)
    {
        return -atan(1 / double(2) * sin(2 * θ));    
    }
    kappa(
        const double & a,
        const double & t,
        const double & θ)
    {
        this->θ = θ;
        this->a = a;
        // x = a cos(t) cot(t)
        x = a * cos(t) * 1/tan(t);
        // y = a cos(t)
        y = a * cos(t);
        // r = atan(θ)
        r = atan(θ);
    }
};

