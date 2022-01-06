#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct first {
    double ϕ;
    double D;
    double γ;
    Eigen::Matrix<double, 3, 3> K;
    Eigen::Matrix<double, 3, 1> ŧ;
    double A;
    double B;
    double C;
    double α_1;
    double α_2;
    double ω(
        const double & x)
    {
        return atan((x - c_x) / double(f_x));    
    }
    double s(
        const double & x)
    {
        return (y - c_y) * cos(ω(x));    
    }
    Eigen::Matrix<double, 3, 3> R(
        const double & α)
    {
        Eigen::Matrix<double, 3, 3> R_0;
        R_0 << -sin(α), 0, -cos(α),
        0, 1, 0,
        cos(α), 0, -sin(α);
        return R_0;    
    }
    Eigen::MatrixXd P(
        const double & α)
    {
        Eigen::MatrixXd P_0(3, 4);
        P_0 << R(α), ŧ;
        return K * P_0;    
    }
    double t(
        const double & α)
    {
        return (α - α_i) / double((α_j - α_i));    
    }
    Eigen::Matrix<double, 2, 1> x(
        const double & α)
    {
        return φ_circumflex_accent_left_curly_bracket_hyphen_minus_1_right_curly_bracket_d((1 - t(α)) * φ_d(x_i) + t(α) * φ_d(x_j));    
    }
    Eigen::VectorXd φ(
        const double & x)
    {
        Eigen::VectorXd φ_0(2);
        φ_0 << ω(x), s(x);
        return φ_0;    
    }
    first(
        const double & f_x,
        const double & f_y,
        const double & c_x,
        const double & c_y,
        const double & r,
        const std::function<Eigen::Matrix<double, 2, 1>(Eigen::Matrix<double, 2, 1>)> & φ_circumflex_accent_left_curly_bracket_hyphen_minus_1_right_curly_bracket_d,
        const std::function<Eigen::Matrix<double, 2, 1>(Eigen::Matrix<double, 2, 1>)> & φ_d,
        const Eigen::Matrix<double, 2, 1> & x_i,
        const Eigen::Matrix<double, 2, 1> & x_j,
        const double & α_i,
        const double & α_j,
        const double & y,
        const double & X,
        const double & Y,
        const double & Z,
        const double & χ)
    {
        Eigen::Matrix<double, 3, 3> K_0;
        K_0 << f_x, 0, c_x,
        0, f_y, c_y,
        0, 0, 1;
        // K = [f_x 0 c_x
        //       0   f_y c_y
        //       0      0    1]
        K = K_0;
        Eigen::Matrix<double, 3, 1> ŧ_0;
        ŧ_0 << 0,
        0,
        -r;
        // ŧ = [0;0;-r]
        ŧ = ŧ_0;
        // A = X ⋅ f_x - Z⋅(χ - c_x )
        A = X * f_x - Z * (χ - c_x);
        // B = Z⋅f_x + X⋅(χ -c_x ) 
        B = Z * f_x + X * (χ - c_x);
        // D = √(A^2 +B^2)
        D = sqrt((pow(A, 2) + pow(B, 2)));
        // γ = tan(B/A)
        γ = tan(B / double(A));
        // C = -r⋅(χ -c_x )
        C = -r * (χ - c_x);
        // ϕ = sin(C/D) 
        ϕ = sin(C / double(D));
        // α_1 = ϕ - γ
        α_1 = ϕ - γ;
        // α_2 = π - ϕ - γ  
        α_2 = M_PI - ϕ - γ;
    }
};

