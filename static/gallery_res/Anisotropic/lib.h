#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct Anisotropic {
    Eigen::Matrix<double, 3, 3> A;
    double I_5;
    Eigen::MatrixXd frac_partial_differential_²I₅_partial_differential_f²;

    Anisotropic(
        const Eigen::Matrix<double, 3, 1> & a,
        const Eigen::Matrix<double, 3, 3> & C)
    {
        // A = a a^T 
        A = a * a.transpose();
        // `$I_5$` = tr(CA)
        I_5 = (C * A).trace();
        Eigen::MatrixXd frac_partial_differential_²I₅_partial_differential_f²_0(9, 9);
        frac_partial_differential_²I₅_partial_differential_f²_0 << A(1-1, 1-1) * Eigen::MatrixXd::Identity(3, 3), A(1-1, 2-1) * Eigen::MatrixXd::Identity(3, 3), A(1-1, 3-1) * Eigen::MatrixXd::Identity(3, 3),
        A(2-1, 1-1) * Eigen::MatrixXd::Identity(3, 3), A(2-1, 2-1) * Eigen::MatrixXd::Identity(3, 3), A(2-1, 3-1) * Eigen::MatrixXd::Identity(3, 3),
        A(3-1, 1-1) * Eigen::MatrixXd::Identity(3, 3), A(3-1, 2-1) * Eigen::MatrixXd::Identity(3, 3), A(3-1, 3-1) * Eigen::MatrixXd::Identity(3, 3);
        // `$\frac{∂²I₅}{∂f²}$` = 2[A₁,₁I₃  A₁,₂I₃  A₁,₃I₃
    //                A₂,₁I₃  A₂,₂I₃  A₂,₃I₃
    //                A₃,₁I₃  A₃,₂I₃  A₃,₃I₃] 
        frac_partial_differential_²I₅_partial_differential_f² = 2 * frac_partial_differential_²I₅_partial_differential_f²_0;
    }
};

