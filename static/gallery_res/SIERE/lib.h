#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct first {
    Eigen::MatrixXd G_left_parenthesis_u_right_parenthesis;
    Eigen::MatrixXd H_left_parenthesis_u_right_parenthesis;
    Eigen::VectorXd v_G;
    Eigen::VectorXd v_H;
    Eigen::VectorXd f_G;
    Eigen::VectorXd f_H;
    Eigen::MatrixXd J_G;
    Eigen::MatrixXd J_H;
    first(
        const Eigen::MatrixXd & U_s,
        const Eigen::MatrixXd & M,
        const Eigen::VectorXd & v,
        const Eigen::VectorXd & f,
        const Eigen::MatrixXd & K)
    {
        const long n = U_s.rows();
        const long s = U_s.cols();
        assert( U_s.rows() == n );
        assert( M.rows() == n );
        assert( M.cols() == n );
        assert( v.size() == n );
        assert( f.size() == n );
        assert( K.rows() == n );
        assert( K.cols() == n );
        // v_G = U_sU_s^T Mv
        v_G = U_s * U_s.transpose() * M * v;
        // v_H =  v - v_G
        v_H = v - v_G;
        // f_G =  MU_sU_s^T f
        f_G = M * U_s * U_s.transpose() * f;
        Eigen::MatrixXd G_left_parenthesis_u_right_parenthesis_0(2*n, 1);
        G_left_parenthesis_u_right_parenthesis_0 << v_G,
        M.colPivHouseholderQr().solve(f_G);
        // G_left_parenthesis_u_right_parenthesis = [v_G
        //           M⁻¹f_G]
        G_left_parenthesis_u_right_parenthesis = G_left_parenthesis_u_right_parenthesis_0;
        // f_H =  f - f_G
        f_H = f - f_G;
        Eigen::MatrixXd H_left_parenthesis_u_right_parenthesis_0(2*n, 1);
        H_left_parenthesis_u_right_parenthesis_0 << v_H,
        M.colPivHouseholderQr().solve(f_H);
        // H_left_parenthesis_u_right_parenthesis = [v_H
        //           M⁻¹f_H]
        H_left_parenthesis_u_right_parenthesis = H_left_parenthesis_u_right_parenthesis_0;
        Eigen::MatrixXd J_G_0(2*n, 2*n);
        J_G_0 << Eigen::MatrixXd::Zero(n, n), U_s * U_s.transpose() * M,
        -U_s * U_s.transpose() * K * U_s * U_s.transpose() * M, Eigen::MatrixXd::Zero(n, n);
        // J_G = [0    U_sU_s^TM
        //       -U_sU_s^TKU_sU_s^TM 0 ]
        J_G = J_G_0;
        Eigen::MatrixXd J_H_0(2*n, 2*n);
        J_H_0 << Eigen::MatrixXd::Zero(n, n), Eigen::MatrixXd::Identity(n, n),
        -M.colPivHouseholderQr().solve(K), Eigen::MatrixXd::Zero(n, n);
        // J_H =  [0     I_n
        //             -M⁻¹K 0] - J_G 
        J_H = J_H_0 - J_G;
    }
};

struct second {
    Eigen::MatrixXd J_G;
    Eigen::MatrixXd Y_1;
    Eigen::MatrixXd Z_1;
    Eigen::MatrixXd Y_2;
    Eigen::MatrixXd Z_2;
    second(
        const Eigen::MatrixXd & U_s,
        const Eigen::MatrixXd & M,
        const Eigen::MatrixXd & K)
    {
        const long n = U_s.rows();
        const long s = U_s.cols();
        assert( U_s.rows() == n );
        assert( M.rows() == n );
        assert( M.cols() == n );
        assert( K.rows() == n );
        assert( K.cols() == n );
        Eigen::MatrixXd Y_1_0(n + 1, s);
        Y_1_0 << U_s,
        Eigen::MatrixXd::Zero(1, s);
        // Y_1 =  [U_s
        //              0]
        Y_1 = Y_1_0;
        Eigen::MatrixXd Z_1_0(n + 1, s);
        Z_1_0 << Eigen::MatrixXd::Zero(1, s),
        M * U_s;
        // Z_1 =  [ 0
        //              MU_s] 
        Z_1 = Z_1_0;
        Eigen::MatrixXd Y_2_0(n + 1, s);
        Y_2_0 << Eigen::MatrixXd::Zero(1, s),
        -U_s * U_s.transpose() * K * U_s;
        // Y_2 =  [0
        //             -U_sU_s^TKU_s]
        Y_2 = Y_2_0;
        Eigen::MatrixXd Z_2_0(n + 1, s);
        Z_2_0 << M * U_s,
        Eigen::MatrixXd::Zero(1, s);
        // Z_2 =  [ MU_s
        //              0] 
        Z_2 = Z_2_0;
        // J_G = Y_1Z_1^T + Y_2Z_2^T 
        J_G = Y_1 * Z_1.transpose() + Y_2 * Z_2.transpose();
    }
};

