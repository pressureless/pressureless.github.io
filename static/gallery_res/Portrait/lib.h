#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct first {
    Eigen::MatrixXd I_left_curly_bracket_out_right_curly_bracket;
    Eigen::MatrixXd I;
    double M_c;
    Eigen::Matrix<double, 3, 1> l_left_curly_bracket_fill_right_curly_bracket;
    Eigen::MatrixXd D;
    Eigen::MatrixXd w;
    double L_left_curly_bracket_feat_right_curly_bracket(
        const double & θ)
    {
        double sum_2 = 0;
        for(int d=1; d<=λ.size(); d++){
            sum_2 += λ.at(d-1) * (Φ.at(d-1)(I_asterisk_operator) - Φ.at(d-1)(f(I_left_curly_bracket_in_right_curly_bracket, θ))).norm();
        }
        return sum_2;    
    }
    double L_left_curly_bracket_pix_right_curly_bracket(
        const double & θ)
    {
        return (I_asterisk_operator - f(I_left_curly_bracket_in_right_curly_bracket, θ)).norm();    
    }
    double L(
        const double & θ)
    {
        return 0.01 * L_left_curly_bracket_feat_right_curly_bracket(θ) + L_left_curly_bracket_pix_right_curly_bracket(θ);    
    }
    first(
        const Eigen::MatrixXd & I_left_curly_bracket_in_right_curly_bracket,
        const Eigen::MatrixXd & A,
        const Eigen::MatrixXd & B,
        const Eigen::MatrixXd & I_l,
        const Eigen::MatrixXd & I_s,
        const Eigen::MatrixXd & M,
        const double & M_left_curly_bracket_in_right_curly_bracket,
        const std::function<double(double)> & G,
        const std::vector<double> & σ_left_curly_bracket_c_comma_right_curly_bracket,
        const std::vector<double> & w_left_curly_bracket_c_comma_right_curly_bracket,
        const Eigen::Matrix<double, 3, 1> & l_left_curly_bracket_key_right_curly_bracket,
        const Eigen::Matrix<double, 3, 1> & n,
        const std::vector<double> & x,
        const std::vector<double> & y,
        const std::vector<double> & u,
        const std::vector<double> & v,
        const std::vector<double> & σ,
        const std::vector<double> & λ,
        const Eigen::MatrixXd & I_asterisk_operator,
        const std::vector<std::function<Eigen::MatrixXd(Eigen::MatrixXd)>> & Φ,
        const std::function<Eigen::MatrixXd(Eigen::MatrixXd, double)> & f)
    {
        const long dim_0 = σ_left_curly_bracket_c_comma_right_curly_bracket.size();
        const long dim_1 = x.size();
        const long dim_2 = u.size();
        const long dim_3 = λ.size();
        const long p = I_left_curly_bracket_in_right_curly_bracket.rows();
        const long q = I_left_curly_bracket_in_right_curly_bracket.cols();
        assert( I_left_curly_bracket_in_right_curly_bracket.rows() == p );
        assert( A.rows() == p );
        assert( A.cols() == q );
        assert( B.rows() == p );
        assert( B.cols() == q );
        assert( I_l.rows() == p );
        assert( I_l.cols() == q );
        assert( I_s.rows() == p );
        assert( I_s.cols() == q );
        assert( M.rows() == p );
        assert( M.cols() == q );
        assert( w_left_curly_bracket_c_comma_right_curly_bracket.size() == dim_0 );
        assert( y.size() == dim_1 );
        assert( v.size() == dim_2 );
        assert( σ.size() == dim_2 );
        assert( I_asterisk_operator.rows() == p );
        assert( I_asterisk_operator.cols() == q );
        assert( Φ.size() == dim_3 );
        // I_left_curly_bracket_out_right_curly_bracket =I_left_curly_bracket_in_right_curly_bracket∘ A+B
        I_left_curly_bracket_out_right_curly_bracket = (I_left_curly_bracket_in_right_curly_bracket).cwiseProduct(A) + B;
        // I =I_l∘ (1_p,q - M)+I_s∘M
        I = (I_l).cwiseProduct((Eigen::MatrixXd::Ones(p, q) - M)) + (I_s).cwiseProduct(M);
        double sum_0 = 0;
        for(int k=1; k<=σ_left_curly_bracket_c_comma_right_curly_bracket.size(); k++){
            sum_0 += M_left_curly_bracket_in_right_curly_bracket * G(σ_left_curly_bracket_c_comma_right_curly_bracket.at(k-1)) * w_left_curly_bracket_c_comma_right_curly_bracket.at(k-1);
        }
        // M_c =sum_k M_left_curly_bracket_in_right_curly_bracket G(σ_left_curly_bracket_c_comma_right_curly_bracket_k)w_left_curly_bracket_c_comma_right_curly_bracket_k
        M_c = sum_0;
        // l_left_curly_bracket_fill_right_curly_bracket = 2(l_left_curly_bracket_key_right_curly_bracket⋅n)n - l_left_curly_bracket_key_right_curly_bracket
        l_left_curly_bracket_fill_right_curly_bracket = 2 * ((l_left_curly_bracket_key_right_curly_bracket).dot(n)) * n - l_left_curly_bracket_key_right_curly_bracket;
        // D_i,j = (x_i - u_j)^2 + (y_i - v_j)^2
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(dim_1, dim_2);
        for( int i=1; i<=dim_1; i++){
            for( int j=1; j<=dim_2; j++){
                D(i-1, j-1) = pow((x.at(i-1) - u.at(j-1)), 2) + pow((y.at(i-1) - v.at(j-1)), 2);
            }
        }
        // w_i,j = exp(-D_i,j/σ_j)/( sum__$j^{\prime}$_ exp(-D_i,_$j^{\prime}$_/σ__$j^{\prime}$_))
        Eigen::MatrixXd w = Eigen::MatrixXd::Zero(dim_1, dim_2);
        for( int i=1; i<=dim_1; i++){
            for( int j=1; j<=dim_2; j++){
                double sum_1 = 0;
                for(int _$j^{\prime}$_=1; _$j^{\prime}$_<=D.cols(); _$j^{\prime}$_++){
                    sum_1 += exp(-D(i-1, _$j^{\prime}$_-1) / double(σ.at(_$j^{\prime}$_-1)));
                }
                w(i-1, j-1) = exp(-D(i-1, j-1) / double(σ.at(j-1))) / double((sum_1));
            }
        }
    }
};

