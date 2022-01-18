#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct portrait {
    Eigen::MatrixXd I_out;
    Eigen::MatrixXd I;
    double M_c;
    Eigen::Matrix<double, 3, 1> vecl_fill;
    Eigen::MatrixXd D;
    Eigen::MatrixXd w;
    Eigen::VectorXd σ;
    std::vector<double> λ;
    Eigen::MatrixXd I_circumflex_accent_asterisk;
    Eigen::MatrixXd I_in;
    double L_feat(
        const double & θ)
    {
        double sum_2 = 0;
        for(int d=1; d<=λ.size(); d++){
            sum_2 += λ.at(d-1) * (Φ.at(d-1)(I_circumflex_accent_asterisk) - Φ.at(d-1)(f(I_in, θ))).norm();
        }
        return sum_2;    
    }
    double L_pix(
        const double & θ)
    {
        return (I_circumflex_accent_asterisk - f(I_in, θ)).norm();    
    }
    double L(
        const double & θ)
    {
        return 0.01 * L_feat(θ) + L_pix(θ);    
    }
    portrait(
        const Eigen::MatrixXd & I_in,
        const Eigen::MatrixXd & A,
        const Eigen::MatrixXd & B,
        const Eigen::MatrixXd & I_l,
        const Eigen::MatrixXd & I_s,
        const Eigen::MatrixXd & M,
        const double & M_in,
        const std::function<double(double)> & G,
        const std::vector<double> & σ_c_comma,
        const std::vector<double> & w_c_comma,
        const Eigen::Matrix<double, 3, 1> & vecl_key,
        const Eigen::Matrix<double, 3, 1> & vecn,
        const std::vector<double> & x,
        const std::vector<double> & y,
        const std::vector<double> & u,
        const std::vector<double> & v,
        const std::function<double(double, double)> & select,
        const double & K_σ,
        const double & j_circumflex_accent_prime,
        const std::vector<double> & λ,
        const Eigen::MatrixXd & I_circumflex_accent_asterisk,
        const std::vector<std::function<Eigen::MatrixXd(Eigen::MatrixXd)>> & Φ,
        const std::function<Eigen::MatrixXd(Eigen::MatrixXd, double)> & f)
    {
        const long dim_0 = σ_c_comma.size();
        const long dim_1 = x.size();
        const long dim_2 = u.size();
        const long dim_3 = λ.size();
        const long p = I_in.rows();
        const long q = I_in.cols();
        assert( I_in.rows() == p );
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
        assert( w_c_comma.size() == dim_0 );
        assert( y.size() == dim_1 );
        assert( v.size() == dim_2 );
        assert( I_circumflex_accent_asterisk.rows() == p );
        assert( I_circumflex_accent_asterisk.cols() == q );
        assert( Φ.size() == dim_3 );
        this->λ = λ;
        this->I_circumflex_accent_asterisk = I_circumflex_accent_asterisk;
        this->I_in = I_in;
        // `$I_{out}$` =`$I_{in}$`∘ A+B
        I_out = (I_in).cwiseProduct(A) + B;
        // I =`$I_l$`∘ (1_p,q - M)+`$I_s$`∘M
        I = (I_l).cwiseProduct((Eigen::MatrixXd::Ones(p, q) - M)) + (I_s).cwiseProduct(M);
        double sum_0 = 0;
        for(int k=1; k<=σ_c_comma.size(); k++){
            sum_0 += M_in * G(σ_c_comma.at(k-1)) * w_c_comma.at(k-1);
        }
        // `$M_c$` =sum_k `$M_{in}$` G(`$σ_{c,}$`_k)`$w_{c,}$`_k
        M_c = sum_0;
        // `$\vec{l}_{fill}$` = 2(`$\vec{l}_{key}$`⋅`$\vec{n}$`)`$\vec{n}$` - `$\vec{l}_{key}$`
        vecl_fill = 2 * ((vecl_key).dot(vecn)) * vecn - vecl_key;
        // D_i,j = (x_i - u_j)^2 + (y_i - v_j)^2
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(dim_1, dim_2);
        for( int i=1; i<=dim_1; i++){
            for( int j=1; j<=dim_2; j++){
                D(i-1, j-1) = pow((x.at(i-1) - u.at(j-1)), 2) + pow((y.at(i-1) - v.at(j-1)), 2);
            }
        }
        // σ_j = select((u_j - u_`$j^{′}$`)^2+(v_j - v_`$j^{′}$`)^2, `$K_σ$`)
        σ.resize(dim_2);
        for( int j=1; j<=dim_2; j++){
            σ[j-1] = select(pow((u.at(j-1) - u.at(j_circumflex_accent_prime-1)), 2) + pow((v.at(j-1) - v.at(j_circumflex_accent_prime-1)), 2), K_σ);
        }
        // w_i,j = exp(-D_i,j/σ_j)/( sum_`$j^{\prime}$` exp(-D_i,`$j^{\prime}$`/σ_`$j^{\prime}$`))
        Eigen::MatrixXd w = Eigen::MatrixXd::Zero(dim_1, dim_2);
        for( int i=1; i<=dim_1; i++){
            for( int j=1; j<=dim_2; j++){
                double sum_1 = 0;
                for(int _$j^{\prime}$_=1; _$j^{\prime}$_<=D.cols(); _$j^{\prime}$_++){
                    sum_1 += exp(-D(i-1, _$j^{\prime}$_-1) / double(σ[_$j^{\prime}$_-1]));
                }
                w(i-1, j-1) = exp(-D(i-1, j-1) / double(σ[j-1])) / double((sum_1));
            }
        }
    }
};

