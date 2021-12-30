#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct first {

    Eigen::Matrix<double, 3, 3> I(
        const Eigen::Matrix<double, 2, 1> & x)
    {
    
    

    
        return (R(x)).cwiseProduct(S(x));    
    }
    first(
        const std::function<Eigen::Matrix<double, 3, 3>(Eigen::Matrix<double, 2, 1>)> & R,
        const std::function<Eigen::Matrix<double, 3, 3>(Eigen::Matrix<double, 2, 1>)> & S)
    {
    
    }
};

struct second {

    Eigen::Matrix<double, 3, 3> I(
        const Eigen::Matrix<double, 2, 1> & x)
    {
    
    

        Eigen::MatrixXd sum_0 = Eigen::MatrixXd::Zero(3, 3);
        for(int k=1; k<=b.size(); k++){
            sum_0 += b.at(k-1) * T.at(k-1)(x);
        }
        return (R(x)).cwiseProduct(sum_0);    
    }
    double 𝐸_left_curly_bracket_data_right_curly_bracket(
        const Eigen::Matrix<double, 2, 1> & X)
    {
    
    

        Eigen::MatrixXd sum_1 = Eigen::MatrixXd::Zero(3, 3);
        for(int k=1; k<=b.size(); k++){
            sum_1 += b.at(k-1) * T.at(k-1)(X);
        }
        return λ_left_curly_bracket_data_right_curly_bracket * pow((I(X) - (R(X)).cwiseProduct(sum_1)).norm(), 2);    
    }
    double 𝐸_left_curly_bracket_clustering_right_curly_bracket(
        const Eigen::Matrix<double, 2, 1> & X)
    {
    
    

    
        return λ_left_curly_bracket_clustering_right_curly_bracket * pow((r(X) - r_left_curly_bracket_cluster_right_curly_bracket(X)).norm(), 2);    
    }
    double 𝐸_left_curly_bracket_r_hyphen_minus_sparsity_right_curly_bracket(
        const Eigen::Matrix<double, 2, 1> & X)
    {
    
    

    
        return λ_left_curly_bracket_r_hyphen_minus_sparsity_right_curly_bracket * pow((r(X)).norm(), 2);    
    }
    double 𝐸_left_curly_bracket_reflectance_right_curly_bracket(
        const Eigen::Matrix<double, 2, 1> & X)
    {
    
    

    
        return 𝐸_left_curly_bracket_clustering_right_curly_bracket(X) + 𝐸_left_curly_bracket_r_hyphen_minus_sparsity_right_curly_bracket(X) + 𝐸_left_curly_bracket_r_hyphen_minus_consistency_right_curly_bracket(X);    
    }
    double 𝐸_left_curly_bracket_decomp_right_curly_bracket(
        const Eigen::Matrix<double, 2, 1> & X)
    {
    
    

    
        return 𝐸_left_curly_bracket_data_right_curly_bracket(X) + 𝐸_left_curly_bracket_reflectance_right_curly_bracket(X) + 𝐸_left_curly_bracket_illumination_right_curly_bracket(X);    
    }
    second(
        const std::function<Eigen::Matrix<double, 3, 3>(Eigen::Matrix<double, 2, 1>)> & R,
        const std::vector<std::function<Eigen::Matrix<double, 3, 3>(Eigen::Matrix<double, 2, 1>)>> & T,
        const std::vector<double> & b,
        const std::function<double(Eigen::Matrix<double, 2, 1>)> & 𝐸_left_curly_bracket_illumination_right_curly_bracket,
        const double & λ_left_curly_bracket_data_right_curly_bracket,
        const std::function<double(Eigen::Matrix<double, 2, 1>)> & 𝐸_left_curly_bracket_r_hyphen_minus_consistency_right_curly_bracket,
        const double & λ_left_curly_bracket_clustering_right_curly_bracket,
        const std::function<Eigen::Matrix<double, 3, 3>(Eigen::Matrix<double, 2, 1>)> & r_left_curly_bracket_cluster_right_curly_bracket,
        const std::function<Eigen::Matrix<double, 3, 3>(Eigen::Matrix<double, 2, 1>)> & r,
        const double & λ_left_curly_bracket_r_hyphen_minus_sparsity_right_curly_bracket)
    {
        const long dim_0 = b.size();
        assert( T.size() == dim_0 );
    
    }
};

