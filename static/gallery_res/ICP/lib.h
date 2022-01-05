#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct icp {
    Eigen::Matrix<double, 3, 1> ã;
    std::vector<Eigen::Matrix<double, 3, 1>> n;
    Eigen::Matrix<double, 3, 1> t̃;
    double reverse_solidus_varepsilon_left_curly_bracket_point_right_curly_bracket;
    double reverse_solidus_varepsilon_left_curly_bracket_plane_right_curly_bracket;
    double reverse_solidus_varepsilon_left_curly_bracket_symm_hyphen_minus_RN_right_curly_bracket;
    double reverse_solidus_varepsilon_left_curly_bracket_symm_right_curly_bracket;
    double reverse_solidus_varepsilon_left_curly_bracket_two_hyphen_minus_plane_right_curly_bracket;
    icp(
        const Eigen::Matrix<double, 3, 3> & R,
        const Eigen::Matrix<double, 3, 1> & a,
        const double & θ,
        const std::vector<Eigen::Matrix<double, 3, 1>> & p,
        const std::vector<Eigen::Matrix<double, 3, 1>> & q,
        const std::vector<Eigen::Matrix<double, 3, 1>> & n_q,
        const std::vector<Eigen::Matrix<double, 3, 1>> & n_p,
        const Eigen::Matrix<double, 3, 1> & t)
    {
        const long dim_0 = p.size();
        assert( q.size() == dim_0 );
        assert( n_q.size() == dim_0 );
        assert( n_p.size() == dim_0 );
        // ã = a tan(θ)
        ã = a * tan(θ);
        // n_i = n_q_i + n_p_i
        std::vector<Eigen::Matrix<double, 3, 1>> n(dim_0);
        for( int i=1; i<=dim_0; i++){
            n.at(i-1) = n_q.at(i-1) + n_p.at(i-1);
        }
        // t̃ = t/cos(θ)
        t̃ = t / double(cos(θ));
        double sum_0 = 0;
        for(int i=1; i<=p.size(); i++){
            sum_0 += (R * p.at(i-1) + t - q.at(i-1)).lpNorm<2>();
        }
        // reverse_solidus_varepsilon_left_curly_bracket_point_right_curly_bracket = ∑_i ||R p_i + t - q_i||
        reverse_solidus_varepsilon_left_curly_bracket_point_right_curly_bracket = sum_0;
        double sum_1 = 0;
        for(int i=1; i<=q.size(); i++){
            sum_1 += pow((((R * p.at(i-1) + t - q.at(i-1))).dot(n_q.at(i-1))), 2);
        }
        // reverse_solidus_varepsilon_left_curly_bracket_plane_right_curly_bracket = ∑_i ((R p_i + t - q_i) ⋅ n_q_i)^2
        reverse_solidus_varepsilon_left_curly_bracket_plane_right_curly_bracket = sum_1;
        double sum_2 = 0;
        for(int i=1; i<=q.size(); i++){
            sum_2 += pow((((R * p.at(i-1) + R.colPivHouseholderQr().solve(q.at(i-1)) + t)).dot((R * n_p.at(i-1) + R.colPivHouseholderQr().solve(n_q.at(i-1))))), 2);
        }
        // reverse_solidus_varepsilon_left_curly_bracket_symm_hyphen_minus_RN_right_curly_bracket = ∑_i ((R p_i + R⁻¹ q_i + t) ⋅ (Rn_p_i + R⁻¹n_q_i))^2
        reverse_solidus_varepsilon_left_curly_bracket_symm_hyphen_minus_RN_right_curly_bracket = sum_2;
        double sum_3 = 0;
        for(int i=1; i<=n.size(); i++){
            sum_3 += pow(cos(θ), 2) * pow((((p.at(i-1) - q.at(i-1))).dot(n.at(i-1)) + ((((p.at(i-1) + q.at(i-1))).cross(n.at(i-1)))).dot(ã) + (n.at(i-1)).dot(t̃)), 2);
        }
        // reverse_solidus_varepsilon_left_curly_bracket_symm_right_curly_bracket = ∑_i cos²(θ)((p_i - q_i)⋅n_i +((p_i+q_i)×n_i)⋅ã+n_i⋅t̃)² 
        reverse_solidus_varepsilon_left_curly_bracket_symm_right_curly_bracket = sum_3;
        double sum_4 = 0;
        for(int i=1; i<=q.size(); i++){
            sum_4 += (pow((((R * p.at(i-1) + R.colPivHouseholderQr().solve(q.at(i-1)) + t)).dot((R * n_p.at(i-1)))), 2) + pow((((R * p.at(i-1) + R.colPivHouseholderQr().solve(q.at(i-1)) + t)).dot((R.colPivHouseholderQr().solve(n_q.at(i-1))))), 2));
        }
        // reverse_solidus_varepsilon_left_curly_bracket_two_hyphen_minus_plane_right_curly_bracket = ∑_i(((R p_i + R⁻¹ q_i + t) ⋅ (R n_p_i))^2 + ((R p_i + R⁻¹ q_i + t) ⋅ (R⁻¹n_q_i))^2)
        reverse_solidus_varepsilon_left_curly_bracket_two_hyphen_minus_plane_right_curly_bracket = sum_4;
    }
};

