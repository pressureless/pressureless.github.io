#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct first {

    double E_left_curly_bracket_pose_right_curly_bracket(
        const Eigen::Matrix<double, 3, 3> & S)
    {
    
    

    
        return E_left_curly_bracket_2D_right_curly_bracket(S) + E_left_curly_bracket_3D_right_curly_bracket(S) + E_left_curly_bracket_silhouette_right_curly_bracket(S) + E_left_curly_bracket_temporal_right_curly_bracket(S) + E_left_curly_bracket_anatomic_right_curly_bracket(S);    
    }
    double Ψ(
        const double & x)
    {
    
    

        double Ψ_ret;
        if(x > θ_left_curly_bracket_max_comma_i_right_curly_bracket){
            Ψ_ret = pow((x - θ_left_curly_bracket_max_comma_i_right_curly_bracket), 2);
        }
        else if(x < θ_left_curly_bracket_min_comma_i_right_curly_bracket){
            Ψ_ret = pow((θ_left_curly_bracket_min_comma_i_right_curly_bracket - x), 2);
        }
        else{
            Ψ = 0;
        }
        return Ψ_ret;    
    }
    first(
        const std::function<double(Eigen::Matrix<double, 3, 3>)> & E_left_curly_bracket_2D_right_curly_bracket,
        const std::function<double(Eigen::Matrix<double, 3, 3>)> & E_left_curly_bracket_3D_right_curly_bracket,
        const std::function<double(Eigen::Matrix<double, 3, 3>)> & E_left_curly_bracket_silhouette_right_curly_bracket,
        const std::function<double(Eigen::Matrix<double, 3, 3>)> & E_left_curly_bracket_temporal_right_curly_bracket,
        const std::function<double(Eigen::Matrix<double, 3, 3>)> & E_left_curly_bracket_anatomic_right_curly_bracket,
        const double & θ_left_curly_bracket_max_comma_i_right_curly_bracket,
        const double & θ_left_curly_bracket_min_comma_i_right_curly_bracket)
    {
    
    }
};

