#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct perceptual {
    double increment_Q_P;
    double increment_Q_O;
    double m(
        const double & ω,
        const double & σ)
    {
    
    

    
        return exp(-2 * pow(M_PI, 2) * pow(ω, 2) * pow(σ, 2));    
    }
    double m̃(
        const double & ω,
        const double & σ)
    {
    
    

    
        return CSF(ω) * m(ω, σ);    
    }
    double E_b(
        const double & σ)
    {
    
    

        double sum_0 = 0;
        for(int i=1; i<=ω.size(); i++){
            sum_0 += pow(((m̃(ω.at(i-1), σ)) / double(m_combining_tilde_left_curly_bracket_t_comma_b_right_curly_bracket)), β_b);
        }
        return sum_0;    
    }
    perceptual(
        const std::function<double(double)> & CSF,
        const std::vector<double> & ω,
        const double & m_combining_tilde_left_curly_bracket_t_comma_b_right_curly_bracket,
        const double & β_b,
        const double & σ_P_circumflex_accent_A,
        const double & σ_P_circumflex_accent_B,
        const double & σ_O_circumflex_accent_A,
        const double & σ_O_circumflex_accent_B)
    {
        const long dim_0 = ω.size();
        // increment_Q_P = E_b(σ_P_circumflex_accent_A) - E_b(σ_P_circumflex_accent_B)
        increment_Q_P = E_b(σ_P_circumflex_accent_A) - E_b(σ_P_circumflex_accent_B);
        // increment_Q_O = E_b(σ_O_circumflex_accent_A) - E_b(σ_O_circumflex_accent_B)
        increment_Q_O = E_b(σ_O_circumflex_accent_A) - E_b(σ_O_circumflex_accent_B);
    }
};

