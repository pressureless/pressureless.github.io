#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct elastic {
    double E;
    double E_q;
    double E_a;
    double E_n;
    double E_p;
    elastic(
        const double & E_r,
        const Eigen::VectorXd & q_g,
        const Eigen::VectorXd & q_h,
        const Eigen::VectorXd & m_g,
        const Eigen::VectorXd & m_h,
        const std::function<Eigen::VectorXd(Eigen::VectorXd, Eigen::VectorXd)> & reverse_solidus_angle,
        const double & λ_left_curly_bracket_q_comma_1_right_curly_bracket,
        const double & λ_left_curly_bracket_q_comma_2_right_curly_bracket,
        const double & t,
        const double & λ_left_curly_bracket_a_comma_1_right_curly_bracket,
        const double & λ_left_curly_bracket_a_comma_2_right_curly_bracket,
        const Eigen::VectorXd & q,
        const Eigen::VectorXd & q_a,
        const Eigen::VectorXd & m,
        const Eigen::VectorXd & m_a,
        const double & δ_circumflex_accent_left_curly_bracket_left_parenthesis_minus_sign_right_parenthesis_right_curly_bracket,
        const double & δ_circumflex_accent_left_curly_bracket_left_parenthesis_plus_sign_right_parenthesis_right_curly_bracket,
        const double & β_q,
        const double & β_circumflex_accent_left_curly_bracket_left_parenthesis_minus_sign_right_parenthesis_right_curly_bracket,
        const double & β_circumflex_accent_left_curly_bracket_left_parenthesis_plus_sign_right_parenthesis_right_curly_bracket,
        const double & μ,
        const double & ε)
    {
        const long n = q_g.size();
        assert( q_h.size() == n );
        assert( m_g.size() == n );
        assert( m_h.size() == n );
        assert( q.size() == n );
        assert( q_a.size() == n );
        assert( m.size() == n );
        assert( m_a.size() == n );
        // E_q = λ_left_curly_bracket_q_comma_1_right_curly_bracket||q_g-q_h+tm_g||^2 + λ_left_curly_bracket_q_comma_1_right_curly_bracket||q_h-q_g+tm_h||^2 + λ_left_curly_bracket_q_comma_2_right_curly_bracket||reverse_solidus_angle(m_g,m_h)||^2 
        E_q = λ_left_curly_bracket_q_comma_1_right_curly_bracket * pow((q_g - q_h + t * m_g).lpNorm<2>(), 2) + λ_left_curly_bracket_q_comma_1_right_curly_bracket * pow((q_h - q_g + t * m_h).lpNorm<2>(), 2) + λ_left_curly_bracket_q_comma_2_right_curly_bracket * pow((reverse_solidus_angle(m_g, m_h)).lpNorm<2>(), 2);
        // E_a = λ_left_curly_bracket_a_comma_1_right_curly_bracket||q-q_a||^2 + λ_left_curly_bracket_a_comma_2_right_curly_bracket||reverse_solidus_angle(m,m_a)||^2 
        E_a = λ_left_curly_bracket_a_comma_1_right_curly_bracket * pow((q - q_a).lpNorm<2>(), 2) + λ_left_curly_bracket_a_comma_2_right_curly_bracket * pow((reverse_solidus_angle(m, m_a)).lpNorm<2>(), 2);
        // E_n = δ_circumflex_accent_left_curly_bracket_left_parenthesis_minus_sign_right_parenthesis_right_curly_bracket(1/10 log((β_q-β_circumflex_accent_left_curly_bracket_left_parenthesis_minus_sign_right_parenthesis_right_curly_bracket)))^2 + δ_circumflex_accent_left_curly_bracket_left_parenthesis_plus_sign_right_parenthesis_right_curly_bracket(1/10 log((β_circumflex_accent_left_curly_bracket_left_parenthesis_plus_sign_right_parenthesis_right_curly_bracket-β_q)))^2
        E_n = δ_circumflex_accent_left_curly_bracket_left_parenthesis_minus_sign_right_parenthesis_right_curly_bracket * pow((1 / double(10) * log((β_q - β_circumflex_accent_left_curly_bracket_left_parenthesis_minus_sign_right_parenthesis_right_curly_bracket))), 2) + δ_circumflex_accent_left_curly_bracket_left_parenthesis_plus_sign_right_parenthesis_right_curly_bracket * pow((1 / double(10) * log((β_circumflex_accent_left_curly_bracket_left_parenthesis_plus_sign_right_parenthesis_right_curly_bracket - β_q))), 2);
        // E_p = (μ log((ε + β_q)))^2 + (μ log((ε + 1 - β_q)))^2
        E_p = pow((μ * log((ε + β_q))), 2) + pow((μ * log((ε + 1 - β_q))), 2);
        // E = E_r + E_q + E_a + E_n + E_p
        E = E_r + E_q + E_a + E_n + E_p;
    }
};

