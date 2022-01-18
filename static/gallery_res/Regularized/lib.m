function output = Regularized(f, s, F, r_epsilon, a, b, epsilon)
% output = Regularized(f, s, F, `$r_ε$`, a, b, ε)
%
%    placeholder = 1
%    
%    `$row_ε$`(r) = (15`$r_ε$`/8 + 1/`$r_ε$`³ )  where r ∈ ℝ^3
%    
%    
%    `$u_ε$`(r) = ((a-b)/`$r_ε$`I_3 + b/`$r_ε$`³ r r^T + a/2 ε²/`$r_ε$`³ I_3 ) f where r ∈ ℝ^3
%    
%    where
%    
%    f ∈ ℝ^3
%    
%    `$t_ε$`(r) = -a(1/`$r_ε$`³ + 3ε²/(2`$r_ε$`⁵) ) Fr where r ∈ ℝ^3
%    
%    
%    `$s_ε$`(r) = (2b-a)(1/`$r_ε$`³ + 3ε²/(2`$r_ε$`⁵))(sr) where r ∈ ℝ^3
%    
%    where
%    
%    s ∈ ℝ
%    
%    `$p_ε$`(r) = (2b-a)/`$r_ε$`³ Fr - 3/(2`$r_ε$`⁵)(2b(rᵀFr)I_3+aε²F)r where r ∈ ℝ^3
%    
%    where
%    
%    F ∈ ℝ^(3×3)
%    `$r_ε$` ∈ ℝ
%    a ∈ ℝ 
%    b ∈ ℝ
%    ε ∈ ℝ
%    
    if nargin==0
        warning('generating random input data');
        [f, s, F, r_epsilon, a, b, epsilon] = generateRandomData();
    end
    function [f, s, F, r_epsilon, a, b, epsilon] = generateRandomData()
        s = randn();
        r_epsilon = randn();
        a = randn();
        b = randn();
        epsilon = randn();
        f = randn(3,1);
        F = randn(3, 3);
    end

    f = reshape(f,[],1);

    assert( numel(f) == 3 );
    assert(numel(s) == 1);
    assert( isequal(size(F), [3, 3]) );
    assert(numel(r_epsilon) == 1);
    assert(numel(a) == 1);
    assert(numel(b) == 1);
    assert(numel(epsilon) == 1);

    % placeholder = 1
    placeholder = 1;
    function ret = row_epsilon(r)
        r = reshape(r,[],1);
        assert( numel(r) == 3 );

        ret = (15 * r_epsilon / 8 + 1 / r_epsilon.^3);
    end

    function ret = u_epsilon(r)
        r = reshape(r,[],1);
        assert( numel(r) == 3 );

        ret = ((a - b) / r_epsilon * speye(3) + reshape(b / r_epsilon.^3 * r, [3, 1]) * r' + a / 2 * epsilon.^2 / r_epsilon.^3 * speye(3)) * f;
    end

    function ret = t_epsilon(r)
        r = reshape(r,[],1);
        assert( numel(r) == 3 );

        ret = -a * (1 / r_epsilon.^3 + 3 * epsilon.^2 / (2 * r_epsilon.^5)) * F * r;
    end

    function ret = s_epsilon(r)
        r = reshape(r,[],1);
        assert( numel(r) == 3 );

        ret = (2 * b - a) * (1 / r_epsilon.^3 + 3 * epsilon.^2 / (2 * r_epsilon.^5)) * (s * r);
    end

    function ret = p_epsilon(r)
        r = reshape(r,[],1);
        assert( numel(r) == 3 );

        ret = (2 * b - a) / r_epsilon.^3 * F * r - 3 / (2 * r_epsilon.^5) * (2 * b * (r' * F * r) * speye(3) + a * epsilon.^2 * F) * r;
    end

    output.placeholder = placeholder;
    output.row_epsilon = @row_epsilon;
    output.u_epsilon = @u_epsilon;
    output.t_epsilon = @t_epsilon;
    output.s_epsilon = @s_epsilon;
    output.p_epsilon = @p_epsilon;
output.r_epsilon = r_epsilon;    
output.a = a;    
output.b = b;    
output.epsilon = epsilon;    
output.f = f;    
output.F = F;    
output.s = s;
end

