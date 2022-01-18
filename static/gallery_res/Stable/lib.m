function output = Stable(F, mu, I_C, _lambda, f_0, f_1, f_2)
% output = Stable(F, μ, `$I_C$`, λ, `$f_0$`, `$f_1$`, `$f_2$`)
%
%    α = 1 + μ/λ - μ/(4λ)
%    J =`$f_0$`⋅(`$f_1$`×`$f_2$`)
%    `$P(F)$` = μ(1-1/(`$I_C$`+1))F + λ(J- α)`$\frac{\partial J}{\partial F}$`
%    
%    where 
%    F: ℝ^(3×3)
%    μ: ℝ 
%    `$I_C$`: ℝ 
%    λ: ℝ  
%    
%    `$\frac{\partial J}{\partial F}$` = [`$f_1$`×`$f_2$` `$f_2$`×`$f_0$` `$f_0$`×`$f_1$`]
%    
%    where 
%    `$f_0$`: ℝ^3 
%    `$f_1$`: ℝ^3 
%    `$f_2$`: ℝ^3 
%    
%    
    if nargin==0
        warning('generating random input data');
        [F, mu, I_C, _lambda, f_0, f_1, f_2] = generateRandomData();
    end
    function [F, mu, I_C, _lambda, f_0, f_1, f_2] = generateRandomData()
        mu = randn();
        I_C = randn();
        _lambda = randn();
        F = randn(3, 3);
        f_0 = randn(3,1);
        f_1 = randn(3,1);
        f_2 = randn(3,1);
    end

    f_0 = reshape(f_0,[],1);
    f_1 = reshape(f_1,[],1);
    f_2 = reshape(f_2,[],1);

    assert( isequal(size(F), [3, 3]) );
    assert(numel(mu) == 1);
    assert(numel(I_C) == 1);
    assert(numel(_lambda) == 1);
    assert( numel(f_0) == 3 );
    assert( numel(f_1) == 3 );
    assert( numel(f_2) == 3 );

    % α = 1 + μ/λ - μ/(4λ)
    alpha = 1 + mu / _lambda - mu / (4 * _lambda);
    % J =`$f_0$`⋅(`$f_1$`×`$f_2$`)
    J = dot(f_0,(cross(f_1, f_2)));
    % `$\frac{\partial J}{\partial F}$` = [`$f_1$`×`$f_2$` `$f_2$`×`$f_0$` `$f_0$`×`$f_1$`]
    fracpartialJpartialF_0 = [[reshape(cross(f_1, f_2), [3, 1]), reshape(cross(f_2, f_0), [3, 1]), reshape(cross(f_0, f_1), [3, 1])]];
    fracpartialJpartialF = fracpartialJpartialF_0;
    % `$P(F)$` = μ(1-1/(`$I_C$`+1))F + λ(J- α)`$\frac{\partial J}{\partial F}$`
    PF = mu * (1 - 1 / (I_C + 1)) * F + _lambda * (J - alpha) * fracpartialJpartialF;
    output.alpha = alpha;
    output.J = J;
    output.PF = PF;
    output.fracpartialJpartialF = fracpartialJpartialF;
end

