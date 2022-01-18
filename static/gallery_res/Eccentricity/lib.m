function output = eccentricity(m, italic_p, italic_f_italic_s_0, a)
% output = eccentricity(m, 𝑝, `$𝑓_{𝑠0}$`, a)
%
%     𝑙(𝐿) = π𝑑(𝐿)^2/4 ⋅ 𝐿 where 𝐿 : ℝ
%    `$𝑙_0$` = 1488
%    q = (5.71 ⋅ 10^(-6), -1.78 ⋅ 10^(-4), 0.204)
%    cos, sin from trigonometry
%    
%    𝑔(x,`$x_0$`, 𝜃,𝜎,`$𝑓_𝑠$`) = exp(-||x-`$x_0$`||^2/(2𝜎^2)) cos(2π`$𝑓_𝑠$`x ⋅[cos(𝜃) sin(𝜃)]) where x: ℝ^n,`$x_0$`: ℝ^n,`$𝑓_𝑠$`: ℝ, 𝜎 : ℝ, 𝜃 : ℝ
%    
%    
%    Ψ(𝑒, `$𝑓_𝑠$`)= m(0, 𝑝₀ 𝜏(`$𝑓_𝑠$`) +𝑝₁ 𝜏(`$𝑓_𝑠$`)+𝑝₂ + (𝑝₃ 𝜏(`$𝑓_𝑠$`)^2 + 𝑝₄ 𝜏(`$𝑓_𝑠$`) +𝑝₅)⋅ 𝜁(`$𝑓_𝑠$`)𝑒 + (𝑝₆ 𝜏(`$𝑓_𝑠$`)^2 +𝑝₇ 𝜏(`$𝑓_𝑠$`) + 𝑝₈)⋅𝜁(`$𝑓_𝑠$`)𝑒^2) where `$𝑓_𝑠$` : ℝ, 𝑒: ℝ
%    
%    𝜁(`$𝑓_𝑠$`) = exp(𝑝₉ 𝜏(`$𝑓_𝑠$`)) - 1 where `$𝑓_𝑠$` : ℝ
%    𝜏(`$𝑓_𝑠$`) = m(log(`$𝑓_𝑠$`)-log(`$𝑓_{𝑠0}$`), 0) where `$𝑓_𝑠$` : ℝ
%    where
%    m: ℝ, ℝ -> ℝ
%    𝑝: ℝ^10
%    `$𝑓_{𝑠0}$`: ℝ
%    
%    
%    𝐴(𝑒)= ln(64) 2.3/(0.106⋅(𝑒+2.3) )  where 𝑒 : ℝ
%    
%    
%    
%    𝑑(𝐿)= 7.75-5.75((𝐿a/846)^0.41/((𝐿a/846)^0.41 + 2))  where 𝐿 : ℝ
%    where 
%    a : ℝ
%    
%    
%     `$\hat{Ψ}$`(𝑒, `$𝑓_𝑠$`, 𝐿) = (𝑠(𝑒, `$𝑓_𝑠$`) ⋅ (log_10(𝑙(𝐿)/`$𝑙_0$`)) + 1) Ψ(𝑒, `$𝑓_𝑠$`) where `$𝑓_𝑠$` : ℝ, 𝑒: ℝ, 𝐿 : ℝ
%    
%    
%    𝑠(𝑒,`$𝑓_𝑠$`) = 𝜁(`$𝑓_𝑠$`)(q_0 𝑒^2 + q_1 𝑒) + q_2  where `$𝑓_𝑠$` : ℝ, 𝑒: ℝ
%    
    if nargin==0
        warning('generating random input data');
        [m, italic_p, italic_f_italic_s_0, a] = generateRandomData();
    end
    function [m, italic_p, italic_f_italic_s_0, a] = generateRandomData()
        italic_f_italic_s_0 = randn();
        a = randn();
        m = @mFunc;
        rseed = randi(2^32);
        function tmp =  mFunc(p0, p1)
            rng(rseed);
            tmp = randn();
        end

        italic_p = randn(10,1);
    end

    italic_p = reshape(italic_p,[],1);

    assert( numel(italic_p) == 10 );
    assert(numel(italic_f_italic_s_0) == 1);
    assert(numel(a) == 1);

    % `$𝑙_0$` = 1488
    italic_l_0 = 1488;
    % q = (5.71 ⋅ 10^(-6), -1.78 ⋅ 10^(-4), 0.204)
    q = [5.71 * 10.^(-6); -1.78 * 10.^(-4); 0.204];
    function ret = italic_g(x, `$x_0$`, 𝜃, 𝜎, `$𝑓_italic_s$`)
        x = reshape(x,[],1);
        `$x_0$` = reshape(`$x_0$`,[],1);
        n = size(x, 1);
        assert( numel(x) == n );
        assert( numel(`$x_0$`) == n );
        assert(numel(𝜃) == 1);
        assert(numel(𝜎) == 1);
        assert(numel(`$𝑓_italic_s$`) == 1);

        italic_g_0 = zeros(1, 2);
        italic_g_0(1,:) = [cos(𝜃), sin(𝜃)];
        ret = exp(-norm(x - `$x_0$`, 2).^2 / (2 * 𝜎.^2)) * cos(reshape(2 * pi * `$𝑓_italic_s$` * x, [n, 1]) * italic_g_0);
    end

    function ret = italic_tau(`$𝑓_italic_s$`)
        assert(numel(`$𝑓_italic_s$`) == 1);

        ret = m(log(`$𝑓_italic_s$`) - log(italic_f_italic_s_0), 0);
    end

    function ret = italic_zeta(`$𝑓_italic_s$`)
        assert(numel(`$𝑓_italic_s$`) == 1);

        ret = exp(italic_p(9) * italic_tau(`$𝑓_italic_s$`)) - 1;
    end

    function ret = Psi(𝑒, `$𝑓_italic_s$`)
        assert(numel(𝑒) == 1);
        assert(numel(`$𝑓_italic_s$`) == 1);

        ret = m(0, italic_p(0) * italic_tau(`$𝑓_italic_s$`) + italic_p(1) * italic_tau(`$𝑓_italic_s$`) + italic_p(2) + (italic_p(3) * italic_tau(`$𝑓_italic_s$`).^2 + italic_p(4) * italic_tau(`$𝑓_italic_s$`) + italic_p(5)) * italic_zeta(`$𝑓_italic_s$`) * 𝑒 + (italic_p(6) * italic_tau(`$𝑓_italic_s$`).^2 + italic_p(7) * italic_tau(`$𝑓_italic_s$`) + italic_p(8)) * italic_zeta(`$𝑓_italic_s$`) * 𝑒.^2);
    end

    function ret = italic_A(𝑒)
        assert(numel(𝑒) == 1);

        ret = log(64) * 2.3 / (0.106 * (𝑒 + 2.3));
    end

    function ret = italic_d(𝐿)
        assert(numel(𝐿) == 1);

        ret = 7.75 - 5.75 * ((𝐿 * a / 846).^0.41 / ((𝐿 * a / 846).^0.41 + 2));
    end

    function ret = italic_l(𝐿)
        assert(numel(𝐿) == 1);

        ret = pi * italic_d(𝐿).^2 / 4 * 𝐿;
    end

    function ret = italic_s(𝑒, `$𝑓_italic_s$`)
        assert(numel(𝑒) == 1);
        assert(numel(`$𝑓_italic_s$`) == 1);

        ret = italic_zeta(`$𝑓_italic_s$`) * (q(0) * 𝑒.^2 + q(1) * 𝑒) + q(2);
    end

    function ret = hatPsi(𝑒, `$𝑓_italic_s$`, 𝐿)
        assert(numel(𝑒) == 1);
        assert(numel(`$𝑓_italic_s$`) == 1);
        assert(numel(𝐿) == 1);

        ret = (italic_s(𝑒, `$𝑓_italic_s$`) * (log10(italic_l(𝐿) / italic_l_0)) + 1) * Psi(𝑒, `$𝑓_italic_s$`);
    end

    output.italic_l_0 = italic_l_0;
    output.q = q;
    output.italic_l = @italic_l;
    output.italic_g = @italic_g;
    output.Psi = @Psi;
    output.italic_zeta = @italic_zeta;
    output.italic_tau = @italic_tau;
    output.italic_A = @italic_A;
    output.italic_d = @italic_d;
    output.hatPsi = @hatPsi;
    output.italic_s = @italic_s;
output.italic_f_italic_s_0 = italic_f_italic_s_0;    
output.italic_p = italic_p;    
output.a = a;
end

