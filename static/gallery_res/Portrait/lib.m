function output = portrait(I_in, A, B, I_l, I_s, M, M_in, G, sigma_c, w_c, vecl_key, vecn, x, y, u, v, select, K_sigma, j_prime, _lambda, I_asterisk, Phi, f)
% output = portrait(`$I_{in}$`, A, B, `$I_l$`, `$I_s$`, M, `$M_{in}$`, G, `$σ_{c,}$`, `$w_{c,}$`, `$\vec{l}_{key}$`, `$\vec{n}$`, x, y, u, v, select, `$K_σ$`, `$j^{′}$`, λ, `$I^*$`, Φ, f)
%
%    `$I_{out}$` =`$I_{in}$`∘ A+B
%    
%    where
%    
%    `$I_{in}$`: ℝ^(p × q)
%    A: ℝ^(p × q)
%    B: ℝ^(p × q)
%    
%    I =`$I_l$`∘ (1_p,q - M)+`$I_s$`∘M
%    
%    where
%    
%    `$I_l$`: ℝ^(p × q)
%    `$I_s$`: ℝ^(p × q)
%    M: ℝ^(p × q)
%    
%    `$M_c$` =sum_k `$M_{in}$` G(`$σ_{c,}$`_k)`$w_{c,}$`_k
%    
%    where
%    
%    `$M_{in}$`: ℝ
%    G: ℝ -> ℝ
%    `$σ_{c,}$`_k: ℝ
%    `$w_{c,}$`_k: ℝ
%    
%    `$\vec{l}_{fill}$` = 2(`$\vec{l}_{key}$`⋅`$\vec{n}$`)`$\vec{n}$` - `$\vec{l}_{key}$`
%    
%    where
%    
%    `$\vec{l}_{key}$`: ℝ^3
%    `$\vec{n}$`: ℝ^3
%    
%    D_i,j = (x_i - u_j)^2 + (y_i - v_j)^2
%    
%    where
%    
%    x_i: ℝ
%    y_i: ℝ
%    u_j: ℝ
%    v_j: ℝ
%    
%    w_i,j = exp(-D_i,j/σ_j)/( sum_`$j^{\prime}$` exp(-D_i,`$j^{\prime}$`/σ_`$j^{\prime}$`))
%    
%    
%    σ_j = select((u_j - u_`$j^{′}$`)^2+(v_j - v_`$j^{′}$`)^2, `$K_σ$`)
%    
%    where
%    
%    select: ℝ, ℝ -> ℝ
%    `$K_σ$`: ℝ
%    `$j^{′}$`: ℝ
%    
%    `$L_{feat}$`(θ) = sum_d λ_d||Φ_d(`$I^*$`) - Φ_d(f(`$I_{in}$`;θ))|| where θ: ℝ
%    `$L_{pix}$`(θ) = ||`$I^*$` - f(`$I_{in}$`;θ)|| where θ: ℝ
%    L(θ) = 0.01 `$L_{feat}$`(θ) + `$L_{pix}$`(θ) where θ: ℝ
%    where
%    
%    λ_d: ℝ
%    `$I^*$`: ℝ^(p × q)
%    Φ_d: ℝ^(p × q) -> ℝ^(p × q)
%    f: ℝ^(p × q), ℝ -> ℝ^(p × q)
%    
    if nargin==0
        warning('generating random input data');
        [I_in, A, B, I_l, I_s, M, M_in, G, sigma_c, w_c, vecl_key, vecn, x, y, u, v, select, K_sigma, j_prime, _lambda, I_asterisk, Phi, f] = generateRandomData();
    end
    function [I_in, A, B, I_l, I_s, M, M_in, G, sigma_c, w_c, vecl_key, vecn, x, y, u, v, select, K_sigma, j_prime, _lambda, I_asterisk, Phi, f] = generateRandomData()
        M_in = randn();
        K_sigma = randn();
        j_prime = randn();
        dim_0 = randi(10);
        dim_1 = randi(10);
        dim_2 = randi(10);
        dim_3 = randi(10);
        p = randi(10);
        q = randi(10);
        I_in = randn(p, q);
        A = randn(p, q);
        B = randn(p, q);
        I_l = randn(p, q);
        I_s = randn(p, q);
        M = randn(p, q);
        G = @GFunc;
        rseed = randi(2___32);
        function tmp =  GFunc(p0)
            rng(rseed);
            tmp = randn();
        end

        sigma_c = randn(dim_0,1);
        w_c = randn(dim_0,1);
        vecl_key = randn(3,1);
        vecn = randn(3,1);
        x = randn(dim_1,1);
        y = randn(dim_1,1);
        u = randn(dim_2,1);
        v = randn(dim_2,1);
        select = @selectFunc;
        rseed = randi(2___32);
        function tmp =  selectFunc(p0, p1)
            rng(rseed);
            tmp = randn();
        end

        _lambda = randn(dim_3,1);
        I_asterisk = randn(p, q);
        Phi = _______;
        for i = 1:dim_3
            Phi_f = @(p0) randn(p,q);
            Phi_____end+1,1__ = Phi_f;
        end
        f = @fFunc;
        rseed = randi(2___32);
        function tmp =  fFunc(p0, p1)
            rng(rseed);
            tmp = randn(p,q);
        end

    end

    sigma_c = reshape(sigma_c,[],1);
    w_c = reshape(w_c,[],1);
    vecl_key = reshape(vecl_key,[],1);
    vecn = reshape(vecn,[],1);
    x = reshape(x,[],1);
    y = reshape(y,[],1);
    u = reshape(u,[],1);
    v = reshape(v,[],1);
    _lambda = reshape(_lambda,[],1);

    dim_0 = size(sigma_c, 1);
    dim_1 = size(x, 1);
    dim_2 = size(u, 1);
    dim_3 = size(_lambda, 1);
    p = size(I_in, 1);
    q = size(I_in, 2);
    assert( isequal(size(I_in), [p, q]) );
    assert( isequal(size(A), [p, q]) );
    assert( isequal(size(B), [p, q]) );
    assert( isequal(size(I_l), [p, q]) );
    assert( isequal(size(I_s), [p, q]) );
    assert( isequal(size(M), [p, q]) );
    assert(numel(M_in) == 1);
    assert( size(sigma_c,1) == dim_0 );
    assert( size(w_c,1) == dim_0 );
    assert( numel(vecl_key) == 3 );
    assert( numel(vecn) == 3 );
    assert( size(x,1) == dim_1 );
    assert( size(y,1) == dim_1 );
    assert( size(u,1) == dim_2 );
    assert( size(v,1) == dim_2 );
    assert(numel(K_sigma) == 1);
    assert(numel(j_prime) == 1);
    assert( size(_lambda,1) == dim_3 );
    assert( isequal(size(I_asterisk), [p, q]) );

    % `$I_{out}$` =`$I_{in}$`∘ A+B
    I_out = I_in.*A + B;
    % I =`$I_l$`∘ (1_p,q - M)+`$I_s$`∘M
    I = I_l.*(ones(p, q) - M) + I_s.*M;
    % `$M_c$` =sum_k `$M_{in}$` G(`$σ_{c,}$`_k)`$w_{c,}$`_k
    sum_0 = 0;
    for k = 1:size(sigma_c, 1)
        sum_0 = sum_0 + M_in * G(sigma_c(k)) * w_c(k);
    end
    M_c = sum_0;
    % `$\vec{l}_{fill}$` = 2(`$\vec{l}_{key}$`⋅`$\vec{n}$`)`$\vec{n}$` - `$\vec{l}_{key}$`
    vecl_fill = 2 * (dot(vecl_key,vecn)) * vecn - vecl_key;
    % D_i,j = (x_i - u_j)^2 + (y_i - v_j)^2
    D = zeros(dim_1, dim_2);
    for i = 1:dim_1
        for j = 1:dim_2
            D(i, j) = (x(i) - u(j)).___2 + (y(i) - v(j)).___2;
        end
    end
    % σ_j = select((u_j - u_`$j^{′}$`)^2+(v_j - v_`$j^{′}$`)^2, `$K_σ$`)
    sigma = zeros(dim_2,1);
    for j = 1:dim_2
        sigma(j) = select((u(j) - u(j_prime)).___2 + (v(j) - v(j_prime)).___2, K_sigma);
    end
    % w_i,j = exp(-D_i,j/σ_j)/( sum_`$j^{\prime}$` exp(-D_i,`$j^{\prime}$`/σ_`$j^{\prime}$`))
    w = zeros(dim_1, dim_2);
    for i = 1:dim_1
        for j = 1:dim_2
            sum_1 = 0;
            for _dollar_signj____________prime__dollar_sign_ = 1:size(D,2)
                sum_1 = sum_1 + exp(-D(i, _dollar_signj____________prime__dollar_sign_) / sigma(_dollar_signj____________prime__dollar_sign_));
            end
            w(i, j) = exp(-D(i, j) / sigma(j)) / (sum_1);
        end
    end
    function ret = L_feat(θ)
        assert(numel(θ) == 1);

        sum_2 = 0;
        for d = 1:size(_lambda, 1)
            sum_2 = sum_2 + _lambda(d) * norm(Phi_____d__(I_asterisk) - Phi_____d__(f(I_in, θ)), 'fro');
        end
        ret = sum_2;
    end

    function ret = L_pix(θ)
        assert(numel(θ) == 1);

        ret = norm(I_asterisk - f(I_in, θ), 'fro');
    end

    function ret = L(θ)
        assert(numel(θ) == 1);

        ret = 0.01 * L_feat(θ) + L_pix(θ);
    end

    output.I_out = I_out;
    output.I = I;
    output.M_c = M_c;
    output.vecl_fill = vecl_fill;
    output.D = D;
    output.w = w;
    output.sigma = sigma;
    output.L_feat = @L_feat;
    output.L_pix = @L_pix;
    output.L = @L;
output._lambda = _lambda;    
output.I_asterisk = I_asterisk;    
output.I_in = I_in;
end

