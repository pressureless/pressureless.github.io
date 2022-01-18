function output = kappa(a, t, theta)
% output = kappa(a, t, θ)
%
%    cos, cot, tan, sin, arctan, atan from trigonometry
%    x = a cos(t) cot(t)
%    
%    where 
%    a: ℝ 
%    t: ℝ 
%    
%    y = a cos(t)
%    
%    r = atan(θ)
%    where 
%    θ: ℝ 
%    
%    κ(θ) = (8(3-sin^2(θ))sin^4(θ))/(a(sin^2(2θ)+4)^(3/2)) where θ: ℝ 
%    
%    `$\phi$`(θ) = -arctan(1/2sin(2θ)) where θ: ℝ
%    
    if nargin==0
        warning('generating random input data');
        [a, t, theta] = generateRandomData();
    end
    function [a, t, theta] = generateRandomData()
        a = randn();
        t = randn();
        theta = randn();
    end

    assert(numel(a) == 1);
    assert(numel(t) == 1);
    assert(numel(theta) == 1);

    % x = a cos(t) cot(t)
    x = a * cos(t) * 1./tan(t);
    % y = a cos(t)
    y = a * cos(t);
    % r = atan(θ)
    r = atan(theta);
    function ret = kappa(theta)
        assert(numel(theta) == 1);

        ret = (8 * (3 - sin(theta).^2) * sin(theta).^4) / (a * (sin(2 * theta).^2 + 4).^(3 / 2));
    end

    function ret = phi(theta)
        assert(numel(theta) == 1);

        ret = -atan(1 / 2 * sin(2 * theta));
    end

    output.x = x;
    output.y = y;
    output.r = r;
    output.kappa = @kappa;
    output.phi = @phi;
output.theta = theta;    
output.a = a;
end

