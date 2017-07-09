function [Bd, Xc, H, D, delta, alpha, Pcell] = unifsq_gauss(R, sigma2)
% [Bd, Xc, H, D, delta, alpha, Pcell] = unifsq_gauss(R, sigma2)
%
% Entropy constrained uniform (infinite-level) scalar quantizer
% 
% INPUT
%   sigma2          Variance of the Gaussian input source
%   R               Required output entropy
%
% OUTPUT
%   Bd              Array of boundary points
%   Xc              Array of reconstruction points
%   Hp              Output entropy
%
% Author: Junting Chen (eejtchen@connect.ust.hk)
% Date: Feb 21, 2017.

sigma = sqrt(sigma2);

if R > 0.45
    delta_max = 4 * sigma;
    delta_min = 0;
else
    P = 0.01:0.01:0.5;
    Hb = - P .* log2(P) - (1 - P) .* log2(1 - P);
    [~, I] = min(abs(R - Hb));
    pout = P(I);
    delta_max = 4 * qfuncinv(pout) * sigma;
    delta_min = qfuncinv(pout) * sigma;
end


MAXLOOP = 100;
eps = 0.02;

i = 0;
H = -1;
while (R - H > eps || R < H) && i < MAXLOOP
    i = i + 1;
    
    delta = (delta_max + delta_min) / 2;
    alpha = 0;
    
    lambda = delta / sigma;
    
    Kmax = ceil(3 * sigma / delta);     % suffcieint to calcualte the [-3 sigma, 3 sigma] region
    Pcell = zeros(1, 2 * Kmax + 1);
    Pcell(1) = qfunc(- alpha * lambda) - qfunc((1 - alpha) * lambda);
    for k = 1:Kmax - 1
        Pcell(2 * k) = qfunc((k - alpha) * lambda) - qfunc((k + 1 - alpha) * lambda);
        Pcell(2 * k + 1) = qfunc((- k - alpha) * lambda) - qfunc((- k + 1 - alpha) * lambda);
    end
    k = Kmax;
    Pcell(2 * k) = qfunc((k - alpha) * lambda) - 0;
    Pcell(2 * k + 1) = 1 - qfunc((- k + 1 - alpha) * lambda);
    
    H = - sum(Pcell .* log2(Pcell));
    
    if H > R %&& R > 2
        H0 = H;
        
        alpha_max = 1 / 2;
        alpha_eps = alpha_max / 100;
        alpha_st = alpha_max / 10;
        
        alpha0 = alpha;
        
        while alpha_st > alpha_eps && alpha < 1 / 2 
            alpha = alpha0 + alpha_st;

            Pcell(1) = qfunc(- alpha * lambda) - qfunc((1 - alpha) * lambda);
            for k = 1:Kmax - 1
                Pcell(2 * k) = qfunc((k - alpha) * lambda) - qfunc((k + 1 - alpha) * lambda);
                Pcell(2 * k + 1) = qfunc((- k - alpha) * lambda) - qfunc((- k + 1 - alpha) * lambda);
            end
            k = Kmax;
            Pcell(2 * k) = qfunc((k - alpha) * lambda) - 0;
            Pcell(2 * k + 1) = 1 - qfunc((- k + 1 - alpha) * lambda);
            H1 = - sum(Pcell .* log2(Pcell));

            if H1 > H0 || H1 < R
                alpha_st = alpha_st / 2;
            else
                alpha0 = alpha;
                H0 = H1;
            end
        end
        if H1 > H0 && H1 > R
            H = H0;
            alpha = alpha0;
        else
            H = H1;
        end
    end
    
    if H > R
        delta_min = delta;
    elseif H < R
        delta_max = delta;
        
        % delta1 = delta;
        % alpha1 = alpha;
    end
end

% -> output alpha and delta

Xc0 = zeros(1, Kmax * 2 + 1);

lambda = delta / sigma;
Pcell(1) = qfunc(- alpha * lambda) - qfunc((1 - alpha) * lambda);
for k = 1:Kmax - 1
    Pcell(2 * k) = qfunc((k - alpha) * lambda) - qfunc((k + 1 - alpha) * lambda);
    Pcell(2 * k + 1) = qfunc((- k - alpha) * lambda) - qfunc((- k + 1 - alpha) * lambda);
end
k = Kmax;
Pcell(2 * k) = qfunc((k - alpha) * lambda) - 0;
Pcell(2 * k + 1) = 1 - qfunc((- k + 1 - alpha) * lambda);

Xc0(1) = sigma / sqrt(2 * pi) * (exp( - ( - alpha)^2 * lambda^2 / 2) ...
                        - exp( - (1 - alpha)^2 * lambda^2 / 2)) / Pcell(1);
for k = 1:Kmax - 1
    Xc0(2 * k) = sigma / sqrt(2 * pi) * (exp( - (k - alpha)^2 * lambda^2 / 2) ...
                            - exp( - (k + 1 - alpha)^2 * lambda^2 / 2)) / Pcell(2 * k);
    Xc0(2 * k + 1) = sigma / sqrt(2 * pi) * (exp( - (- k - alpha)^2 * lambda^2 / 2) ...
                            - exp( - (- k + 1 - alpha)^2 * lambda^2 / 2)) / Pcell(2 * k + 1);
end
k = Kmax;
Xc0(2 * k) = sigma / sqrt(2 * pi) * (exp( - (k - alpha)^2 * lambda^2 / 2) ...
                        - 0) / Pcell(2 * k);
Xc0(2 * k + 1) = sigma / sqrt(2 * pi) * (0 ...
                        - exp( - (- k + 1 - alpha)^2 * lambda^2 / 2)) / Pcell(2 * k + 1);
%  
% Boundary and reconstruction points of the quantizer
Xc = sort(Xc0, 'ascend');    
Bd = ((- Kmax:Kmax) - alpha) * delta;
Bd(1) = -inf;

D = 0;
% Distortion
K = length(Xc);
d = zeros(1, K);
for k = 1:K
    a = Bd(k);
    if k < K
        b = Bd(k + 1);
    else
        b = Inf;
    end
    d(k) = integral(@(x) (x - Xc(k)).^2 .* exp(- x.^2 / (2 * sigma2)) / (sqrt(2 * pi) * sigma), ...
        a, b);
end

D = sum(d);