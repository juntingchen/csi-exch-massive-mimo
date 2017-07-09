function [R, Rm] = genChannCorr_2010a(Nt, ...         % # of TX antennas
                                theta, ...      % The azimuth angle of the signal (in degree)
                                delta_deg, ...  % Angular spread in degree
                                rk)             % Expected rank of the channel
%
% Description:
%   Generate the channel spatial correlaiton on the TX side. The pathloss
%   is assumed to be always normalized to 1.
%
% MODEL:
%   Truncated Gaussian distribution for power angular spectrum density.
%
% OUTPUT
%   R              Nt * Nt matrix
%   Rm             Square root of R
%

% if strcmp(version, '7.10.0.499 (R2010a)')
%     return
% end

% DEBUG MODE
if nargin < 1
    Nt = 40;
    theta = 0;
    delta_deg = 15;
    rk = 10;
end
% END DEBUG MODE 

delta = delta_deg / 180 * pi;
sigma = delta / sqrt(8 * log(2));   % By this, the AS is defined as the 
                                    % width of the beam, for which at the
                                    % boundary, the PAS drops by half,
                                    % i.e., f(delta/2) = f(0) / 2.
%

c = 1/(sqrt(2 * pi) * sigma * erf(pi/(sqrt(8) * sigma)));

Corr = zeros(1, Nt);
for d = 0:Nt - 1
    my_fun = @(x) exp(-(x - theta).^2/(2 * sigma^2) + 1i * pi * d * sin(x));
    Corr(d + 1) = c * int(my_fun, -pi/2, pi/2);
end

R = zeros(Nt);
for p = 1:Nt
    for q = p:Nt
        R(p, q) = real(Corr(q-p+1));
        if p ~= q
            R(q, p) = R(p, q);
        end
    end
end
     
Rm = sqrtm(R);
% [U, D] = eig(R);
% d = diag(D);
% d(d < 1e-9) = 0;

if nargin < 1
    % Debug ----------------------
    figure,
    stem(sort(eig(R), 'descend'));
    DEBUG_MODE_PAUSE_HERE = 1;
end
