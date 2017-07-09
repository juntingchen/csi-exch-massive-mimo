function xhat = sq_decoder(x, S) 
% 
% Author: Junting Chen (eejtchen@connect.ust.hk)
% Date: Feb 21, 2017.

N = length(x);

y = S.U' * x;
yhat = zeros(N, 1);

m = length(find(S.Bits > 0));

for i = 1:m
    z = y(i);
    zr = real(z);
    zi = imag(z);
    
    I = find(zr > S.Bd{i}, 1, 'last');
    zr_hat = S.Xc{i}(I);
    
    I = find(zi > S.Bd{i}, 1, 'last');
    zi_hat = S.Xc{i}(I);
    
    z_hat = zr_hat + 1i * zi_hat;
    yhat(i) = z_hat;
    
end
xhat = S.U * yhat;
