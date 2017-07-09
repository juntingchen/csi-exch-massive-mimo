function [R, Rm] = genChannCorr(Nt, ...         % # of TX antennas
                                theta_deg, ...  % The azimuth angle of the signal (in degree)
                                delta_deg, ...  % Angular spread in degree
                                DB)             % Whether to show test information
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

if strcmp(version, '7.10.0.499 (R2010a)')
    [R, Rm] = genChannCorr_2010a(Nt, ...         % # of TX antennas
                                 theta, ...      % The azimuth angle of the signal (in degree)
                                 delta_deg);     % Angular spread in degree
    return
end

% DEBUG MODE
if nargin < 1
    Nt = 40;
    theta = 0;
    delta_deg = 15;
    DB = 1;
elseif nargin < 4
    DB = 0;
end
% END DEBUG MODE 

theta = theta_deg / 180 * pi;
delta = delta_deg / 180 * pi;
sigma = delta / sqrt(8 * log(2));   % By this, the AS is defined as the 
                                    % width of the beam, for which at the
                                    % boundary, the PAS drops by half,
                                    % i.e., f(delta/2) = f(0) / 2.
%

% c = 1/(sqrt(2 * pi) * sigma * erf(pi/(sqrt(8) * sigma)));
c = 1 / erf(pi / (sqrt(2) * sigma));

Corr = zeros(1, Nt);
for d = 0:Nt - 1
    my_fun = @(x) exp(-(x - theta).^2/(2 * sigma^2) + 1i * pi * d * sin(x)) / sqrt(2 * pi * sigma^2);
    Corr(d + 1) = c * integral(my_fun, theta - pi, theta + pi);
end

R = zeros(Nt);
for p = 1:Nt
    for q = p:Nt
%         R(p, q) = real(Corr(q-p+1));
        R(p, q) = Corr(q-p+1);
        if p ~= q
            R(q, p) = conj(R(p, q));
        end
    end
end
     
Rm = sqrtm(R);
% [U, D] = eig(R);
% d = diag(D);
% d(d < 1e-9) = 0;

if DB
    % Debug ----------------------
    figure,
    stem(cumsum(sort(eig(R), 'descend')));
    title('Eigenvalue distribution');
    DEBUG_MODE_PAUSE_HERE = 1;

    figure,
    my_detpath(R, 2, 8, 1);
    title('Spatial channel response');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------ My subfunctions ------------------------------ % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Aod_id, Path] = my_detpath(R, numPath, numVec, ShowFigure)
% Path detection
% [Aod_id, Path] = detpath(R, numPath, numVec, ShowFigure)

% AS = 15;        % Default angular spread
if nargin < 4
	ShowFigure = 0;
    if nargin < 3
        numVec = 3 * numPath; 
        if nargin < 2
            numPath = 2;
        end
    end
end

n = size(R, 1);
Ut = fft(eye(n));

Ra = Ut' * R * Ut;  % R is the covariance matrix, must real on diagonal
Gain = real(diag(Ra)).';

%%%% REMARK: WE MAY NOT NEED THE AVERAGE
% N = 3;  % Average window
% % N = min(round(AS * n / 180), floor(n/2));
% 
% L_offset = floor(N/2);
% R_offset = - L_offset + N - 1;
% Pwr_3vec = repmat(Gain, 1, 3);
% Pa = zeros(1, n);
% for i = n + 1:n + n
%     Pa(i - n) = sum(Pwr_3vec(i - L_offset:i + R_offset)) / N;
% end

%%%% REMARK: DEFFERENTIATION METHOD
% % dPa = Pa - [Pa(end) Pa(1:end-1)];
% % ddPa = [dPa(2:end) dPa(1)] - dPa;
% % 
% % [~, I_pulse] = sort(ddPa, 'ascend');
% 
% PathPwr = zeros(1, n);
% PathId = zeros(1, n);
% 
% PathPwr(1) = Pa(I_pulse(1));
% PathId(1) = I_pulse(1);
% 
% i_path = 1;
% while i_path < maxPath && Pa(I_pulse(i_path + 1)) > PathPwr(1) * 0.5
%     i_path = i_path + 1;
%     PathPwr(i_path) = Pa(I_pulse(i_path));
%     PathId(i_path) = I_pulse(i_path);
% end
% 
% PathPwr = PathPwr(1:i_path);
% PathId = PathId(1:i_path);

% Average the data / NEEDED?
N = 3;  
L_offset = floor(N/2);
R_offset = - L_offset + N - 1;
Gain_rep = repmat(Gain, 1, 3);
Gain_av = zeros(1, n);
for i = n + 1:n + n
    Gain_av(i - n) = sum(Gain_rep(i - L_offset:i + R_offset)) / N;
end

maxGain = max(Gain_av);
gain_td = maxGain * 0.1;

Gain1 = Gain_av;
Gain1(Gain1 < gain_td) = 0;

Path.pwr = zeros(1, n);
Path.id = zeros(1, n);
Path.Loffset = zeros(1, n);
Path.Roffset = zeros(1, n);

i_path = 0;
while i_path < numPath && ~isempty(find(Gain1 > 0, 1))
    [pathpwr, pathid] = max(Gain1);
    
    Gain1(pathid) = 0;
    % Find the path cluster from the left
    i = pathid - 1;
    L_offset = 0;
    if i < 1
        i = n;
    end
    while Gain1(i) > 0
        L_offset = L_offset + 1;
        Gain1(i) = 0;
        i = i - 1;
        if i < 1
            i = n;
        end
    end
    
    % Find the path cluster from the right
    i = pathid + 1;
    R_offset = 0;
    if i > n
        i = 1;
    end
    while Gain1(i) > 0
        R_offset = R_offset + 1;
        Gain1(i) = 0;
        i = i + 1;
        if i > n
            i = 1;
        end
    end
    
    if R_offset + L_offset > 0
        i_path = i_path + 1;
        Path.pwr(i_path) = pathpwr;
        Path.id(i_path) = pathid;
        Path.Loffset(i_path) = L_offset;
        Path.Roffset(i_path) = R_offset;
    end
        
end
num_found_path = i_path;

Path.pwr = Path.pwr(1:num_found_path);
Path.id = Path.id(1:num_found_path);
Path.Loffset = Path.Loffset(1:num_found_path);
Path.Roffset = Path.Roffset(1:num_found_path);



numVec = max(3, numVec);
Aod_id = zeros(1, numVec);
Gain2 = Gain_av;
cnt = 0;
i = 0;
while cnt < numVec && i < num_found_path
    i = i + 1;
    
    if isempty(find(Aod_id == Path.id(i), 1))
        cnt = cnt + 1;
        Aod_id(cnt) = Path.id(i);
        Gain2(Path.id(i)) = 0;
    end
    
    vecid = Path.id(i) - 1;
    if vecid < 1
        vecid = n;
    end
    if isempty(find(Aod_id == vecid, 1))
        cnt = cnt + 1;
        Aod_id(cnt) = vecid;
        Gain2(vecid) = 0;
    end
    
    vecid = Path.id(i) + 1;
    if vecid > n
        vecid = 1;
    end
    if isempty(find(Aod_id == vecid, 1))
        cnt = cnt + 1;
        Aod_id(cnt) = vecid;
        Gain2(vecid) = 0;
    end
end

while cnt < numVec && ~isempty(find(Gain2 > 0, 1))
    [~, vecid] = max(Gain2);
    cnt = cnt + 1;
    Aod_id(cnt) = vecid;
    Gain2(vecid) = 0;
end
    
Aod_id = Aod_id(1:cnt);

stem(Gain);
hold on
stem(Aod_id, Gain(Aod_id), 'r');
hold off

end

