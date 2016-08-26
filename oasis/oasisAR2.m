function [c, s, Aset] = oasisAR2(y, g1, g2, lam, smin, T_over_ISI, jitter)
%% Infer the most likely discretized spike train underlying an AR(2) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g1*c_{t-1}-g2*c_{t-2} >=s_min or =0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities 
        %withone entry per time-bin.
%   g1:  scalar, first parameter of the AR(2) process that models the fluorescence ...
        %impulse response.   
%   g2:  scalar, second parameter of the AR(2) process that models the fluorescence ...
        %impulse response.   
%   lam:  scalar, sparsity penalty parameter lambda. 
%   smin: scalar, optional, default 0
        %miniumal non-zero activity within each bin (minimal 'spike size').
%   T_over_ISI: scalar, ratio of recording duration T and maximumal
%       inter-spike-interval. default: 1
%   jitter: bool, perform correction step by jittering spike times to
%       minimize RSS. it helps to avoid delayed spike detection. default:
%       false; 

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes) 
%   Aset: npool * 4 matrix, active set 

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References 
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

%% initialization
y = reshape(y, [], 1);
if ~exist('g', 'var') || isempty(g)
    g = estimate_time_constant(y); 
end
if ~exist('lam', 'var') || isempty(lam);   lam = 0; end
if ~exist('smin', 'var') || isempty(smin);   smin = 0; end
if ~exist('T_over_ISI', 'var') || isempty(T_over_ISI)  
    T_over_ISI = 1; 
end
if ~exist('jitter', 'var') || isempty(jitter)
    jitter = false; 
end

%% initialization 
T = length(y); 
yp = y - lam * (1-g1-g2); 
yp(end-1) = y(end-1) - lam*(1-g1); 
yp(end) = y(end) - lam; 

% active set 
len_active_set = length(y);
Aset = [yp, yp, (1:T)', ones(T,1)]; 

% precompute 
len_g = len_active_set / T_over_ISI; 
temp = roots([1, -g1, -g2]); 
d = max(temp); r = min(temp); 
g11 = (exp(log(d)*(1:len_g)') - exp(log(r)*(1:len_g)')) / (d-r); 
g12 = [0; g2*g11(1:(end-1))]; 
g11g11 = cumsum(g11.*g11); 
g11g12 = cumsum(g11.*g12); 

%% run OASIS
ii = 2;
while ii < len_active_set
    % find the active set
    while (ii<len_active_set) && (g11(Aset(ii,4))*Aset(ii,1) + g12(Aset(ii,4))...
            *Aset(ii-1,2) + smin <= Aset(ii+1,1))
        ii = ii + 1;
    end   
    if ii == len_active_set; break; end
    
    %% merge pools
    Aset(ii, 4) = Aset(ii, 4) + Aset(ii+1, 4);
    ti = Aset(ii,3); 
    li = Aset(ii,4); 
    Aset(ii,1) = (g11(1:li)'*yp(ti+(1:li)-1) - g11g12(li)*Aset(ii-1,2))/g11g11(li); 
    Aset(ii,2) = g11(li)*Aset(ii,1) + g12(li)*Aset(ii-1,2); 
    Aset(ii+1, :) = [];
    len_active_set = len_active_set - 1;
    
    %% backtrack until violations fixed 
    while ii>2 &&  (g11(Aset(ii-1,4))*Aset(ii-1,1) + g12(Aset(ii-1,4))...
            *Aset(ii-2,2) + smin > Aset(ii,1))
        ii = ii - 1;
        Aset(ii, 4) = Aset(ii, 4) + Aset(ii+1, 4);
        ti = Aset(ii,3);
        li = Aset(ii,4);
        Aset(ii,1) = (g11(1:li)'*yp(ti+(1:li)-1) - g11g12(li)*Aset(ii-1,2))/g11g11(li);
        Aset(ii,2) = g11(li)*Aset(ii,1) + g12(li)*Aset(ii-1,2);
        Aset(ii+1, :) = [];
        len_active_set = len_active_set - 1;
    end
end

%% jitter 
a_s = Aset; 
if jitter
    pause; 
end

%% construct solution for all t
c = zeros(T,1);
for ii=1:len_active_set
    vi = Aset(ii,1); 
    ti = Aset(ii,3); 
    li = Aset(ii,4); 
    c(ti) = vi; 
    if li>1
        c(ti+1) = g1*c(ti) + g2*a_s(ii-1,2); 
        for m=2:li
            c(ti+m-1) = g1*c(ti+m-2) + g2*c(ti+m-3); 
        end
    end
end 
c(c<0) = 0; 

s = [0;0;0; c(4:end)-g1*c(3:(end-1))-g2*c(2:(end-2))]; 
s(s<smin) = 0; 
