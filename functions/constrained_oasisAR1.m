function [c, s, Aset] = constrained_oasisAR1(y, g, sn, optimize_b,...
    optimize_g, decimate, maxIter)
%% Infer the most likely discretized spike train underlying an AR(1) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g c_{t-1} >=s_min or =0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities 
        %withone entry per time-bin.
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
        %impulse response.   
%   sn:  scalar, standard deviation of the noise distribution 
%   optimize_b: bool, optimize baseline if True
%   optimize_g: integer, number of large, isolated events to consider for
%       optimizing g 
%   decimate: int, decimation factor for estimating hyper-parameters faster
%       on decimated data 
%   maxIter:  int, maximum number of iterations 
%   Aset: npool x 4 matrix, warm stared active sets 
       
%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes) 
%   b: scalar, fluorescence baseline 
%   g: scalar, parameter of the AR(1) process 
%   lam: scalar, sparsity penalty parameter 
%   Aset: npool x 4 matrix, active sets 

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References 
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

%% initialization
y = reshape(y, [], 1);
if ~exist('g', 'var') || isempty(g);   g = 0.95; end
if ~exist('lam', 'var') || isempty(lam);   lam = 0; end
if ~exist('smin', 'var') || isempty(smin);   smin = 0; end
if ~exist('Aset', 'var') || isempty(Aset)
    len_active_set = length(y);
    Aset = ones(len_active_set, 4); % each row is one pool: (vi, wi, t, l)
    Aset(:,1) = y - lam*(1-g);            % vi
    Aset(:,3) = (1:len_active_set);  % ti
    Aset(end, 1) = y(end) - lam;
    Aset(end, 3) = len_active_set;
else
    len_active_set = size(Aset,1);
end

%% run OASIS
ii = 1;
while ii < len_active_set
    % find the active set
    while (ii<len_active_set) && ...
            (Aset(ii+1,1)/Aset(ii+1,2)>=Aset(ii,1)/Aset(ii,2)*g^(Aset(ii,4))+smin)
        ii = ii + 1;
    end
    
    if ii == len_active_set; break; end
    
    %% merge pools
    Aset(ii,1) = Aset(ii,1) + Aset(ii+1,1)* (g^(Aset(ii,4)));
    Aset(ii,2) = Aset(ii,2) + Aset(ii+1,2)*(g^(2*Aset(ii,4)));
    Aset(ii, 4) = Aset(ii, 4) + Aset(ii+1, 4);
    Aset(ii+1, :) = [];
    len_active_set = len_active_set - 1;
    
    %% backtrack until violations fixed 
    while ii>1 && (Aset(ii,1)/Aset(ii,2)<=Aset(ii-1,1)/Aset(ii-1,2)*g^(Aset(ii-1,4))+smin)
        ii = ii - 1;
        Aset(ii,1) = Aset(ii,1) + Aset(ii+1,1)*(g^(Aset(ii,4)));
        Aset(ii,2) = Aset(ii,2) + Aset(ii+1,2)*(g^(2*Aset(ii,4)));
        Aset(ii, 4) = Aset(ii, 4) + Aset(ii+1, 4);
        Aset(ii+1, :) = [];
        len_active_set = len_active_set - 1;
    end  
end

%% construct solution for all t
c = zeros(size(y));
s = c; 
% first pool 
c(Aset(1,3):Aset(1,4)) = max(0,Aset(1,1)/Aset(1,2)) * (g.^((Aset(1,3):Aset(1,4))-1));
s(1) = 0; 
% other pools 
for ii=2:len_active_set
    t0 = Aset(ii,3); 
    tau = Aset(ii, 4); 
    c(t0:(t0+tau-1)) = max(0,Aset(ii,1)/Aset(ii,2)) * (g.^(0:(tau-1)));
    s(t0) = c(t0) - g*c(t0-1); 
end

