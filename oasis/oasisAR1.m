function [c, s, active_set] = oasisAR1(y, g, lam, smin, active_set)
%% Infer the most likely discretized spike train underlying an AR(1) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g c_{t-1} >=s_min or =0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities 
        %withone entry per time-bin.
      % OR %%
      % len_active_set*4 matrix, active set 
      
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
        %impulse response.   
%   lam:  scalar, sparsity penalty parameter lambda. 
%   smin: scalar, optional, default 0
        %miniumal non-zero activity within each bin (minimal 'spike size').
%   active_set: npool x 4 matrix, warm stared active sets 
       
%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes) 
%   active_set: npool x 4 matrix, active sets 

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
if ~exist('active_set', 'var') || isempty(active_set)
    len_active_set = length(y);
    active_set = ones(len_active_set, 4); % each row is one pool: (vi, wi, t, l)
    active_set(:,1) = y - lam*(1-g);            % vi
    active_set(:,3) = (1:len_active_set);  % ti
    active_set(end, 1) = y(end) - lam;
    active_set(end, 3) = len_active_set;
else
    len_active_set = size(active_set,1);
end

%% run OASIS
ii = 1;
while ii < len_active_set
    % find the active set
    while (ii<len_active_set) && ...
            (active_set(ii+1,1)/active_set(ii+1,2)>=active_set(ii,1)/active_set(ii,2)*g^(active_set(ii,4))+smin)
        ii = ii + 1;
    end
    
    if ii == len_active_set; break; end
    
    %% merge pools
    active_set(ii,1) = active_set(ii,1) + active_set(ii+1,1)* (g^(active_set(ii,4)));
    active_set(ii,2) = active_set(ii,2) + active_set(ii+1,2)*(g^(2*active_set(ii,4)));
    active_set(ii, 4) = active_set(ii, 4) + active_set(ii+1, 4);
    active_set(ii+1, :) = [];
    len_active_set = len_active_set - 1;
    
    %% backtrack until violations fixed 
    while ii>1 && (active_set(ii,1)/active_set(ii,2)<=active_set(ii-1,1)/active_set(ii-1,2)*g^(active_set(ii-1,4))+smin)
        ii = ii - 1;
        active_set(ii,1) = active_set(ii,1) + active_set(ii+1,1)*(g^(active_set(ii,4)));
        active_set(ii,2) = active_set(ii,2) + active_set(ii+1,2)*(g^(2*active_set(ii,4)));
        active_set(ii, 4) = active_set(ii, 4) + active_set(ii+1, 4);
        active_set(ii+1, :) = [];
        len_active_set = len_active_set - 1;
    end  
end

%% construct solution for all t
c = zeros(size(y));
s = c; 
for ii=1:len_active_set
    t0 = active_set(ii,3); 
    tau = active_set(ii, 4); 
    c(t0:(t0+tau-1)) = max(0,active_set(ii,1)/active_set(ii,2)) * (g.^(0:(tau-1)));
end

s(active_set(2:end,3)) = c(active_set(2:end,3)) - g*c(active_set(2:end,3)-1);
