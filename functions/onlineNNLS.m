function [c, s] = onlineNNLS(y, g, lam, mid, post, tol, maxIter)
%% Infer the most likely discretized spike train underlying an AR(2) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g c_{t-1} >=s_min or =0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities 
        %withone entry per time-bin.
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
        %impulse response.   
%   lam:  scalar, sparsity penalty parameter lambda. 

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes) 

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References 
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

%% input arguments  
if ~exist('lam', 'var') || isempty(lam)
    lam = 0; 
end
if ~exist('mid', 'var') || isempty(mid)
    mid = 100; 
end
if ~exist('post', 'var') || isempty(post)
    post = 100; 
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-9; 
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = []; 
end

%% get the response kernel of AR(2) model 
temp = roots([1, -g(1), -g(2)]); 
d = max(temp); 
r = min(temp); 
w = mid + post;     % window size 
h = (exp(log(d)*(1:w)) - exp(log(r)*(1:w))) / (d-r); % convolution kernel 
[u, t] = meshgrid(1:w, 1:w);  
K = zeros(w); 
ind = 1+t-u; 
K(ind>0) = h(ind(ind>0));   % convolution matrix 
A = K'*K; 

%% initialization 
T = length(y); 
y = reshape(y, [], 1); 
yp = y - lam*(1-sum(g)); 
yp(end-1) = y(end-1) - lam*(1-g(1));
yp(end) = y(end) - lam; 
s = zeros(T, 1); 
c = zeros(T, 1); 

%% run online deconvolution 
t = 1; 
while t <= T-w+1
    ind = t:(t+w-1); 
    s(ind) = nnls(A, K'*yp(ind), s(ind), tol, maxIter); 
    yp(ind) = yp(ind) - K(:, 1:mid)*s(t:(t+mid-1)); 
    c(ind) = c(ind) + K(:, 1:mid)*s(t:(t+mid-1)); 
    t = t + mid; 
end 
s(t:T) = nnls(A((t+w-T):w, (t+w-T):w), K(1:(T-t+1), 1:(T-t+1))'*yp(t:T), s(t:T), tol, maxIter); 
c(t:T) = c(t:T) + K((t+w-T):w, (t+w-T):w) * s(t:T); 

function s = nnls(A, b, s, tol, maxIter)
%% fast algorithm for solving nonnegativity constrained least squared
% problem minize norm(y-K*s, 2), s.t. s>=0. 

%% inputs: 
%   A: n x p matrix, K'*K
%   b: n x 1 vector, K'*y
%   s: p x 1 vector, warm started s 
%   tol: scalar, smallest nonzero values 
%   maxIter: scalar, maximum nonzero values 

%% outputs: 
%   s: p x 1 vector, solution 

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References 
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging
% Bro R & Jong S, Journal of Chemometrics 1997, A FAST NON-NEGATIVITY-CONSTRAINED LEAST SQUARES ALGORITHM


%% input arguments 
p = size(A,2);      % dimension of s 
if ~exist('s', 'var') || isempty(s)
    s = zeros(p, 1); 
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-9; 
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = p; 
end

for miter=1:maxIter
l = b - A*s;            % negative gradient
    Pset = (s>0);       % passive set

    if max(l) < tol         % no more passive set 
        break; 
    end
    
    [~, temp] = max(l);         % choose the one with the largest gradient
    Pset(temp) = true;         % add it to the passive set
    
    % correct nonnegativity violations
    while any(Pset)
        % run unconstrained least squares for variables in passive sets
        try
            mu = A(Pset, Pset) \ b(Pset);
        catch
            mu = (A(Pset, Pset) + tol*eye(sum(Pset))) \ b(Pset);
        end
                
        if all(mu>tol)
            break; 
        end
        
        temp = s(Pset) ./ (s(Pset)-mu);
        temp(mu>tol) = [];
        a = min(temp);
        s(Pset) = s(Pset) + a*(mu-s(Pset));
        Pset(s<=tol) = false;
    end 
    
    s(Pset) = mu;
end
