function [c, s, b, g, lam, active_set] = constrained_oasisAR2(y, g, sn,...
    optimize_b, optimize_g, T_over_ISI, maxIter, smin)
%% Infer the most likely discretized spike train underlying an AR(1) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g c_{t-1} >=s_min or =0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   g1:  scalar, the first parameter of the AR(2) process that models the fluorescence ...
%impulse response.%
%   g2:  scalar, the second parameter of the AR(2) process that models the fluorescence ...
%impulse response.
%   sn:  scalar, standard deviation of the noise distribution
%   optimize_b: bool, optimize baseline if True
%   optimize_g: integer, number of large, isolated events to consider for
%       optimizing g
%   T_over_ISI: int, T/ISI
%   maxIter:  int, maximum number of iterations
%   smin: scalar, minimize spike size

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes)
%   b: scalar, fluorescence baseline
%   g: scalar, parameter of the AR(1) process
%   lam: scalar, sparsity penalty parameter
%   active_set: npool x 4 matrix, active sets

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging


%% input arguments
y = reshape(y, [], 1);
T = length(y);

if ~exist('g', 'var') || isempty(g)
    g = estimate_time_constant(y, 2);
end
if ~exist('sn', 'var') || isempty(sn)
    sn = GetSn(y);
end
if ~exist('lam', 'var') || isempty(lam);   lam = 0; end
if ~exist('optimize_b', 'var') || isempty(optimize_b)
    optimize_b = false;
end
if ~exist('optimize_g', 'var') || isempty(optimize_g)
    optimize_g = 0;
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 10;
end

thresh = sn * sn * T;
lam = 0;

g_converged = false;

%% initialize the active_set
b = 0;
[solution, spks, active_set] = oasisAR2(y, g, lam);
len_active_set = size(active_set,1);

if optimize_b
    % update b and lam
    b_new = mean(y-solution);
    db = b_new - b;
    lam = -db / (1-g);
    b = b_new;
    % correct the last pool
    active_set(end,1) = active_set(end,1) - lam*g^(active_set(end,4));
    ti = active_set(end,3); li = active_set(end,4); idx = 0:(li-1);
    solution(ti+idx) = max(0, active_set(end,1)/active_set(end,2)) * (g.^idx);
end

%% optimize parameters
tol = 1e-4;
flag_lam = true;
if ~optimize_b
    %% don't optimize b nor g, just the dual variable lambda
    for miter=1:maxIter
        res = y - solution;
        RSS = res' * res;
        if or(abs(RSS-thresh) < tol, ~flag_lam)  % constrained form has been found, stop
            break;
        else
            % update lam
            [solution, active_set, lam, spks, flag_lam] = update_lam(y, solution, ...
                active_set, g, lam, thresh);
            
            % update g
            if and(optimize_g, ~g_converged);
                g0 = g;
                [solution, active_set, g, spks] = update_g(y, active_set, g, lam);
                if abs(g-g0)/g0 < 1e-4; g_converged = true; end
            end
        end
    end
else
    %% optimize the baseline b and dependends on the optimized g too
    % run iterations until the noise constraint is tight or spike train is
    % empty or maxIter reached.
    for miter=1:maxIter
        res = y - solution;
        RSS = res' * res;
        if or(abs(RSS-thresh) < tol, sum(solution)<1e-9)
            break;
        end
        
        % update b and g
        update_b_lam();
        
        % update b and g
        if and(optimize_g, ~g_converged);
            g0 = g;
            update_b_g();
            if abs(g-g0)/g0 < 1e-4;
                g_converged = true;
            end
        end
    end
end

c = solution;
s = spks;

%% nested functions
    function update_b_lam()
        %% update total shift dphi due to contribution of baseline and labda
        temp = zeros(size(solution));
        maxl = max(active_set(:, 4));
        h = g.^(0:maxl);
        tmp_sum = 0;
        for ii=1:len_active_set
            ti = active_set(ii, 3);
            li = active_set(ii, 4);
            idx = 0:(li-1);
            if ii<len_active_set
                temp(ti+idx) = (1-g^li)/ active_set(ii,2) * h(1:li);
                tmp_sum = tmp_sum + (1-g^li).^2 / active_set(ii,2);
            else
                temp(ti+idx) = 1/active_set(ii,2) * h(1:li);
                tmp_sum = tmp_sum + 1.0 / active_set(ii,2);
            end
        end
        temp = temp - tmp_sum/T/((1-g));
        
        res = y - solution - mean(y-solution);
        aa = temp'*temp;
        bb = res'*temp;
        cc = RSS-thresh;
        dphi = (-bb + sqrt(bb^2-aa*cc)) / aa;
        if imag(dphi)~=0
            return;
        end
        b = b + dphi * (1-g);
        
        % perform shift
        active_set(:,1) = active_set(:,1) - dphi*(1-g.^active_set(:,4));
        
        % run OASIS
        [solution, spks, active_set] = oasisAR1(y, g, [], [], active_set);
        len_active_set = size(active_set,1);
        
        % update b and lam
        b_new = mean(y-solution);
        db = b_new - b;
        dlam = -db / (1-g);
        lam = lam + dlam;
        b = b_new;
        
        %% correct the last pool
        % correct the last pool
        active_set(end,1) = active_set(end,1) - lam*g^(active_set(end,4));
        ti = active_set(end,3); li = active_set(end,4); idx = 0:(li-1);
        solution(ti+idx) = max(0, active_set(end,1)/active_set(end,2)) * (g.^idx);
        len_active_set = size(active_set,1);
    end

    function update_b_g()
        %% update  b and g
        %% initialization
        maxl = max(active_set(:, 4));   % maximum ISI
        c = zeros(size(y));     % the optimal denoised trace
        
        %% find the optimal g and get the warm started active_set
        g = fminbnd(@rss_b_g, 0, 1);
        h = exp(log(g)*(0:maxl)');   % response kernel
        hs = cumsum(h);
        hh = cumsum(h.*h);        % hh(k) = h(1:k)'*h(1:k)
        yh = zeros(len_active_set,1);
        for ii=1:len_active_set
            li = active_set(ii, 4);
            ti = active_set(ii, 3);
            idx = ti:(ti+li-1);
            if ii<len_active_set
                yh(ii) = (y(idx)-lam*(1-g))' * h(1:li);
            else
                yh(ii) = (y(idx)-lam)' * h(1:li);
            end
        end
        aa = hs(active_set(:, 4));
        bb = hh(active_set(:, 4));
        b = (sum(yh) - sum(yh.*aa./bb)) / (length(y)-sum(aa.*aa./bb));
        active_set(:,1) = (yh-b*aa)./bb;
        
        [c,s,active_set] = oasisAR1(y, g, lam, [], active_set);
        len_active_set = size(active_set,1);
        
        function rss = rss_b_g(g)
            h = exp(log(g)*(0:maxl)');   % response kernel
            hs = cumsum(h);
            hh = cumsum(h.*h);        % hh(k) = h(1:k)'*h(1:k)
            yh = zeros(len_active_set,1);
            for iii=1:len_active_set
                li = active_set(iii, 4);
                ti = active_set(iii, 3);
                idx = ti:(ti+li-1);
                
                if iii<len_active_set
                    yh(iii) = (y(idx)-lam*(1-g))' * h(1:li);
                else
                    yh(iii) = (y(idx)-lam)' * h(1:li);
                end
            end
            aa = hs(active_set(:, 4));
            bb = hh(active_set(:, 4));
            tmp_b = (sum(yh) - sum(yh.*aa./bb)) / (length(y)-sum(aa.*aa./bb));
            v = max(0, (yh-tmp_b*aa)./bb);
            for jj=1:len_active_set
                li = active_set(jj, 4);
                ti = active_set(jj, 3);
                idx = ti:(ti+li-1);
                c(idx) = v(jj) * h(1:li);
            end
            res = y-c;
            rss = res'*res;     % residual sum of squares
        end
        
    end
end





