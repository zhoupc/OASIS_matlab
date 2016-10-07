%% test modulate for all oasis functions. 
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]}; % colors

%% example 1:  foopsi, AR1 model. This model is used when the sampling rate is low
g = 0.95;         % AR coefficient 
noise = .3; 
T = 3000; 
framerate = 30;     
firerate = 0.5; 
b = 0;              % baseline 
N = 1;              % number of trials 
seed = 13;          % seed for genrating random variables 
[y, true_c, true_s] = gen_data(g, noise, T, framerate, firerate, b, N, seed); 

% case 1: all parameters are known 
lambda = 2.4; 
[c_oasis, s_oasis] = deconvolveCa(y, 'ar1', g, 'foopsi', 'lambda', lambda);  %#ok<*ASGLU>
[c_cvx, s_cvx] = foopsi(y, g, lambda); 

figure('name', 'FOOPSI, AR1, known: g, lambda', 'papersize', [15, 4]); 
plot_cvx = true; 
show_results; 
plot_cvx = false; 

% case 2: know lambda
lambda = 2.4; 
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', 'foopsi', 'lambda', lambda); 

fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 

figure('name', 'FOOPSI, AR1, known:lambda, estimated: g', 'papersize', [15, 4]); 
show_results; 

% case 3: know lambda, fit g
lambda = 2.4; 
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', 'foopsi', 'lambda', lambda, ...
    'optimize_pars'); 

fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 

figure('name', 'FOOPSI, AR1, known:lambda, estimated:b, g, updated: b, g', 'papersize', [15, 4]); 
show_results; 
%%%%%%%%%%%%%%  END %%%%%%%%%%%%%%%%%%

%% example 2: foopsi, AR2 model 
g = [1.7, -0.712];         % AR coefficient 
noise = 1; 
T = 3000; 
framerate = 30;     
firerate = 0.5; 
b = 0;              % baseline 
N = 20;              % number of trials 
seed = 3;          % seed for genrating random variables 
[Y, trueC, trueS] = gen_data(g, noise, T, framerate, firerate, b, N, seed); 
y = Y(1,:); 
true_c = trueC(1,:);  %#ok<*NASGU>
true_s = trueS(1,:); 
% case 1: all parameters are known 
lambda = 25; 
[c_oasis, s_oasis] = deconvolveCa(y, 'ar2', g, 'foopsi', 'lambda', lambda);  %#ok<*ASGLU>

figure('name', 'FOOPSI, AR2, known: g, lambda', 'papersize', [15, 4]); 
show_results; 

% case 2: know lambda
lambda = 2.5; 
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar2', 'sn', noise, 'foopsi', 'lambda',...
    lambda); 
fprintf('true gamma:        %.3f\t %.3f\n', g(1), g(2)); 
fprintf('estimated gamma:   %.3f\t %.3f\n', options.pars(1),  options.pars(2)); 

figure('name', 'FOOPSI, AR2, known:lambda, estimated: g', 'papersize', [15, 4]); 
show_results; 

%%%%%%%%%%%%%%  END %%%%%%%%%%%%%%%%%%

%% example 3: foopsi, convolution kernel 
g = [1.7, -0.712];         % AR coefficient 
noise = 1; 
T = 3000; 
framerate = 30;     
firerate = 0.5; 
b = 0;              % baseline 
N = 20;              % number of trials 
seed = 3;          % seed for genrating random variables 
[Y, trueC, trueS] = gen_data(g, noise, T, framerate, firerate, b, N, seed); 
y = Y(1,:); 
true_c = trueC(1,:);  %#ok<*NASGU>
true_s = trueS(1,:);
temp = roots([1, -g(1), -g(2)]);
d = max(temp);
r = min(temp);
w = 200;
ht = (exp(log(d)*(1:w)) - exp(log(r)*(1:w))) / (d-r); % convolution kernel
  
% case 1: use the difference of two exponential functions to construct a
% kernel 
lambda = 25; 
pars = [d, r]; 
[c_oasis, s_oasis] = deconvolveCa(y, 'exp2', pars, 'foopsi', 'lambda', lambda, ...
    'shift', 100, 'window', 200);  %#ok<*ASGLU>

figure('name', 'FOOPSI, exp2, known: g, lambda', 'papersize', [15, 4]); 
show_results; 

% case 2: use the kernel directly 
lambda = 25; 
[c_oasis, s_oasis] = deconvolveCa(y, 'kernel', ht, 'foopsi', 'lambda', ...
    lambda, 'shift', 100, 'window', 200);  %#ok<*ASGLU>

figure('name', 'FOOPSI, kernel, known: g, lambda', 'papersize', [15, 4]); 
show_results; %%%%%%%%%%%%%%  END %%%%%%%%%%%%%%%%%%

%% example 4: hard threshold, AR1 model
g = 0.95;         % AR coefficient 
noise = .3; 
T = 3000; 
framerate = 30;     
firerate = 0.5; 
b = 0;              % baseline 
N = 1;              % number of trials 
seed = 13;          % seed for genrating random variables 
[y, true_c, true_s] = gen_data(g, noise, T, framerate, firerate, b, N, seed); 

% case 1: all parameters are known 
smin = 0.5; 
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', g, 'thresholded', 'smin', smin);  %#ok<*ASGLU>

figure('name', 'threshold, AR1, known: g, lambda, smin', 'papersize', [15, 4]); 
show_results; 

% case 2: know lambda
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', 'thresholded', 'smin', smin); 

fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 

figure('name', 'threshold, AR1, known:lambda, estimated: g', 'papersize', [15, 4]); 
show_results; 

% case 5: optimize the thershold, g, and the baseline
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', g,  ...
    'thresholded', 'optimize_smin', 'optimize_pars', 'thresh_factor', 0.99);  %#ok<*ASGLU>

figure('name', 'threshold, AR1, known: g, sn, estimate: smin', 'papersize', [15, 4]); 
show_results; 

% case 6: optimize the thershold, g, and the baseline
[c_oasis, s_oasis, options] = deconvolveCa(y+1, 'ar1', g,  ...
    'thresholded', 'optimize_smin', 'optimize_pars', 'thresh_factor', .97, ...
    'optimize_b');  %#ok<*ASGLU>

figure('name', 'threshold, AR1, known: g, sn, estimate: smin', 'papersize', [15, 4]); 
show_results; 


%%%%%%%%%%%%%%  END %%%%%%%%%%%%%%%%%%


%% example 5: threshold, AR2 model 
g = [1.7, -0.712];         % AR coefficient 
noise = 1; 
T = 3000; 
framerate = 30;     
firerate = 0.5; 
b = 0;              % baseline 
N = 20;              % number of trials 
seed = 3;          % seed for genrating random variables 
[Y, trueC, trueS] = gen_data(g, noise, T, framerate, firerate, b, N, seed); 
y = Y(1,:); 
true_c = trueC(1,:);  %#ok<*NASGU>
true_s = trueS(1,:); 
% case 1: all parameters are known 
smin = 0.5; 
[c_oasis, s_oasis] = deconvolveCa(y, 'ar2', g, 'thresholded', 'smin', smin);  %#ok<*ASGLU>
figure('name', 'threshold, AR2, known: g, smin', 'papersize', [15, 4]); 
show_results; 

% case 2: know smin
smin = 0.5; 
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar2', 'sn', noise, 'thresholded',...
    'smin', smin); 
fprintf('true gamma:        %.3f\t %.3f\n', g(1), g(2)); 
fprintf('estimated gamma:   %.3f\t %.3f\n', options.pars(1),  options.pars(2)); 

figure('name', 'threshold, AR2, known:smin, estimated: g', 'papersize', [15, 4]); 
show_results; 

%% case 3: estimate smin 
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar2', 'sn', noise, 'thresholded',...
    'optimize_smin','optimize_pars', 'thresh_factor', 1); 
% fprintf('true gamma:        %.3f\t %.3f\n', g(1), g(2)); 
% fprintf('estimated gamma:   %.3f\t %.3f\n', options.pars(1),  options.pars(2)); 
fprintf('estimated smin:    %.3f\n', options.smin); 
figure('name', 'threshold, AR2, known:smin, estimated: g', 'papersize', [15, 4]); 
show_results; 
%%%%%%%%%%%%%%  END %%%%%%%%%%%%%%%%%%

%% example 6: threshold, convolution kernel  
g = [1.7, -0.712];         % AR coefficient 
noise = 1; 
T = 3000; 
framerate = 30;     
firerate = 0.5; 
b = 0;              % baseline 
N = 20;              % number of trials 
seed = 3;          % seed for genrating random variables 
[Y, trueC, trueS] = gen_data(g, noise, T, framerate, firerate, b, N, seed); 
y = Y(1,:); 
true_c = trueC(1,:);  %#ok<*NASGU>
true_s = trueS(1,:); 
temp = roots([1, -g(1), -g(2)]);
d = max(temp); 
r = min(temp);
w = 200;
ht = (exp(log(d)*(1:w)) - exp(log(r)*(1:w))) / (d-r); % convolution kernel
  
% case 1: all parameters are known 
smin = 0.5; 
pars = [d, r]; 
[c_oasis, s_oasis] = deconvolveCa(y, 'exp2', pars, 'thresholded', 'smin', smin);  %#ok<*ASGLU>
figure('name', 'threshold, exp2, known: taur, taud, smin', 'papersize', [15, 4]); 
show_results; 

% case 1: all parameters are known 
smin = 0.5; 
[c_oasis, s_oasis] = deconvolveCa(y, 'kernel', ht, 'thresholded', 'smin', smin);  %#ok<*ASGLU>
figure('name', 'threshold, kernel, known: kernel, smin', 'papersize', [15, 4]); 
show_results; 
%%%%%%%%%%%%%%  END %%%%%%%%%%%%%%%%%%

%% example 7:  constrained-foopsi, AR1
g = 0.95;         % AR coefficient 
noise = .3; 
T = 3000; 
framerate = 30;     
firerate = 0.5; 
b = 0;              % baseline 
N = 1;              % number of trials 
seed = 13;          % seed for genrating random variables 
[y, true_c, true_s] = gen_data(g, noise, T, framerate, firerate, b, N, seed); 

% cvx solution 
[c_cvx, s_cvx] = constrained_foopsi(y, g, noise); 
% case 1: all parameters are known 
[c_oasis, s_oasis] = deconvolveCa(y, 'ar1', g, 'constrained', 'sn', noise);  %#ok<*ASGLU>

figure('name', 'constrained-FOOPSI, AR1, known: g, sn', 'papersize', [15, 4]);
plot_cvx = true; 
show_results; 
plot_cvx = false; 

% case 2: nothing is known, estimate g with auto-correlation method
[c_oasis, s_oasis,options] = deconvolveCa(y, 'ar1', 'constrained'); 

fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 

figure('name', 'FOOPSI, AR1, estimated: g, sn', 'papersize', [15, 4]); 
show_results; 

% case 3: nothing is know, estimate g with auto-correlation method first
% and then update it to minimize the RSS
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', 'constrained', ...
    'optimize_pars'); 

fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 

figure('name', 'FOOPSI, AR1, estimated: g, sn, update:g', 'papersize', [15, 4]); 
show_results; 

% case 4: nothing is know, estimate g with auto-correlation method first
% and then update it to minimize the RSS, the baseline is also unknown
true_b = 0.5; 
[c_oasis, s_oasis, options] = deconvolveCa(y+true_b, 'ar1', g,...
    'constrained','optimize_b', 'sn', noise); 
fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 
fprintf('true b:       %.3f\n', true_b); 
fprintf('estimated b:       %.3f\n', options.b); 
fprintf('tuning parameter:  %.3f\n', options.lambda); 

figure('name', 'FOOPSI, AR1, estimated: g, sn, lambda', 'papersize', [15, 4]); 
show_results; 

% case 5: nothing is know, estimate g with auto-correlation method first
% and then update it to minimize the RSS, the baseline is also unknown
true_b = 0.5; 
[c_oasis, s_oasis, options] = deconvolveCa(y+true_b, 'ar1',...
    'constrained','optimize_b', 'optimize_pars'); 
fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 
fprintf('estimated b:       %.3f\n', options.b); 
fprintf('tuning parameter:  %.3f\n', options.lambda); 

figure('name', 'FOOPSI, AR1, estimated: g, sn, lambda, update:g', 'papersize', [15, 4]); 
show_results; 

%% 