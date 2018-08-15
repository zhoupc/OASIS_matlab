oasis_folder = fileparts(mfilename('fullpath'));

if ~exist('oasis_loaded', 'var') || ~oasis_loaded
    addpath(oasis_folder);
    addpath(fullfile(oasis_folder, 'functions'));
    addpath(fullfile(oasis_folder, 'packages', 'oasis'));
    addpath(fullfile(oasis_folder, 'packages', 'oasis_kernel'));
    addpath(fullfile(oasis_folder, 'packages', 'constrained-foopsi'));
    addpath(fullfile(oasis_folder, 'packages', 'MCMC'));
    addpath(fullfile(oasis_folder, 'packages', 'MCMC', 'utilities'));
    oasis_loaded = true;
end

%% install convex optimization solvers
% install cvx
if isempty(which('cvx_begin.m'))
    if ~exist(fullfile(oasis_folder, 'packages', 'cvx'), 'dir')
        %install cvx
        if ismac
            cvx_url = 'http://web.cvxr.com/cvx/cvx-maci64.zip';
        elseif isunix
            cvx_url = 'http://web.cvxr.com/cvx/cvx-a64.zip';
        elseif ispc
            cvx_url = 'http://web.cvxr.com/cvx/cvx-w64.zip';
        else
            fprints('Your platform is not supported by CVX\n');
            return;
        end
        fprintf('Downloading CVX...\n');
        unzip(cvx_url, fullfile(oasis_folder, 'packages'));
    end
    run(fullfile(oasis_folder, 'packages', 'cvx', 'cvx_setup.m'));
end
%% save the current path
%savepath();
