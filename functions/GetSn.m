function sn = GetSn(Y, range_ff, method)
%% Estimate noise standard deviation

%% inputs:
%   Y: N X T matrix, fluorescence trace
%   range_ff : 1 x 2 vector, nonnegative, max value <= 0.5, range of frequency (x Nyquist rate) over which the spectrum is averaged
%   method: string, method of averaging: Mean, median, exponentiated mean of logvalues (default)

%% outputs:
%   sn: scalar, std of the noise

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% adapted from the MATLAB implemention by Eftychios Pnevmatikakis and the
% Python implementation from Johannes Friedrich

%% References
% Pnevmatikakis E. et.al., Neuron 2016, Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data

%% input arguments
if ~exist('range_ff', 'var') || isempty(range_ff)
    range_ff = [.25, .5];
end
if ~exist('method', 'var') || isempty(method)
    method = 'logmexp';
end

block_size = [64, 64];
split_data = false;

dims = ndims(Y);
sizY = size(Y);
N = sizY(end);
Fs = 1;
ff = 0:Fs/N:Fs/2;
indf=ff>range_ff(1);
indf(ff>range_ff(2))=0;
if dims > 1
    d = prod(sizY(1:dims-1));
    Y = reshape(Y,d,N);
    Nb = prod(block_size);
    SN = cell(ceil(d/Nb),1);
    PSDX = cell(ceil(d/Nb),1);
    if ~split_data
        for ind = 1:ceil(d/Nb);
            xdft = fft(Y((ind-1)*Nb+1:min(ind*Nb,d),:),[],2);
            xdft = xdft(:,1: floor(N/2)+1); % FN: floor added.
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(:,2:end-1) = 2*psdx(:,2:end-1);
            %SN{ind} = mean_psd(psdx(:,indf),method);
            switch method
                case 'mean'
                    SN{ind}=sqrt(mean(psdx(:,indf)/2,2));
                case 'median'
                    SN{ind}=sqrt(median(psdx(:,indf)/2),2);
                case 'logmexp'
                    SN{ind} = sqrt(exp(mean(log(psdx(:,indf)/2),2)));
            end
            PSDX{ind} = psdx;
        end
    else
        nc = ceil(d/Nb);
        Yc = mat2cell(Y,[Nb*ones(nc-1,1);d-(nc-1)*Nb],N);
        parfor ind = 1:ceil(d/Nb);
            xdft = fft(Yc{ind},[],2);
            xdft = xdft(:,1:floor(N/2)+1);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(:,2:end-1) = 2*psdx(:,2:end-1);
            Yc{ind} = [];
            switch method
                case 'mean'
                    SN{ind}=sqrt(mean(psdx(:,indf)/2,2));
                case 'median'
                    SN{ind}=sqrt(median(psdx(:,indf)/2),2);
                case 'logmexp'
                    SN{ind} = sqrt(exp(mean(log(psdx(:,indf)/2),2)));
            end
            
        end
    end
    sn = cell2mat(SN);
else
    xdft = fft(Y);
    xdft = xdft(:,1:floor(N/2)+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(:,2:end-1) = 2*psdx(:,2:end-1);
    switch method
        case 'mean'
            sn = sqrt(mean(psdx(:,indf)/2,2));
        case 'median'
            sn = sqrt(median(psdx(:,indf)/2),2);
        case 'logmexp'
            sn = sqrt(exp(mean(log(psdx(:,indf)/2),2)));
    end
end
if dims > 2
    sn = reshape(sn,sizY(1:dims-1));
end
