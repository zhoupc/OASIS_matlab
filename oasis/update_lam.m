function [c, Aset, lam, s] = update_lam(y, c, Aset, g, lam, thresh)
%% update the tuning parameter lambda to reduce |s|_1 while |y-c|_2^2 <= sn^2*T

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities 
        %withone entry per time-bin.
%   c: T*1 vector, previous c 
%   Aset: npools*4 matrix, previous active sets 
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
        %impulse response.   
%   lam:  scalar, curret value of sparsity penalty parameter lambda. 
%   thresh: maximum residual sn*T^2
       
%% outputs
%   c: T*1 vector
%   s: T*1 vector, spike train 
%   Aset: npool x 4 matrix, active sets 
%   lam: scalar, new tuning parameter 

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References 
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

c = reshape(c, [],1); 
y = reshape(y, [],1); 
len_active_set = size(Aset,1); 

res = y - c; 
RSS = res'*res; 

temp = zeros(size(c)); 
for ii=1:len_active_set
    ti = Aset(ii, 3); 
    li = Aset(ii, 4); 
    idx = 0:(li-1); 
    temp(ti+idx) = (1-g^li)/ Aset(ii,2) * g.^(idx);   
end

aa = temp'*temp; 
bb = res'*temp; 
cc = RSS-thresh; 
ll = (-bb + sqrt(bb^2-aa*cc)) / aa; 
lam = lam + ll; 

Aset(:,1) = Aset(:,1) - ll*(1-g.^Aset(:,4)); 
[c, s, Aset] = oasisAR1(y, g, lam, 0, Aset); 

























