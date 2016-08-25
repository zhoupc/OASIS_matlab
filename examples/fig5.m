%%
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]};
figure; 
init_fig; 

%% 
gam = 0.95; 
sn = .3; 
T = 1500; 
[Y, trueC, trueSpikes] = gen_sinusoidal_data(gam, sn, T); 