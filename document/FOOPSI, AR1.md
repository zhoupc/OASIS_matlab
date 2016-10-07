### FOOPSI, AR1

#### 1: know AR coefficient g, lambda 

`[c_oasis, s_oasis] = deconvolveCa(y, 'ar1', g, 'foopsi', 'lambda', lambda);`


#### 2: know lambda, estimate g from autocorrelation function

`[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', 'foopsi', 'lambda', lambda);` 



#### 3: know lambda, estimate g from autocorrelation function, then optimize it 

`[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', 'foopsi', 'lambda', lambda, 'optimize_pars');` 


