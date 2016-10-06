# OASIS: Fast online deconvolution of calcium imaging data

Code accompanying the paper "Fast Active Set Method for Online Spike Inference from Calcium Imaging". [arXiv, 2016]

Python: https://github.com/j-friedrich/OASIS

MATLAB: https://github.com/zhoupc/OASIS_matlab

### Installation
Run setup.m to add OASIS function to the search path of MATLAB

`>> setup`

### Examples
The scripts to reproduce the figures are in the subfolder 'examples' with names obvious from the arXiv paper. They can be run with 

`>> run examples/paper/fig*.m ` 

There is also a function **deconvolveCa.m** wrapping all formulations of the problem. You might only need this one for your problem. See a list of examples  in the demo script **oasis_test.m**, 

`>> edit examples/oasis_test.m`

### Benchmarks