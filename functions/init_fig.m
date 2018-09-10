if ~exist('fig_x0', 'var')
    fig_x0 = 0;
end
if ~exist('fig_y0', 'var')
    fig_y0 = 0;
end
set(gcf, 'color', 'w', ...
    'defaultAxesFontSize', 20, ...)
    'defaultlinelinewidth',2, ...
    'position', [fig_x0, fig_y0, 100*get(gcf, 'papersize')], ...
    'paperposition', [0, 0, get(gcf, 'papersize')]); 
clear fig_x0 fig_y0; 
