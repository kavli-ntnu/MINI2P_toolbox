%% Perceptually Uniform Colormaps from MatPlotLib
% The <www.mathworks.com/matlabcentral/fileexchange/62729
% MatPlotLib Perceptually Uniform Colormaps> submission includes the
% default colormap family and default line colororder family from
% MatPlotLib 2 and 3. This document shows examples of their usage.
%% Overview
matplotlib_plot
%% |VIRIDIS| Default Colormap
close()
load spine
image(X)
colormap(viridis)
%% |CIVIDIS|
colormap(cividis)
%% |INFERNO|
colormap(inferno)
%% |MAGMA|
colormap(magma)
%% |PLASMA|
colormap(plasma)
%% |TWILIGHT|  (cyclical)
colormap(twilight)
%% |TAB10| Default Line ColorOrder
N = 20;
clf()
axes('ColorOrder',tab10(N),'NextPlot','replacechildren')
X = linspace(0,pi*3,1000);
Y = bsxfun(@(x,n)n*sin(x+2*n*pi/N), X(:), 1:N);
plot(X,Y, 'linewidth',4)
%% |TAB20|
clf()
axes('ColorOrder',tab20(N),'NextPlot','replacechildren')
plot(X,Y, 'linewidth',4)
%% |TAB20B|
clf()
axes('ColorOrder',tab20b(N),'NextPlot','replacechildren')
plot(X,Y, 'linewidth',4)
%% |TAB20C|
clf()
axes('ColorOrder',tab20c(N),'NextPlot','replacechildren')
plot(X,Y, 'linewidth',4)