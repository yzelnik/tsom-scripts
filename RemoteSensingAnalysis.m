%% Scripts for remote sensing analysis, comparing it to the vegetation rehabilitation model

% choose which image to show and analyze
yearchs = 2019;
areanum = 117;

% load image data
rgbimg = imread(sprintf('area%d_rgb%d.png',areanum,yearchs));
gsimg  = imread(sprintf('area%d_gs%d.tif',areanum,yearchs));

% plot it out
subplot(1,2,1);
imagesc(rgbimg);
subplot(1,2,2);
imagesc(gsimg);
colormap gray

%% Clean up rgb image and plot the FFT analysis
fftsz=7; % number of frequencies in FFT

[fftimg,gsimg2]=AnalyzeWithFFT(rgbimg,fftsz,1);
% Note that at this point, gsimg and gsimg2 are virtually identical.
% This is just to show how we got our gsimg (gray-scale) representation

% Plot out FFT results next to the gray-scale one
subplot(1,2,1);
imagesc(fftimg);
subplot(1,2,2);
imagesc(gsimg2);
