function [fftimg,gsimg]=AnalyzeWithFFT(img,fftsz,channel)
if(nargin<2) channel=1; end; % default red channel
if(channel)
    tmpgsimg=img(:,:,channel);
else
    tmpgsimg=rgb2gray(img);
end;
    
% adapt histogram
mask = adapthisteq(tmpgsimg);

% erode
se = strel('disk',5);
marker = imerode(mask,se);

% reconstruct better image
gsimg = imreconstruct(marker,mask);

% get fft image
fftimg=plotfft(gsimg',-fftsz);

end