function fftimg=plotfft(img,showsz)
if nargin<2 showsz=0; end;
if showsz(1)<0 showsz=-showsz; plotflag=0; else plotflag=1; end;

img=img(:,:,1); 
if(size(img,1) || size(img,2)==1) % reshape to square
    newsz = sqrt(size(img,1)*size(img,2));
    img=reshape(img,newsz,newsz);
end; 
if showsz(1)==0 
    showsz = round(0.2*size(img));
elseif showsz(1)<1
    showsz = round(showsz(:)'.*size(img));
end;
if(length(showsz)<2) 
    showsz=[showsz showsz]; 
end;

%showsz = round(showsz/2)*2;
zoom = [round(size(img)/2)-showsz+1 round(size(img)/2)+showsz+1];

tmpfft=abs(fft2(img));
tmpfft(1)=0;

tmpfft=fftshift(tmpfft);

fftimg=tmpfft(zoom(1):zoom(3),zoom(2):zoom(4))';

if(plotflag)
    imagesc(fftimg);
    set(gca,'xTick',[],'yTick',[]);
end;

end