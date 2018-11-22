cd c:/shared/badger_mri/russell/nii ; ls 
t1 = load_untouch_nii('t1.nii') ; 
img = double(squeeze(t1.img(:,:,75))) ; 

imagesc(log(abs(fftshift(fft2(img))))) ; 
fimg = fftshift(fft2(img)) ; 

fimg(111:130,:) = fimg(1:20,:) ; 
ifimg = ifft2(fftshift(fimg)) ; 

for j=1:10 ; 
for i=1:240 
    
    ifi = fftshift(fft2(img)) ; 
    zimg = zeros(size(ifi)) ; 
    zimg(:,1:i) = ifi(:,1:i) ; 
    ifi(:,1:i)= 0 ; 
    invi = real(ifft2(fftshift(zimg))) ; 
    
    pres = log(abs(zimg)) ; 
    
    subplot(1,2,1) ; imagesc(pres) ; colormap gray ; title('frequency space') ; 
    subplot(1,2,2) ; imagesc(invi) ; colormap gray ; title('normal space') ;
    getframe ; 
end
end


% corrupted acquisition

for j=1:10 ; 
    
    pzimg = fftshift(fft2(img)) ; 
    
for i=1:239 
    
    ifi = fftshift(fft2(img)) ; 
    %zimg = zeros(size(ifi)) ; 
    pzimg(:,i:i+1) = ifi(:,i:i+1) ; 
    randinds = ceil(rand(1,240)*239) ; 
    pzimg(randinds(1:100),i:i+1)= pzimg(randinds(101:200),i:i+1) ; 
    invi = real(ifft2(fftshift(pzimg))) ; 
    
    pres = log(abs(pzimg)) ; 
    
    subplot(1,2,1) ; imagesc(pres) ; colormap gray ; title('frequency space') ; 
    subplot(1,2,2) ; imagesc(invi) ; colormap gray ; title('normal space') ;
    getframe ; 
end
end

