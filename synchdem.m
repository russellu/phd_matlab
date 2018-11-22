
clear img1 img2 img3
for i=1:100 ;
    for j=1:100
        randij1 = 0 ; 
        randij2 = rand*2 ; 
        randij3 = rand*20 ; 
        img1(i,j,:) = sin(randij1:.1:100+randij1) ;
        img2(i,j,:) = sin(randij2:.1:100+randij2) ;
        img3(i,j,:) = sin(randij3:.1:100+randij3) ;
    end
end
disp3d(img2)
a = squeeze(sum(sum(img1,1),2)) ;
sig1 = squeeze(sum(sum(img1,1),2)) ; 
sig2 = squeeze(sum(sum(img2,1),2)) ; 
sig3 = squeeze(sum(sum(img3,1),2)) ; 
sigmin = min(sig1) ; sigmax = max(sig1) ;

for i=100:size(img1,3)
    subplot(2,3,1) ; imagesc(squeeze(img1(:,:,i)),[-1,1]) ; 
    subplot(2,3,2) ; imagesc(squeeze(img2(:,:,i)),[-1,1]) ; 
    subplot(2,3,3) ; imagesc(squeeze(img3(:,:,i)),[-1,1]) ; 
    subplot(2,3,4) ; plot(sig1(i-99:i)) ; ylim([sigmin,sigmax]) ; 
    subplot(2,3,5) ; plot(sig2(i-99:i)) ; ylim([sigmin,sigmax]) ; 
    subplot(2,3,6) ; plot(sig3(i-99:i)) ; ylim([sigmin,sigmax]) ; 
    getframe ; 
end