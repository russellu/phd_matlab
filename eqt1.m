t1im = t1.img ; 
allimgs = zeros(8,240,240,150) ; 
stepcount = 1 ;
for step = 5:12
newimg = zeros(size(t1im)) ; 
stepcount
for i=1:step:size(t1im,1)-step
    for j=1:step:size(t1im,2)-step
        for k=1:step:size(t1im,3)-step
            imgi = reshape(t1im(i:i+step-1,j:j+step-1,k:k+step-1),[step,step*step]) ;
            imgi = localhist(imgi) ; 
            newimg(i:i+step-1,j:j+step-1,k:k+step-1) = reshape(imgi,[step,step,step]) ; 
        end
    end
end
allimgs(stepcount,:,:,:) = newimg ; 
stepcount = stepcount + 1 ;
end
beep



%for i=1:150 ; eqt1z(:,:,i) = localhist(squeeze(t1im(:,:,i))) ; i,end
%for i=1:240 ; eqt1x(i,:,:) = localhist(squeeze(t1im(i,:,:))) ; i,end
%for i=1:240 ; eqty(:,i,:) = localhist(squeeze(t1im(:,i,:))) ; i,end


