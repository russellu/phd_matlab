cd c:/shared/
a = load_nii('a.nii.gz' ) ;
%a = load_untouch_nii('a.nii.gz') ; 
aim = a.img ; 
x = reshape(aim,[size(aim,1),size(aim,2)*size(aim,3)]) ;
x = norm255(x) ; 





%{
nbhd=20 ; 
newim = zeros(size(aim)) ; 
for i=1:nbhd:size(aim,1)-nbhd
    disp(i)
    for j=1:nbhd:size(aim,2)-nbhd
        for k=1:nbhd:size(aim,3)-nbhd
            box = aim(i:i+nbhd,j:j+nbhd,k:k+nbhd) ;
            if sum(sum(sum(box))) > 2000
                flatbox = reshape(box,[size(box,1),size(box,2)*size(box,3)]) ; 
                [~,kbox] = kmeans(flatbox,5) ; 
                newbox = reshape(kbox,size(box)) ; 
                newim(i:i+nbhd,j:j+nbhd,k:k+nbhd) = newbox ; 
            end
        end
    end
end
%}





fim = imfilter(aim,fspecial3('log',[7,7,7])) ; 

aimr = reshape(aim,[size(aim,1),size(aim,2)*size(aim,3)]) ; 
[~,klusts] = kmeans(norm255(aimr),5) ; 
klusts(klusts==5) = 0 ; 
klusts(klusts==1) = 0 ; 
klusts(klusts~=0) = 1 ; 
klusts = double(klusts).*double(aimr) ; 
[~,k2] = kmeans(klusts,5) ; 


%%% histeq only in regions with high STDV??? avoid applying histeq in the
%%% white matter...
%%% kmeans the gradient image?



ak = reshape(klusts,size(aim)) ; 























%{
c = (squeeze(aim(60:100,100:140,90:110))) ;
for i=1:21 ; 
    subplot(3,7,i) ; 
    img = squeeze(c(:,:,i)) ; 
    [~,y] = kmeans(squeeze(c(:,:,i)),4) ; 
    y(y==4) = 0 ;
    y(y==1) = 0 ;
    y(y~=0) = 1 ; 
    imagesc(double(y).*double(img))
    x(:,:,i) = double(y).*double(img) ; 
end
%}