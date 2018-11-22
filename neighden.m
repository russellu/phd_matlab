cd c:/shared/nbhden ;
ls
t1 = load_nii('T1.nii.gz') ; 
t1im = t1.img ; 
chunk = t1im(50:120,100:150,50:80) ; 

nsize = 5 ; 
floorsize = floor(nsize/2) ; 
[x,y,z] = ndgrid(-floorsize:floorsize,-floorsize:floorsize,-floorsize:floorsize) ; 
xs = zeros(nsize,nsize,nsize) ; 
ys = zeros(nsize,nsize,nsize) ; 
zs = zeros(nsize,nsize,nsize) ;  
xs(x==0) = 1 ;
ys(y==0) = 1 ; 
zs(z==0) = 1 ; 

stdxchunk = zeros(size(chunk)) ; 
stdychunk = zeros(size(chunk)) ; 
stdzchunk = zeros(size(chunk)) ; 

for i=floorsize+1:size(chunk,1)-floorsize
    for j=floorsize+1:size(chunk,2)-floorsize
        for k=floorsize+1:size(chunk,3)-floorsize
            c = chunk(i-floorsize:i+floorsize,j-floorsize:j+floorsize,k-floorsize:k+floorsize) ;
            stdxchunk(i,j,k) = std(nonzeros(c.*xs)) ;  
            stdychunk(i,j,k) = std(nonzeros(c.*ys)) ;  
            stdzchunk(i,j,k) = std(nonzeros(c.*zs)) ;  
        end
    end
end
stdxchunk = norm255(stdxchunk) ; stdxflat = reshape(stdxchunk,[size(stdxchunk,1),size(stdxchunk,2)*size(stdxchunk,3)]) ; 
stdychunk = norm255(stdychunk) ; stdyflat = reshape(stdychunk,[size(stdychunk,1),size(stdychunk,2)*size(stdychunk,3)]) ; 
stdzchunk = norm255(stdzchunk) ; stdzflat = reshape(stdzchunk,[size(stdzchunk,1),size(stdzchunk,2)*size(stdzchunk,3)]) ; 
[~,xk] = kmeans(stdxflat,6) ; 
[~,yk] = kmeans(stdyflat,6) ; 
[~,zk] = kmeans(stdzflat,6) ; 





for i=1:31
    subplot(6,6,i) ; 
    imagesc(squeeze(stdxchunk(:,:,i))) ;
end
figure,
for i=1:31
    subplot(6,6,i) ; 
    imagesc(squeeze(chunk(:,:,i)),[0,0.5]) ; 
end