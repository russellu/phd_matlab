cd c:/shared/denoise2
ls
t1 = load_nii('fslice.nii.gz')  ;
t1im = t1.img ; 
stdt1 = load_nii('fslice_stdev.nii.gz') ; 
stdim = stdt1.img ; 

%t1chunk = (squeeze(t1im(80:160,40:120,50:89))) ; 
t1chunk = (squeeze(t1im(:,:,:))) ; 
%stdchunk = (squeeze(stdim(80:160,40:120,50:89))) ; 
stdchunk = (squeeze(stdim(:,:,:))) ; 
stdchunk = norm255(stdchunk) ;

flatstd = reshape(stdchunk,[size(stdchunk,1),size(stdchunk,2)*size(stdchunk,3)]) ; 
[~,clusts] = kmeans(flatstd,6) ; 
stdclusts = reshape(clusts,size(stdchunk)) ; 

nbhds = [5,3,1,1,1,1] ;  
%maxnbhd = floor(nbhds(size(nbhds,2))/2) ; 
maxnbhd = floor(nbhds(1)/2) ; 

cleaned = zeros(size(t1chunk)) ;

for i=maxnbhd+1:size(t1chunk,1)-maxnbhd
    for j=maxnbhd+1:size(t1chunk,2)-maxnbhd
        for k=maxnbhd+1:size(t1chunk,3)-maxnbhd
            vali = stdclusts(i,j,k) ; 
            csize = floor(nbhds(vali)/2) ; 
            vol = t1chunk(i-csize:i+csize,j-csize:j+csize,k-csize:k+csize) ; 
            linevol = reshape(vol,[1,size(vol,1)*size(vol,2)*size(vol,3)]) ;
            med = median(double(linevol)) ; 
            cleaned(i,j,k) = med ;
        end
    end
end



    

%%% find some way to maximize the "noise removal" measure and then iterate






