clear all ; close all; 
subs = {'badger_alex','badger_dina','badger_genevieve','badger_jeremie','badger_russell','badger_sukhman','badger_tegan','badger_valerie'};

for sb=1%:length(subs)
    cd(['e:/fmris/',subs{sb},'/mel']);
    ls; 
    
    mix = load('melodic_mix'); 
    [pxx,f] = pwelch(double(mix),40,30,20,1/0.693); 
    pxx = log(pxx); 
    atlas = load_untouch_nii('../epi_aseg.nii.gz'); 
    atlasinds = unique(atlas.img); atlasinds = atlasinds(2:end); 
    
    rois = load_untouch_nii('melodic_IC.nii.gz'); 
    res_rois = mat2gray(reshape(rois.img,[numel(rois.img(:,:,:,1)),size(rois.img,4)])); 
    clear atlasvals; 
    for i=1:length(atlasinds)
       roi_inds = find(atlas.img == atlasinds(i));  
       for j=1:size(rois.img,4)
          atlasvals(i,j) = log(mean(res_rois(roi_inds,j)));           
       end
    end
    clear features;
    %features(:,1:size(atlasvals,1)) = mat2gray(atlasvals'); 
    %features(:,size(atlasvals,1)+1:size(pxx,1)+size(atlasvals,1)) = mat2gray(pxx'); 
    %features(:,size(pxx,1)+1:size(pxx,1)*2-1) = mat2gray(diff(pxx,1,1))'; 
    %features(:,size(pxx,1)+1:size(pxx,1)+size(atlasvals,1)) = mat2gray(atlasvals'); 
    features(:,1:size(pxx,1)) = (pxx(1:end,:)'); 
    %features(:,size(pxx,1)+1:size(pxx,1)+size(atlasvals,1)) = mat2gray(atlasvals');
    
    k = kmeans(features,4); 
    
    %for i=1:10
    %    subplot(4,5,i*2-1) ; plot(mix(end-449:end,(k==i))); 
    %end
    
    
    for i=1:4
       inds = find(k==i);
       figure;
       for j=1:length(inds)
           subplottight(ceil(sqrt(length(inds))),ceil(sqrt(length(inds))),j); 
           imagesc(imrotate(squeeze(max(rois.img(:,:,:,inds(j)),[],1)),90)); colormap jet; 
           
       end

    end
    
    
end


badcs = [3,1,4]; 
badinds = zeros(1,length(k)); 
for i=1:length(badcs)
    inds = find(k==badcs(i));
    badinds(inds) = 1; 
end

fnames = {'retino_allstims_01','retino_allstims_02','retino_gamma_01','retino_gamma_02','retino_movie','retino_rest'}; 
goodcs = mix(:,(badinds==0)); 
badcs = mix(:,(badinds==1)); clear inds; 
inds{1} = 1:735; inds{2} = 736:1470; inds{3} = 1471:2205; inds{4} = 2206:2940; inds{5} = 2941:3675; inds{6} = 3676:3676+449; 

for i=1:size(badcs,2) ; normbadcs(:,i) = mat2gray(badcs(:,i)); end;

kbad = kmeans(normbadcs',15); 
for i=1:15
   subplot(3,5,i); 
   plot(mean(badcs(1:735,find(kbad==i)),2)); 
   nbadcs(:,i) = mean(badcs(:,find(kbad==i)),2); 
end



for i=1:length(inds)
dlmwrite(['good_',fnames{i},'.txt'],goodcs(inds{i},:)); 
dlmwrite(['bad_',fnames{i},'.txt'],nbadcs(inds{i},:)); 
end

badinds = find(badinds==1); 
dlmwrite('badinds.txt',badinds); 




