clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};
subdates = {'07/06/2018','07/07/2018','01/06/2018','01/06/2018','17/05/2018','17/05/2018','11/05/2018','11/05/2018','12/04/2018','26/04/2018','26/04/2018','12/04/2018','03/05/2018','17/05/2018'};
subcomps = {[31,21],[21],[20],[23,40],[60,33],[43,26,18],[27,52],[13,50],[62,43],[5,69],[44],[81],[54],[18,26]}; 
cd e:/saved;  filtmix = load('mean_mix'); filtmix = (filtmix.mean_mix); 
    tr = 0.68; 

quad_design = zeros(1,round(750*tr)); 
for i=1:2:24; quad_design((i-1)*10+1:(i-1)*10+10) = 1; end
for i=1:2:24; quad_design((i-1)*10+1+250:(i-1)*10+10+250) = 1; end
quad_design = circshift(quad_design,6); quad_design = smooth(imresize(quad_design,[1,750])); 
hemis_design = quad_design; 
fullfov_design = quad_design; 
h1_fullfov_design = fullfov_design(1:end/2);
h2_fullfov_design = fullfov_design(376:end); 

for sb=1:length(subs) 

    cd(['e:\orientation_retinotopy\',subs{sb}]);
    lowcontrast = load_untouch_nii('bp_mc_orientation_2.nii');
    highcontrast = load_untouch_nii('bp_mc_orientation_1.nii'); 
    
    cd(['e:\orientation_retinotopy\',subs{sb},'\mel']);
    mix = load('melodic_mix'); 
    mix = eegfiltfft(mix',1/0.68,0.01,2)'; 
    clear corrs; 
    for i=1:size(mix,2)
       corrs(i) = corr(mix(20:375-20,i),h1_fullfov_design(20:375-20));  
        
    end
    
    [sv,si] = sort(corrs,'descend'); 
    
    res_lowcontrast = reshape(lowcontrast.img,[numel(lowcontrast.img(:,:,:,1)),size(lowcontrast.img,4)]); 
    res_highcontrast = reshape(highcontrast.img,[numel(highcontrast.img(:,:,:,1)),size(highcontrast.img,4)]); 
    
    fovlow(sb,:,:) = mix(751:751+749,si(1)); 
    fovhigh(sb,:,:) = mix(1501:1501+749,si(1)); 

end


startind = 10/0.68 + 6; 
inds = round(startind:startind + (8*60)/0.68); 

orvalshigh = (squeeze((fovhigh(:,inds)))');
orvalslow = (squeeze((fovlow(:,inds)))');
tincr = length(orvalshigh)/8; 
icount=1;
clear eporvalshigh eporvalslow
for i=1:tincr:length(orvalshigh)
   eporvalshigh(icount,:,:) = orvalshigh(floor(i):floor(i)+floor(tincr),:);  
   eporvalslow(icount,:,:) = orvalslow(floor(i):floor(i)+floor(tincr),:);  
   icount=icount+1; 
end


plot(squeeze(mean(mean(eporvalshigh,1),3))); hold on ;
plot(squeeze(mean(mean(eporvalslow,1),3))); 




