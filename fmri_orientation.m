clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};
subdates = {'07/06/2018','07/07/2018','01/06/2018','01/06/2018','17/05/2018','17/05/2018','11/05/2018','11/05/2018','12/04/2018','26/04/2018','26/04/2018','12/04/2018','03/05/2018','17/05/2018'};
subcomps = {[31,21],[21],[20],[23,40],[60,33],[43,26,18],[27,52],[13,50],[62,43],[5,69],[44],[81],[54],[18,26]}; 
cd e:/saved;  filtmix = load('mean_mix'); filtmix = (filtmix.mean_mix); 

for sb=1:length(subs) ; disp(subs{sb}); 
    tr = 0.68; 
    cd(['e:\orientation_retinotopy\',subs{sb}]);    
    fullfov = load_untouch_nii('bp_fullfov.nii.gz');     
    lowcontrast = load_untouch_nii('blur_tproj_topup_mc_orientation_2.nii.gz');
    highcontrast = load_untouch_nii('blur_tproj_topup_mc_orientation_1.nii.gz'); 
    quad_design = zeros(1,round(750*tr)); 
    for i=1:2:24; quad_design((i-1)*10+1:(i-1)*10+10) = 1; end
    for i=1:2:24; quad_design((i-1)*10+1+250:(i-1)*10+10+250) = 1; end
    quad_design = circshift(quad_design,6); quad_design = smooth(imresize(quad_design,[1,750])); 
    hemis_design = quad_design; 
    fullfov_design = quad_design; 
    h1_fullfov = fullfov.img(:,:,:,1:end/2); 
    h2_fullfov = fullfov.img(:,:,:,376:end); 
    h1_fullfov_design = fullfov_design(1:end/2);
    h2_fullfov_design = fullfov_design(376:end); 
    
    corrs1_fullfov = voxcorr(h1_fullfov(:,:,:,20:end-20),h1_fullfov_design(20:end-20)); 
    corrs2_fullfov = voxcorr(h2_fullfov(:,:,:,20:end-20),h2_fullfov_design(20:end-20)); 
    
    fovcorrs{sb} = corrs1_fullfov; 
    
    [sv,si] = sort(corrs1_fullfov(:),'ascend'); 

    
    [sv,si2] = sort(corrs2_fullfov(:),'ascend'); 

    res_fullfov = reshape(fullfov.img,[numel(fullfov.img(:,:,:,1)),size(fullfov.img,4)]); 
      
    res_lowcontrast = reshape(lowcontrast.img,[numel(lowcontrast.img(:,:,:,1)),size(lowcontrast.img,4)]); 
    res_highcontrast = reshape(highcontrast.img,[numel(highcontrast.img(:,:,:,1)),size(highcontrast.img,4)]); 
    
    fovlow(sb,:,:) = res_lowcontrast(si(1:500),:); 
    fovhigh(sb,:,:) = res_highcontrast(si(1:500),:); 
    fulllow(sb,:,:) = res_lowcontrast(si2(1:500),:); 
    fullhigh(sb,:,:) = res_highcontrast(si2(1:500),:); 
       
end

startind = 10/0.68 + 6; 
inds = round(startind:startind + (8*60)/0.68); 

orvalshigh = (squeeze(mean(fovhigh(:,1:50,inds),2))');
orvalslow = (squeeze(mean(fovlow(:,1:50,inds),2))');
tincr = length(orvalshigh)/8; 
icount=1;
clear eporvals
for i=1:tincr:length(orvalshigh)
   eporvalshigh(icount,:,:) = orvalshigh(floor(i):floor(i)+floor(tincr),:);  
   eporvalslow(icount,:,:) = orvalslow(floor(i):floor(i)+floor(tincr),:);  
   icount=icount+1; 
end


plot(squeeze(mean(mean(eporvalshigh,1),3))); hold on ;
plot(squeeze(mean(mean(eporvalslow,1),3))); 
%goodsubs = [1:6,8,10:14];
goodsubs = 1:14; 
clear res_orvalshigh res_orvalslow 
for i=1:8
   res_orvalshigh(i,:,:) = (imresize(squeeze(circshift(eporvalshigh(i,:,goodsubs),-22)),[360,length(goodsubs)])); 
   res_orvalslow(i,:,:) = (imresize(squeeze(circshift(eporvalslow(i,:,goodsubs),-22)),[360,length(goodsubs)])); 
end

mhigh = (res_orvalshigh(2:end,1:180,:) + res_orvalshigh(2:end,181:end,:)) /2; 
mlow = (res_orvalslow(2:end,1:180,:) + res_orvalslow(2:end,181:end,:)) /2; 

subplot(2,2,1); 
shadedErrorBar(1:180,squeeze(mean(mean(mhigh,1),3)),squeeze(std(mean(mhigh,1),0,3))/sqrt(14),'r'); hold on ; 
shadedErrorBar(1:180,squeeze(mean(mean(mlow,1),3)),squeeze(std(mean(mlow,1),0,3))/sqrt(14),'b'); xlabel('orientation (deg\circ)') ; vline(90,'k') ; vline([45,135]); 
oblique = [35:55,125:145];
vertical = 80:100; 

oblique_bars = [squeeze(mean(mean(mhigh(:,oblique,:),1),2)), squeeze(mean(mean(mlow(:,oblique,:),1),2))];
vertical_bars = [squeeze(mean(mean(mhigh(:,vertical,:),1),2)), squeeze(mean(mean(mlow(:,vertical,:),1),2))];
subplot(2,2,3); 
[h,p,ci,stats] = ttest(oblique_bars(:,1),oblique_bars(:,2)); 
barwitherr(squeeze(std(oblique_bars,0,1))/sqrt(14),squeeze(mean(oblique_bars,1))); title(['@oblique, ',format_t(stats.tstat),' ',format_p(p)]); 
subplot(2,2,4); 
barwitherr(squeeze(std(vertical_bars,0,1))/sqrt(14),squeeze(mean(vertical_bars,1)));
[h,p,ci,stats] = ttest(vertical_bars(:,1),vertical_bars(:,2));  title(['@vertical, ',format_t(stats.tstat),' ',format_p(p)]);

% todo -get the components, and compute the hrf delay, try to get it with
% no jaggies


%{
for i=1:8
   res_orvalshigh(i,:,:) = circshift(imresize(squeeze(eporvalshigh(i,:,:)),[360,14]),-90); 
   res_orvalslow(i,:,:) = circshift(imresize(squeeze(eporvalslow(i,:,:)),[360,14]),-90); 

end
%}


res_orvalshigh = res_orvalshigh(2:end,1:180,:)/2 + res_orvalshigh(2:end,181:end,:)/2;
res_orvalslow = res_orvalslow(2:end,1:180,:)/2 + res_orvalslow(2:end,181:end,:)/2;

res_orvalshigh = zscore(res_orvalshigh,[],2); 
res_orvalslow = zscore(res_orvalslow,[],2); 

subplot(1,2,1); 
plot((squeeze(mean(mean(res_orvalshigh,1),3)))) ; 
subplot(1,2,2); 
plot((squeeze(mean(mean(res_orvalslow,1),3)))) ; 

bold_lowvals = squeeze(mean(res_orvalslow,1)); 
bold_highvals = squeeze(mean(res_orvalshigh,1)); 

cd e:/saved; 
save('bold_lowvals','bold_lowvals');
save('bold_highvals','bold_highvals'); 


allhigh = mhigh; save('allhigh','allhigh');
alllow = mlow; save('alllow','alllow'); 
figure,subplot(1,2,1); 
plot(squeeze(mean(res_orvalslow,1)),'r')
subplot(1,2,2); 
plot(squeeze(mean(res_orvalshigh,1)),'b')

