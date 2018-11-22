clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','valerie'}; 
scans = {'bp_reg_topup_mc_retino_allstims_01','bp_reg_topup_mc_retino_allstims_02','bp_reg_topup_mc_retino_gamma_01','bp_reg_topup_mc_retino_gamma_02','bp_reg_topup_mc_retino_movie','bp_reg_topup_mc_retino_rest'}; 

for sb=1:length(subs)
    cd(['E:\rawbadger\badger_mri\',subs{sb},'\nii']);   
    fcorrs1 = dir('warp_fcorr1*gz');
    fcorrs2 = dir('warp_fcorr2*gz'); 
    
    for fc=1:length(fcorrs1)
        fc1_nii = load_untouch_nii(fcorrs1(fc).name); 
        allcorrs1(sb,fc,:,:,:,:) = fc1_nii.img; 
        fc2_nii = load_untouch_nii(fcorrs2(fc).name); 
        allcorrs2(sb,fc,:,:,:,:) = fc2_nii.img; 
    end
    disp(sb); 
    
    
    
    
end


cd e:/rawbadger/badger_mri/russell/nii ; ls
f41 = load_untouch_nii('f_41.nii.gz'); 
cd e:/meanepis;

mcorrs1 = squeeze(mean(mean(allcorrs1,1),2)); 
mcorrs2 = squeeze(mean(mean(allcorrs2,1),2)); 

f41.img = mcorrs1; 
save_untouch_nii(f41,'mcorrs1.nii.gz'); 
f41.img = mcorrs2; 
save_untouch_nii(f41,'mcorrs2.nii.gz'); 

mean_nii = load_untouch_nii('c:/shared/epireg/mean.nii.gz'); 
tps=18:27; 
times=(-20:20)*0.693; 
for tp=1:length(tps)
m1 = squeeze(mean(mean(allcorrs1(:,:,:,:,:,tps(tp)),1),2));
m2 = squeeze(mean(mean(allcorrs2(:,:,:,:,:,tps(tp)),1),2));

m1(isnan(m1)) = 0 ; m2(isnan(m2)) = 0; 



subplot(2,10,tp);
plotoverlayIntensity2D(squeeze(mean_nii.img(32,:,:)),squeeze(mat2gray(abs(mean(m1(30:34,:,:),1)))),squeeze(mean(m1(30:34,:,:),1)),90); title(times(tps(tp))); 
subplot(2,10,10+tp); 
plotoverlayIntensity2D(squeeze(mean_nii.img(32,:,:)),squeeze(mat2gray(abs(mean(m2(30:34,:,:),1)))),squeeze(mean(m2(30:34,:,:),1)),90); 

end




