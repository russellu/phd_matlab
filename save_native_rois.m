clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};
for sb=2:length(subs) 
    cd(['e:\orientation_retinotopy\',subs{sb}])
    disp(subs{sb}); 
    % quads: br, bl, tl, tr, rep 10s*12, break, rep 10s*12
    tr = 0.68; 
    quad_design = zeros(1,round(750*tr)); 
    for i=1:2:24; quad_design((i-1)*10+1:(i-1)*10+10) = 1; end
    for i=1:2:24; quad_design((i-1)*10+1+250:(i-1)*10+10+250) = 1; end
    quad_design = circshift(quad_design,6); quad_design = smooth(imresize(quad_design,[1,750])); 
    hemis_design = quad_design; 
    fullfov_design = quad_design; 
    % hemis
    % contrast
    % fullfov
    % orientation

    clear allcorrs
    
    fullfov = load_untouch_nii('reg_bp_mc_fullfov.nii'); 
    h1_fullfov = fullfov.img(:,:,:,1:end/2); 
    h2_fullfov = fullfov.img(:,:,:,376:end);
    h1_fullfov_design = fullfov_design(1:end/2);
    h2_fullfov_design = fullfov_design(376:end); 
    corrs1_fullfov = voxcorr(h1_fullfov(:,:,:,20:end-20),h1_fullfov_design(20:end-20)); 
    corrs2_fullfov = voxcorr(h2_fullfov(:,:,:,20:end-20),h2_fullfov_design(20:end-20)); 

    allcorrs(:,:,:,5)  =corrs1_fullfov;
    allcorrs(:,:,:,6) = corrs2_fullfov; 

    quads = load_untouch_nii('reg_bp_mc_quads.nii'); 
    h1_quads = quads.img(:,:,:,1:end/2); 
    h2_quads = quads.img(:,:,:,376:end);
    h1_quad_design = quad_design(1:end/2);
    h2_quad_design = quad_design(376:end); 
    corrs1_quads = voxcorr(h1_quads(:,:,:,20:end-20),h1_quad_design(20:end-20)); 
    corrs2_quads = voxcorr(h2_quads(:,:,:,20:end-20),h2_quad_design(20:end-20)); 

    allcorrs(:,:,:,1) = corrs1_quads;
    allcorrs(:,:,:,2) = corrs2_quads; 
    
    hemis = load_untouch_nii('reg_bp_mc_hemis.nii'); 
    h1_hemis = hemis.img(:,:,:,1:end/2); 
    h2_hemis = hemis.img(:,:,:,376:end);
    h1_hemis_design = hemis_design(1:end/2);
    h2_hemis_design = hemis_design(376:end); 
    corrs1_hemis = voxcorr(h1_hemis(:,:,:,20:end-20),h1_hemis_design(20:end-20)); 
    corrs2_hemis = voxcorr(h2_hemis(:,:,:,20:end-20),h2_hemis_design(20:end-20)); 
    
    allcorrs(:,:,:,3) = corrs1_hemis;
    allcorrs(:,:,:,4) = corrs2_hemis; 
    
    ors1 = load_untouch_nii('reg_bp_mc_orientation_1.nii'); 
    ors2 = load_untouch_nii('reg_bp_mc_orientation_2.nii'); 
    ors1.img = (ors2.img + ors1.img)/2; 
    res_ors1 = reshape(ors1.img,[numel(ors1.img(:,:,:,1)),750]); 
    rotangles = mod(1:360*8,360) ;  rot_time = 480/0.68; rotangles = imresize(rotangles,[1,rot_time],'nearest'); 
    % convert rotangles to TR indices and shift by 6s
    %ts = zeros(1,round(750*0.68)); ts(11:11+479) = rotangles; ts = circshift(ts,6); ts = imresize(ts,[1,750],'nearest'); %
    ts = zeros(1,750); startt = round(10/0.68); ts(startt:startt+length(rotangles)-1) = rotangles; ts = circshift(ts,round(6/0.68)); 
    %f1 = load_untouch_nii('f1.nii.gz'); 

    maskinds = find(corrs2_fullfov(:) < -.5); 
    %maskinds = find(f1.img > 150); 

    ors1_ts = res_ors1(maskinds,:); 

    inds = find(diff(ts)<-200); 
    clear ors_i; 
    for i=1:length(inds)
        ors_i(i,:,:) = ors1_ts(:,inds(i)-87:inds(i)); 
    end
    or_axis = imresize(1:360,[1,88]); 
    plot(or_axis,squeeze(mean(mean(ors_i(2:end,:,:),1),2))) ; vline([45,135,225,315],'g');
    vline([1,180,360],'r'); vline([90,270],'k');
    
    all_ors(sb,:,:) = squeeze(mean(ors_i,2)); 
   
    %}
    cont1 = load_untouch_nii('reg_bp_mc_contrast.nii'); 
    res_cont1 = reshape(cont1.img,[numel(cont1.img(:,:,:,1)),size(cont1.img,4)]); 
    cont_ts = zeros(1,round(size(cont1.img,4)*0.68)); 
    for i=10:20:length(cont_ts)-10
       cont_ts(i:i+9) = 1; 
    end
    conved_ts= conv(cont_ts,spm_hrf(1),'full');
    conved_ts = conved_ts(1:length(cont_ts)); 
    conved_ts = imresize(conved_ts,[1,size(cont1.img,4)]); 
    corrs = voxcorr(cont1.img(:,:,:,10:end-10),conved_ts(10:end-10));
    corrs(isnan(corrs)) = 0; 
    [sv,si] = sort(corrs(:),'descend'); 

    incr = round(10/0.693);
    ctrigs = [14, 147+14,147*2+14,147*3+14,147*4+14];
    c1 = ctrigs(1):ctrigs(2); c2 = ctrigs(2):ctrigs(3); c3 = ctrigs(3):ctrigs(4); c4 = ctrigs(4):ctrigs(5); c5 = ctrigs(5):ctrigs(5)+147; 
    allcs = [c1;c2;c3;c4;c5];
    clear c_corrs; 
    for i=1:size(allcs,1)
       c_corrs(:,:,:,i) = voxcorr(cont1.img(:,:,:,allcs(i,10:end-10)),conved_ts(allcs(i,10:end-10)));   
    end
    
    subcorrs{sb} = allcorrs; 
   % allcorrs(:,:,:,7:11,sb) = c_corrs; 
   % allcont_ts(:,sb) = squeeze(mean(res_cont1(si(1:200),:),1));  
   
    thresh = 0.25; 
    %clear bot_right bot_left top_right top_left left_hemi right_hemi top_hemi bot_hemi full fov periph
    bot_right = double(allcorrs(:,:,:,1)>thresh);
    bot_left = double(allcorrs(:,:,:,1)<-thresh);
    top_right = double(allcorrs(:,:,:,2)<-thresh);
    top_left = double(allcorrs(:,:,:,2)>thresh);
    left_hemi = double(allcorrs(:,:,:,3)>thresh); 
    right_hemi = double(allcorrs(:,:,:,3)<-thresh); 
    top_hemi = double(allcorrs(:,:,:,4)>thresh); 
    bot_hemi = double(allcorrs(:,:,:,4)<-thresh); 
    full = double(allcorrs(:,:,:,5)<-thresh);
    fov = double(allcorrs(:,:,:,6)<-thresh);
    periph = double(allcorrs(:,:,:,6)>thresh);

    rets(:,:,:,1) = bot_right; 
    rets(:,:,:,2) = bot_left; 
    rets(:,:,:,3) = top_right; 
    rets(:,:,:,4) = top_left; 
    rets(:,:,:,5) = left_hemi; 
    rets(:,:,:,6) = right_hemi; 
    rets(:,:,:,7) = top_hemi; 
    rets(:,:,:,8) = bot_hemi; 
    rets(:,:,:,9) = full; 
    rets(:,:,:,10) = fov; 
    rets(:,:,:,11) = periph; 

    fret = load_untouch_nii('fret.nii.gz'); 
    fret.img = rets ; save_untouch_nii(fret,'native_rets.nii.gz');
    
end

% try to fit nvox * dist to reproduce the EEG scalp amplitude
% find which subjects you have both data for
barwitherr(squeeze(std(mdists,0,2))/sqrt(14),mean(mdists,2)); ylim([13,20]);


for i=1:11
    for j=1:11
     [h,p,ci,stats] = ttest(mdists(i,:),mdists(j,:)); 
     tstats(i,j) = stats.tstat; 
     pvals(i,j) = p; 
    end
end





