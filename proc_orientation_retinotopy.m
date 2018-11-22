clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};
for sb=1:length(subs) 
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

    
    fullfov = load_untouch_nii('reg_bp_mc_fullfov.nii'); 
    h1_fullfov = fullfov.img(:,:,:,1:end/2); 
    h2_fullfov = fullfov.img(:,:,:,376:end);
    h1_fullfov_design = fullfov_design(1:end/2);
    h2_fullfov_design = fullfov_design(376:end); 
    corrs1_fullfov = voxcorr(h1_fullfov(:,:,:,20:end-20),h1_fullfov_design(20:end-20)); 
    corrs2_fullfov = voxcorr(h2_fullfov(:,:,:,20:end-20),h2_fullfov_design(20:end-20)); 
    
    allcorrs(:,:,:,5,sb)  =corrs1_fullfov;
    allcorrs(:,:,:,6,sb) = corrs2_fullfov; 
    
    
    quads = load_untouch_nii('reg_bp_mc_quads.nii'); 
    h1_quads = quads.img(:,:,:,1:end/2); 
    h2_quads = quads.img(:,:,:,376:end);
    h1_quad_design = quad_design(1:end/2);
    h2_quad_design = quad_design(376:end); 
    corrs1_quads = voxcorr(h1_quads(:,:,:,20:end-20),h1_quad_design(20:end-20)); 
    corrs2_quads = voxcorr(h2_quads(:,:,:,20:end-20),h2_quad_design(20:end-20)); 

    allcorrs(:,:,:,1,sb) = corrs1_quads;
    allcorrs(:,:,:,2,sb) = corrs2_quads; 
    
    hemis = load_untouch_nii('reg_bp_mc_hemis.nii'); 
    h1_hemis = hemis.img(:,:,:,1:end/2); 
    h2_hemis = hemis.img(:,:,:,376:end);
    h1_hemis_design = hemis_design(1:end/2);
    h2_hemis_design = hemis_design(376:end); 
    corrs1_hemis = voxcorr(h1_hemis(:,:,:,20:end-20),h1_hemis_design(20:end-20)); 
    corrs2_hemis = voxcorr(h2_hemis(:,:,:,20:end-20),h2_hemis_design(20:end-20)); 
    
    allcorrs(:,:,:,3,sb) = corrs1_hemis;
    allcorrs(:,:,:,4,sb) = corrs2_hemis; 
    
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
    
    allcorrs(:,:,:,7:11,sb) = c_corrs; 
    allcont_ts(:,sb) = squeeze(mean(res_cont1(si(1:200),:),1)); 
    
end

figure,
mconts = squeeze(allcorrs(:,:,:,7:11,:));
for i=1:5 
    rescont = reshape(mconts(:,:,:,i,:),[80*80*42,14]); 
    
    for j=1:length(subs)
        voxconts(i,j) = sum(rescont(:,j)>.3) ; 
    end
    
end
subplot(2,2,1);
barwitherr(std(voxconts,0,2)/sqrt(7),mean(voxconts,2)); ylabel('#voxels'); 

highinds = find(conved_ts(c1)>0.6); 
lowinds = find(conved_ts(c1)<0.1); 
clear ampchange
for i=1:size(allcs,1)    
    ts_i = allcont_ts(allcs(i,:),:); 
    ampchange(i,:) = mean(ts_i(highinds,:),1) - mean(ts_i(lowinds,:),1);    
end
subplot(2,2,2); 
barwitherr(std(ampchange,0,2)/sqrt(7),mean(ampchange,2)); ylabel('BOLD amplitude'); 

subplot(2,2,3);
shadedErrorBar(or_axis,squeeze(mean(mean(all_ors,1),2)),std(mean(all_ors,2),0,1)/sqrt(14));



cd e:/orientation_retinotopy/fatlas ; ls 
mean_nii = load_untouch_nii('mean.nii.gz'); 
titles = {'ret__botquads.nii.gz','ret_topquads.nii.gz','ret_leftright.nii.gz','ret_topbot.nii.gz','ret_full.nii.gz','ret_periphfov.nii.gz'};
for i=1:6 ; mean_nii.img = squeeze(mean(allcorrs(:,:,:,i,:),5)) ; save_untouch_nii(mean_nii,titles{i}) ; end

names = {'bot_right','bot_left','top_right','top_left','left_hemi','right_hemi','top_hemi','bot_hemi','full','fov','periph'}; 
cd E:\orientation_retinotopy\fatlas; 
clear fimg; 
thresh = 0.45; 

fimg(:,:,:,1) = squeeze(mean(allcorrs(:,:,:,1,:)>thresh,5));
fimg(:,:,:,2) = squeeze(mean(allcorrs(:,:,:,1,:)<-thresh,5));
fimg(:,:,:,3) = squeeze(mean(allcorrs(:,:,:,2,:)<-thresh,5));
fimg(:,:,:,4) = squeeze(mean(allcorrs(:,:,:,2,:)>thresh,5));
fimg(:,:,:,5) = squeeze(mean(allcorrs(:,:,:,3,:)>thresh,5));
fimg(:,:,:,6) = squeeze(mean(allcorrs(:,:,:,3,:)<-thresh,5));
fimg(:,:,:,7) = squeeze(mean(allcorrs(:,:,:,4,:)>thresh,5));
fimg(:,:,:,8) = squeeze(mean(allcorrs(:,:,:,4,:)<-thresh,5));
fimg(:,:,:,9) = squeeze(mean(allcorrs(:,:,:,5,:)<-thresh,5));
fimg(:,:,:,10) = squeeze(mean(allcorrs(:,:,:,6,:)<-thresh,5));
fimg(:,:,:,11) = squeeze(mean(allcorrs(:,:,:,6,:)>thresh,5));
meanimg = load_untouch_nii('mean.nii.gz'); 
for i=1:length(names)
    meanimg.img = fimg(:,:,:,i); save_untouch_nii(meanimg,['avg_',names{i},'.nii.gz']); 
end


subfimg(:,:,:,:,1) = squeeze((allcorrs(:,:,:,1,:)>thresh));
subfimg(:,:,:,:,2) = squeeze((allcorrs(:,:,:,1,:)<-thresh));
subfimg(:,:,:,:,3) = squeeze((allcorrs(:,:,:,2,:)<-thresh));
subfimg(:,:,:,:,4) = squeeze((allcorrs(:,:,:,2,:)>thresh));
subfimg(:,:,:,:,5) = squeeze((allcorrs(:,:,:,3,:)>thresh));
subfimg(:,:,:,:,6) = squeeze((allcorrs(:,:,:,3,:)<-thresh));
subfimg(:,:,:,:,7) = squeeze((allcorrs(:,:,:,4,:)>thresh));
subfimg(:,:,:,:,8) = squeeze((allcorrs(:,:,:,4,:)<-thresh));
subfimg(:,:,:,:,9) = squeeze((allcorrs(:,:,:,5,:)<-thresh));
subfimg(:,:,:,:,10) = squeeze((allcorrs(:,:,:,6,:)<-thresh));
subfimg(:,:,:,:,11) = squeeze((allcorrs(:,:,:,6,:)>thresh));

cd e:/saved ; all_fmri_rets = subfimg ; save('all_fmri_rets','all_fmri_rets'); 


for i=1:14
   bot_right(i) = sum(sum(sum(allcorrs(:,:,:,1,i)>thresh)));
   bot_left(i) = sum(sum(sum(allcorrs(:,:,:,1,i)<-thresh)));
   top_right(i) = sum(sum(sum(allcorrs(:,:,:,2,i)<-thresh)));
   top_left(i) = sum(sum(sum(allcorrs(:,:,:,2,i)>thresh)));
   left_hemi(i) = sum(sum(sum(allcorrs(:,:,:,3,i)>thresh))); 
   right_hemi(i) = sum(sum(sum(allcorrs(:,:,:,3,i)<-thresh))); 
   top_hemi(i) = sum(sum(sum(allcorrs(:,:,:,4,i)>thresh))); 
   bot_hemi(i) = sum(sum(sum(allcorrs(:,:,:,4,i)<-thresh))); 
   full(i) = sum(sum(sum(allcorrs(:,:,:,5,i)<-thresh)));
   fov(i) = sum(sum(sum(allcorrs(:,:,:,6,i)<-thresh)));
   periph(i) = sum(sum(sum(allcorrs(:,:,:,6,i)>thresh)));
end
all_configs = [full;right_hemi+left_hemi;top_hemi+bot_hemi;fov+periph;bot_right+bot_left+top_left+top_right];
subplot(2,2,1)
topbot_hemis = [bot_hemi;top_hemi];
barwitherr(std(topbot_hemis,0,2)/sqrt(10),mean(topbot_hemis,2)); ylabel('# voxels');
set(gca,'XTick',[1,2],'XTickLabel',{'bottom visual field','top visual field'}); 
[h,p,ci,stats] = ttest(topbot_hemis(1,:),topbot_hemis(2,:));
title(['bottom vs top, p= ',num2str(p)]); 
subplot(2,2,2);
topbot_quads = [bot_right;bot_left;top_right;top_left];
barwitherr(std(topbot_quads,0,2)/sqrt(10),mean(topbot_quads,2)); ylabel('# voxels');
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'bottom right','bottom left','top right','top left'}); 
[h,p,ci,stats] = ttest(mean(topbot_quads([1,2],:)),mean(topbot_quads([3,4],:)));
title(['quadrants, p= ',num2str(p)]); 
subplot(2,2,3);
leftright_hemis = [left_hemi;right_hemi];
barwitherr(std(topbot_hemis,0,2)/sqrt(10),mean(topbot_hemis,2)); ylabel('# voxels');
set(gca,'XTick',[1,2],'XTickLabel',{'left visual field','right visual field'}); 
[h,p,ci,stats] = ttest(leftright_hemis(1,:),leftright_hemis(2,:));
title(['left vs right, p= ',num2str(p)]); 
subplot(2,2,4); 
periphfov_hemis = [fov;periph];
barwitherr(std(periphfov_hemis,0,2)/sqrt(10),mean(periphfov_hemis,2)); ylabel('# voxels');
set(gca,'XTick',[1,2],'XTickLabel',{'fovea','periphery'}); 
[h,p,ci,stats] = ttest(periphfov_hemis(1,:),periphfov_hemis(2,:));
title(['fovea vs periphery, p= ',num2str(p)]); 

barwitherr(std(all_configs,0,2)/sqrt(10),mean(all_configs,2)); ylabel('# voxels'); 

%{

subplot(2,2,1);
shadedErrorBar(or_axis,squeeze(mean(mean(all_ors,1),2)),squeeze(std(mean(all_ors,2),0,1))/sqrt(4),'k'); 
vline([1,180,360],'r'); vline([90,270],'b');

horizontal = [85:95,267:275];
vertical = [1:5,175:185,355:360];
oblique = [40:50,130:140,220:230,310:322]; 
card = [horizontal,vertical]; 


for i=1:length(horizontal); ind_diffs = abs(horizontal(i) - or_axis); horiz_inds(i) = find(ind_diffs==min(ind_diffs)); end; horiz_inds = unique(horiz_inds); 
for i=1:length(vertical); ind_diffs = abs(vertical(i) - or_axis); vert_inds(i) = find(ind_diffs==min(ind_diffs)); end ; vert_inds = unique(vert_inds); 

for i=1:length(oblique); ind_diffs = abs(oblique(i) - or_axis); ob_inds(i) = find(ind_diffs==min(ind_diffs)); end ; ob_inds = unique(ob_inds); 
for i=1:length(card); ind_diffs = abs(card(i) - or_axis); card_inds(i) = find(ind_diffs==min(ind_diffs)); end ; card_inds = unique(card_inds); 

subplot(2,2,2);
shadedErrorBar(or_axis,squeeze(mean(mean(all_ors,1),2)),squeeze(std(mean(all_ors,2),0,1))/sqrt(4),'k'); 
vline(or_axis(horiz_inds),'r'); vline(or_axis(vert_inds),'b'); xlabel('orientation (deg)') ;ylabel('BOLD intensity (A.U)'); 

subplot(2,2,3);
horiz_vert = [squeeze(mean(mean(all_ors(:,:,horiz_inds),2),3)),squeeze(mean(mean(all_ors(:,:,vert_inds),2),3))] + 9;
barwitherr(squeeze(std(horiz_vert,0,1))/sqrt(4),squeeze(mean(horiz_vert,1))); set(gca,'XTickLabel',{'horizontal orientation','vertical orientation'});
[h,p,ci,stats] = ttest(horiz_vert(:,1),horiz_vert(:,2)); 
title(['p=',num2str(p)]); ylabel('BOLD intensity (A.U)');

subplot(2,2,4); 
oblique_card = [squeeze(mean(mean(all_ors(:,:,ob_inds),2),3)),squeeze(mean(mean(all_ors(:,:,card_inds),2),3))] + 9;
barwitherr(squeeze(std(oblique_card,0,1))/sqrt(4),squeeze(mean(oblique_card,1))); set(gca,'XTickLabel',{'oblique orientation','cardinal orientation'});
[h,p,ci,stats] = ttest(oblique_card(:,1),oblique_card(:,2)); 
title(['p=',num2str(p)]); ylabel('BOLD intensity (A.U)');


% quads: br, bl, tl, tr, rep 10s*12, break, rep 10s*12
% allcorrs = [80,80,42,6,5] = quads bot, quads top, hemis leftright, hemis top bot, full, fovperiph
thresh = 0.4; 
for i=1:4
   bot_right(i) = sum(sum(sum(allcorrs(:,:,:,1,i))>thresh));
   bot_left(i) = sum(sum(sum(allcorrs(:,:,:,1,i))<-thresh));
   top_right(i) = sum(sum(sum(allcorrs(:,:,:,2,i))>thresh));
   top_left(i) = sum(sum(sum(allcorrs(:,:,:,2,i))<-thresh));
   top_hemi(i) = sum(sum(sum(allcorrs(:,:,:,4,i)>thresh))); 
   bot_hemi(i) = sum(sum(sum(allcorrs(:,:,:,4,i)<-thresh))); 
   left_hemi(i) = sum(sum(sum(abs(allcorrs(:,:,:,2,i)))>thresh)); 
   right_hemi(i) = sum(sum(sum(abs(allcorrs(:,:,:,1,i))>thresh))); 
   full(i) = sum(sum(sum(abs(allcorrs(:,:,:,1,i))>thresh)));
   fov(i) = sum(sum(sum(abs(allcorrs(:,:,:,1,i))>thresh)));
   periph(i) = sum(sum(sum(abs(allcorrs(:,:,:,1,i))>thresh)));
end
topbot_hemis = [bot_hemi;top_hemi];
barwitherr(std(topbot_hemis,0,2)/sqrt(5),mean(topbot_hemis,2)); ylabel('# voxels');
set(gca,'XTick',[1,2],'XTickLabel',{'bottom visual field','top visual field'}); 
[h,p,ci,stats] = ttest(topbot_hemis(1,:),topbot_hemis(2,:));
title(['p= ',num2str(p)]); 

quads = [bot_right;bot_left;top_right;top_left]; 


m_topbot = squeeze(mean(allcorrs(:,:,:,4,:),5)); 
for i=1:42 ; subplot(4,11,i) ; imagesc(imrotate(squeeze(m_topbot(i+15,:,:)),90),[-.7,.7]) ; colormap jet; end
cd E:\orientation_retinotopy\sub_russell
f1 = load_untouch_nii('f1.nii.gz') ;
f1.img = m_topbot; save_untouch_nii(f1,'m_topbot.nii.gz'); 

%}

% do some kind of spatial comparison...
%{
meanors = squeeze(mean(ors_i(2:end,:,:),1)); 
res_meanors = imresize(meanors,[size(meanors,1),360]); 
res_zbrain = zeros(size(res_ors1,1),360); 
res_zbrain(maskinds,:) = res_meanors; 
zbrain = reshape(res_zbrain,[size(ors1.img,1),size(ors1.img,2),size(ors1.img,3),360]); 
f360 = load_untouch_nii('f360.nii.gz'); 
f360.img = zbrain ; save_untouch_nii(f360,'ors.nii.gz'); 
%}

% contrast analysis:
%{
cont1 = load_untouch_nii('bp_mc_contrast.nii.gz'); 
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

comps = bwconncomp(conved_ts>1); 
pixlist = comps.PixelIdxList; 
for i=1:length(pixlist)
   pix_i = pixlist{i};
   endpix = pix_i(end); 
   newbrain(:,i) = mean(res_cont1(:,pix_i),2) - mean(res_cont1(:,endpix+3:endpix+7),2);  
    
    
end

res_newbrain = reshape(newbrain,[size(cont1.img,1),size(cont1.img,2),size(cont1.img,3),25]); 
cont_newbrain(:,:,:,1) = mean(res_newbrain(:,:,:,1:5),4); 
cont_newbrain(:,:,:,2) = mean(res_newbrain(:,:,:,6:10),4); 
cont_newbrain(:,:,:,3) = mean(res_newbrain(:,:,:,11:15),4); 
cont_newbrain(:,:,:,4) = mean(res_newbrain(:,:,:,16:20),4); 
cont_newbrain(:,:,:,5) = mean(res_newbrain(:,:,:,21:25),4); 
refcont = load_untouch_nii('fcont.nii.gz');
refcont.img = cont_newbrain ;save_untouch_nii(refcont,'cont_int.nii.gz'); 
alpha = [-1.5,-1.25,-1.25,-1.4,-1.75];

f1 = load_untouch_nii('f1.nii.gz'); 
f1.img = medfilt3(voxcorr(refcont.img,alpha)); 
save_untouch_nii(f1,'cont_corrs.nii.gz'); 

%}


