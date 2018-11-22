clear all ; close all ; 
subs = {'sub_lyes','sub_amal','sub_lyndis','sub_valerie','sub_reihaneh','sub_samuel','sub_cesar','sub_angelina','sub_russell','sub_esteban','sub_raphael','sub_greg','sub_felix','sub_pascal'};
subdates = {'07/06/2018','07/07/2018','01/06/2018','01/06/2018','17/05/2018','17/05/2018','11/05/2018','11/05/2018','12/04/2018','26/04/2018','26/04/2018','12/04/2018','03/05/2018','17/05/2018'};

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
    
    ors1 = load_untouch_nii('reg_bp_mc_orientation_1.nii'); 
    ors1.img = ors1.img;
    res_ors1 = reshape(ors1.img,[numel(ors1.img(:,:,:,1)),750]); 
    rotangles = mod(1:360*8,360) ;  rot_time = 480/0.68; rotangles = imresize(rotangles,[1,rot_time],'nearest'); 
    ts = zeros(1,750); startt = round(10/0.68); ts(startt:startt+length(rotangles)-1) = rotangles; ts = circshift(ts,round(6/0.68)); 

    maskinds = find(corrs1_fullfov(:) < -.5); 
    ors1_ts = res_ors1(maskinds,:); 
    inds = find(diff(ts)<-200); 
    clear ors_i; 
    for i=1:length(inds)
        ors_i(i,:,:) = ors1_ts(:,inds(i)-87:inds(i)); 
    end
    or_axis = imresize(1:360,[1,88]); 
    %plot(or_axis,squeeze(mean(mean(ors_i(2:end,:,:),1),2))) ; vline([45,135,225,315],'g');
    %vline([1,180,360],'r'); vline([90,270],'k');
    
    all_ors_1(sb,:,:) = squeeze(mean(ors_i,2)); 
    
    % ORSS 2
    ors1 = load_untouch_nii('reg_bp_mc_orientation_2.nii'); 
    ors1.img = ors1.img;
    res_ors1 = reshape(ors1.img,[numel(ors1.img(:,:,:,1)),750]); 
    rotangles = mod(1:360*8,360) ;  rot_time = 480/0.68; rotangles = imresize(rotangles,[1,rot_time],'nearest'); 
    ts = zeros(1,750); startt = round(10/0.68); ts(startt:startt+length(rotangles)-1) = rotangles; ts = circshift(ts,round(6/0.68)); 

    maskinds = find(corrs1_fullfov(:) < -.5); 
    ors1_ts = res_ors1(maskinds,:); 
    inds = find(diff(ts)<-200); 
    clear ors_i; 
    for i=1:length(inds)
        ors_i(i,:,:) = ors1_ts(:,inds(i)-87:inds(i)); 
    end
    or_axis = imresize(1:360,[1,88]); 
    plot(or_axis,squeeze(mean(mean(ors_i(2:end,:,:),1),2))) ; vline([45,135,225,315],'g');
    vline([1,180,360],'r'); vline([90,270],'k');
    
    all_ors_2(sb,:,:) = squeeze(mean(ors_i,2)); 
    
    allfov1(sb,:,:,:) = corrs1_fullfov; 
    allfov2(sb,:,:,:) = corrs2_fullfov; 

    
    horizontal = [85:95,267:275];
    vertical = [1:5,175:185,355:360];
    oblique = [40:50,130:140,220:230,310:322]; 
    card = [horizontal,vertical]; 

    for i=1:length(horizontal); ind_diffs = abs(horizontal(i) - or_axis); horiz_inds(i) = find(ind_diffs==min(ind_diffs)); end; horiz_inds = unique(horiz_inds); 
    for i=1:length(vertical); ind_diffs = abs(vertical(i) - or_axis); vert_inds(i) = find(ind_diffs==min(ind_diffs)); end ; vert_inds = unique(vert_inds); 

    for i=1:length(oblique); ind_diffs = abs(oblique(i) - or_axis); ob_inds(i) = find(ind_diffs==min(ind_diffs)); end ; ob_inds = unique(ob_inds); 
    for i=1:length(card); ind_diffs = abs(card(i) - or_axis); card_inds(i) = find(ind_diffs==min(ind_diffs)); end ; card_inds = unique(card_inds); 

    
    m_orsi = squeeze(mean(ors_i,1)); 
    
    obcount = 0 ; cardcount = 0; horizcount = 0; vertcount = 0; 
    for i=1:size(m_orsi,1)
       if mean(m_orsi(i,ob_inds),2) > mean(m_orsi(i,card_inds),2) ; obcount = obcount + 1; else cardcount = cardcount +  1; end
       if mean(m_orsi(i,horiz_inds),2) > mean(m_orsi(i,vert_inds),2) ; horizcount = horizcount + 1; else vertcount = vertcount +  1; end
    end
    
    allcounts(sb,:) = [obcount,cardcount,horizcount,vertcount]; 
    
end



badsubs = [9:12]; goodsubs = zeros(size(badsubs)) ; goodsubs(badsubs) = 1 ; goodsubs = find(goodsubs==0) ;



subplot(2,1,1); 
shadedErrorBar(or_axis,squeeze(mean(mean(all_ors_1,1),2)),squeeze(std(mean(all_ors_1,2),0,1))/sqrt(14),{'b'}); xlabel('orientation') ; ylabel('BOLD amplitude'); 
hold on ; 
shadedErrorBar(or_axis,squeeze(mean(mean(all_ors_2(goodsubs,:,:),1),2)),squeeze(std(mean(all_ors_2(goodsubs,:,:),2),0,1))/sqrt(12),{'r'}); xlabel('orientation') ; ylabel('BOLD amplitude'); 

plot(squeeze(mean(mean(all_ors_1))),squeeze(mean(mean(all_ors_2))),'kd') ; lsline
vline(horizontal,'r'); vline(vertical,'b'); vline(oblique,'g'); 


plot(squeeze(mean(mean(all_ors_1(:,:,1:44),1),2)) + squeeze(mean(mean(all_ors_1(:,:,45:88),1),2))) ;
plot(squeeze(mean(mean(all_ors_2(:,:,1:44),1),2)) + squeeze(mean(mean(all_ors_2(:,:,45:88),1),2))) ;



m_card = squeeze(mean(mean(all_ors_1(:,:,card_inds),2),3)); 
m_ob = squeeze(mean(mean(all_ors_1(:,:,ob_inds),2),3)); 
m_vert = squeeze(mean(mean(all_ors_1(:,:,vert_inds),2),3)); 
m_horiz = squeeze(mean(mean(all_ors_1(:,:,horiz_inds),2),3)); 
cardobs = [m_card,m_ob]; horizverts = [m_horiz,m_vert];
[~,cardobs_p,~,~] = ttest2(cardobs(:,1),cardobs(:,2)); 
[~,horizverts_p,~,~] = ttest2(horizverts(:,1),horizverts(:,2)); 
cardobs = cardobs + 5; horizverts = horizverts+5; 
subplot(2,2,3);
barwitherr(squeeze(std(cardobs,0,1))/sqrt(5),squeeze(mean(cardobs,1))); set(gca,'XTickLabel',{'cardinal','oblique'});
title(['p=',num2str(cardobs_p)]);
subplot(2,2,4); 
barwitherr(squeeze(std(horizverts,0,1))/sqrt(5),squeeze(mean(horizverts,1))); set(gca,'XTickLabel',{'horizontal','vertical'});
title(['p=',num2str(horizverts_p)]);


cd e:/saved; 
subdb = load('subdb') ; subdb = subdb.subdb; 

mors1 = squeeze(mean(mean(all_ors_1,1),2)); 
mors2 = squeeze(mean(mean(all_ors_2,1),2)); 

resmors1 = imresize(mors1,[360,1]) ;
resmors2 = imresize(mors2,[360,1]) ;

clear ocorrs1 ocorrs2; 
for i=1:5
for j=1:60 
    ocorrs1(i,j) = corr(resmors1,squeeze(mean(mean(subdb(:,i,j,:),1),2))) ; 
    ocorrs2(i,j) = corr(resmors2,squeeze(mean(mean(subdb(:,i,j,:),1),2))) ; 
end
end









