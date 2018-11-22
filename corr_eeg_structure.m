clear all ; close all ; 
subs = {'charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','lisa','marc','marie',...
    'mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ;

clear allersp ; 
for s=1:length(subs)
    cd(['c:/shared/allres/',subs{s}]) ; ls 
    ersp = load('ersp.mat') ; 
    stersp = load('stersp') ; stersp = stersp.ersp ; 
    allstersp(s,:,:,:,:) = squeeze(stersp(:,1:3,:,:)) ;
    ersp = ersp.ersp ; allersp(s,:,:) = squeeze(mean(ersp(1:2,:,:),1)) ; 
  
    cd(['c:/shared/all_white_normals/a1_good/sub_',subs{s}]) ; 
    corrs = load_untouch_nii('noeyecorrs.nii.gz') ; 
    leftnorms = dir('*whitesub*lh*gz') ; rightnorms = dir('*whitesub*rh*gz') ; 
    leftnorms = load_untouch_nii(leftnorms.name) ; rightnorms = load_untouch_nii(rightnorms.name) ; 
    bothnorms = leftnorms.img + rightnorms.img ; 
    threshcount = 1 ; 
    for thresh=0.1:.01:0.5
        corrmask = corrs.img > thresh ; corrmask(:,1:end-120,:) = 0 ; 
        inds = find(corrmask==1) ; 
        resnorms = reshape(bothnorms,[numel(bothnorms(:,:,:,1)),3]) ;  
        resleft = reshape(leftnorms.img,[numel(leftnorms.img(:,:,:,1)),3]) ;  
        resright = reshape(rightnorms.img,[numel(rightnorms.img(:,:,:,1)),3]) ;     
        maskednorms = resnorms(inds,:) ;
        maskedleftnorms = resleft(inds,:) ;
        maskedrightnorms = resright(inds,:) ;
        i0 = 1-norm(sum(maskednorms,1))/length(maskednorms) ;  
        i0left = 1-norm(sum(maskedleftnorms,1))/length(maskedleftnorms) ;  
        i0right = 1-norm(sum(maskedrightnorms,1))/length(maskedrightnorms) ;  
        alli0(s,threshcount) = i0 ; alli0left(s,threshcount) = i0left ; alli0rights(s,threshcount) = i0right ; 
        threshcount = threshcount + 1 ; 
    end
end

mersp = squeeze(mean(mean(mean(allstersp(:,:,:,:,40:180),2),3),5)) ; 
[r,p] = corr(alli0,mersp) ; 
[rr,pr] = corr(alli0rights,mersp) ; 
subplot(1,2,1) ; 
imagesc(1:2:120,0.1:0.01:0.5,r) ; colorbar ; ylabel('ROI correlation threshold (rho)') ; xlabel('frequency (hz)') ; 
hline(.3,'m') ; hline(.45,'k') ; 
subplot(1,2,2) ; 
plot(r(21,:),'m','LineWidth',3) ;hold on ; hline(0,'k') ; plot(r(36,:),'k','LineWidth',3) ; legend({'corr thresh=0.30','corr thresh=0.45'}) ; 
ylabel('correlation (rho)') ; xlabel('frequency(hz)') ; 
set(gca,'XTick',1:5:60,'XTickLabel',1:10:120) 




