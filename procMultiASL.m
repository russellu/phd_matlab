clear all ; close all ; 
cd c:/shared/asl2 ; ls 
subjs = {
    'nifti_gina_asl_dicom'   
    'nifti_vincent_asl_dicom'
    'nifti_russell2_asl_dicom' 
    'nifti_julie_asl_dicom' 
    'nifti_jeremie_asl_dicom'
    'nifti_alex_asl_dicom'    
    'nifti_genevieve_asl_dicom'
    'nifti_mingham_asl_dicom' 
} ;
stimtypes = [1,2,6] ; 
sub = 1 ; 
for sub=1:max(size(subjs)) ;
    cd(['c:/shared/asl2/',subjs{sub}]) ;
    asl=dir('conc*') ; 
    bold=dir('reg_mc_bold*') ; 
    
    % get the ASL
    disp(['loading asl images for ',subjs{sub}]) ;
    clear alld
    for i=1:max(size(asl))
        aslnii = load_untouch_nii(asl(i).name) ; 
        aslimg = aslnii.img ; 
        dasl = aslimg(:,:,:,1:60)-aslimg(:,:,:,61:120) ;
        alld(i,:,:,:,:) = dasl ; 
    end
    % get the BOLD
    disp(['loading bold images for ',subjs{sub}]) ;
    clear allb
    for i=1:max(size(bold))
        boldnii = load_untouch_nii(bold(i).name) ;  
        boldimg = boldnii.img ; 
        allb(i,:,:,:,:) = boldimg ; 
    end    
    %%% get the triggers and canonical HRF convolved ideal time series
    disp(['calculating correlations for ',subjs{sub}]) ;
    boldtrigs = [2,4,6] ; asltrigs = [1,3,5] ; 
    bold_hrf = spm_hrf(2) ; % spm hrf with a tr of 2 seconds
    asl_hrf = spm_hrf(8) ; 
    stimTimes = dir('stimTimes*') ; clear times types stimvols ideal 
    for s=1:max(size(stimTimes)) ; 
        st = load(stimTimes(s).name) ;  
        st = st.stimTimes ; 
        for i=1:max(size(st)) ;
            times(s,i) = st{i}(1) ;
            types(s,i) = st{i}(2) ; 
        end
        boldtrvols = ceil(times(1,:)./2) ; allboldtrvols = ceil(times./2) ; 
        asltrvols = ceil(times(1,:)./8) ; allasltrvols = ceil(times./8) ; 
        
        % find the indices closest to the onset of the stimulus
        % ************ BOLD *************
        ideal = zeros(1,240) ; 
        for i=1:max(size(boldtrvols)) % fill the ideal with 1s where there was a stimulus
            ideal(boldtrvols(i):boldtrvols(i)+12) = 1 ; 
        end
        boldidealhrf = conv(ideal,bold_hrf) ; % convolve ideal with HRF    
        boldidealhrf = boldidealhrf(1:240) ; 
        
        
         % find the indices closest to the onset of the stimulus 
         % ********* ASL **********  
        ideal = zeros(1,60) ; 
        for i=1:max(size(asltrvols)) % fill the ideal with 1s where there was a stimulus
            ideal(asltrvols(i):asltrvols(i)+3) = 1 ; 
        end
        aslidealhrf = conv(ideal,asl_hrf) ; % convolve ideal with HRF    
        aslidealhrf = aslidealhrf(1:60) ;        
    end
    
    % correlate the ideal HRF with the asl time course in each voxel
    % ******** ASL ************
    clear aslcorrvols
    for i=1:size(alld,1)  
        aslcorrvols(i,:,:,:) = voxcorr(squeeze(alld(i,:,:,:,:)),aslidealhrf) ; 
    end
    postmask = zeros(size(aslcorrvols)) ; postmask(:,1:end,1:30,1:end) = 1 ; 
    aslcorrvols = aslcorrvols.*postmask ; 
    
    % *********** bold corrvols *************
    clear boldcorrvols
    for i=1:size(allb,1)  
        boldcorrvols(i,:,:,:) = voxcorr(squeeze(allb(i,:,:,:,:)),boldidealhrf) ; 
    end
    postmask = zeros(size(boldcorrvols)) ; postmask(:,1:end,1:30,1:end) = 1 ; 
    boldcorrvols = boldcorrvols.*postmask ; 
    
    bothcorrvols = squeeze(mean(aslcorrvols,1)+(mean(boldcorrvols,1)))./2 ; 
    ref = load_untouch_nii('ref.nii.gz') ; 
    allrefs{sub} = ref.img ; 
    allcorrvols{sub} = bothcorrvols ;
    ref.img = bothcorrvols ; save_untouch_nii(ref,'bothcorrvols.nii.gz') ; 
    ref.img = squeeze(mean(aslcorrvols,1)) ; save_untouch_nii(ref,'aslcorrvols.nii.gz') ; 
    ref.img = squeeze(mean(boldcorrvols,1)) ; save_untouch_nii(ref,'boldcorrvols.nii.gz') ; 

    corrthresh = .35 ;
   
    clear aslepochs 
    for sc=1:size(alld,1) %  for all scans
        inds = find(bothcorrvols>corrthresh) ; 
        [i1,i2,i3] = ind2sub(size(bothcorrvols),inds) ; 
        clear aslcorrvoxels
        for i=1:max(size(i1)) ; % get the r values at each index
            aslcorrvoxels(i) = bothcorrvols(i1(i),i2(i),i3(i)) ; 
            aslcorrts(i,:) = alld(sc,i1(i),i2(i),i3(i),:) ; 
        end     
        [cvrvals,cvsortinds] = sort(aslcorrvoxels,'descend') ; 
        sortvox = aslcorrts(cvsortinds,:) ; 
        % get the stimulus indices
        asltrigvols = allasltrvols(asltrigs(sc),:) ; % stimulus onset for that scan
        asltrigtypes = types(asltrigs(sc),:) ; 
        for i=1:max(size(stimtypes))
           asltypinds = asltrigvols(find(asltrigtypes==stimtypes(i))) ;  
           for j=1:max(size(asltypinds))
                aslepochs(i,(sc*3-3)+j,:,:) = sortvox(:,asltypinds(j)-2:asltypinds(j)+5) ;                
           end
        end
    end
    % baseline correct the single voxel epochs
    clear aslnormepochs
    aslepochs = double(aslepochs) ; 
    baseepochs = squeeze(mean(aslepochs(:,:,:,1:2),4)) ; 
    zbases = zscore(baseepochs,0,3) ;    
    for i=1:size(aslepochs,1) 
        for j=1:size(aslepochs,2) ; kcount = 1 ;
            for k=1:size(aslepochs,3)
                baseijk = squeeze(mean(aslepochs(i,j,k,1:2),4)) ; 
                if abs(zbases(i,j,k))<3  % arbitrary %change at 0.1 (veins) and z score < 3 for low baseline
                    eps = (aslepochs(i,j,k,:))./baseijk ; 
                    eps(isnan(eps)) = 0 ; eps(isinf(eps)) = 0 ; 
                    aslnormepochs{i}{j}(kcount,:) = eps ; 
                    kcount = kcount + 1 ; 
                end
            end            
        end
    end
    clear mtrials
    for i=1:max(size(aslnormepochs)) ; 
        for j=1:max(size(aslnormepochs{i})) ; 
            zs = zscore(sum(aslnormepochs{i}{j},2)) ; 
            mtrials(i,j,:) = squeeze(mean(aslnormepochs{i}{j}(abs(zs)<2,:),1)) ;  
        end
    end    
    aslallmtrials(sub,:,:,:) = mtrials ; 

    % correlate the ideal hrf with the BOLD time course, get 3 correlation
    % *********** BOLD ***************
    clear boldcorrvols
    boldcorrvol = zeros(1,1,1,240) ; boldcorrvol(1,1,1,:) = boldidealhrf ; 
    boldcorrvol = repmat(boldcorrvol,[size(allb,2),size(allb,3),size(allb,4),1]) ; 
    for i=1:size(allb,1)  
        boldcorrvols(i,:,:,:) = voxcorr(squeeze(allb(i,:,:,:,:)),boldidealhrf) ; 
    end
    postmask = zeros(size(boldcorrvols)) ; postmask(:,1:end,1:30,1:end) = 1 ; 
    boldcorrvols = boldcorrvols.*postmask ; 
    
    clear boldepochs 
    for sc=1:size(allb,1) %  for all scans
        inds = find(bothcorrvols>corrthresh) ; 
        [i1,i2,i3] = ind2sub(size(bothcorrvols),inds) ; 
        clear boldcorrvoxels
        for i=1:max(size(i1)) ; % get the r values at each index
            boldcorrvoxels(i) = bothcorrvols(i1(i),i2(i),i3(i)) ; 
            boldcorrts(i,:) = allb(sc,i1(i),i2(i),i3(i),:) ; 
        end     
        [cvrvals,cvsortinds] = sort(boldcorrvoxels,'descend') ; 
        sortvox = boldcorrts(cvsortinds,:) ; 
        % get the stimulus indices
        boldtrigvols = allboldtrvols(boldtrigs(sc),:) ; % stimulus onset for that scan
        boldtrigtypes = types(boldtrigs(sc),:) ; 
        for i=1:max(size(stimtypes))
           boldtypinds = boldtrigvols(find(boldtrigtypes==stimtypes(i))) ;  
           for j=1:max(size(boldtypinds))
                boldepochs(i,(sc*3-3)+j,:,:) = sortvox(:,boldtypinds(j)-8:boldtypinds(j)+20) ;                
           end
        end
    end
    % baseline correct the single voxel epochs
    clear normepochs
    boldepochs = double(boldepochs) ; 
    baseepochs = squeeze(mean(boldepochs(:,:,:,4:8),4)) ; 
    zbases = zscore(baseepochs,0,3) ;    
    for i=1:size(boldepochs,1) 
        for j=1:size(boldepochs,2) ; kcount = 1 ;
            for k=1:size(boldepochs,3)
                baseijk = squeeze(mean(boldepochs(i,j,k,4:8),4)) ; 
                basenorm = (boldepochs(i,j,k,:)-baseijk)./baseijk ; 
                if abs(zbases(i,j,k))<3 & abs(max(basenorm)) < .1 % arbitrary %change at 0.1 (veins) and z score < 3 for low baseline
                    boldnormepochs{i}{j}(kcount,:) = (boldepochs(i,j,k,:)-baseijk)./baseijk ; 
                    kcount = kcount + 1 ; 
                end
            end            
        end
    end
    clear mtrials
    for i=1:max(size(boldnormepochs)) ; 
        for j=1:max(size(boldnormepochs{i})) ; 
            mtrials(i,j,:) = squeeze(mean(boldnormepochs{i}{j},1)) ;  
        end
    end    
    boldallmtrials(sub,:,:,:) = mtrials ; 
end
%{
msubs = squeeze(mean(boldallmtrials,3)) ; 
boldt = 9:23 ; 
subplot(2,2,1),errorbar(squeeze(mean(msubs(:,[1,2,3],:),1))',squeeze(std(msubs(:,[1,2,3],:),0,1))'./sqrt(6),'LineWidth',2) ; ylim([-.01,.04]) ; hline(0,'k') ; vline(9,'k') ; vline(21,'k') ; title('BOLD time series') ;
set(gca,'XTick',1:2:29,'XTickLabel',-16:4:40) ; xlabel('time(s)') ; ylabel('(task-rest)/rest') ; legend({'unperturbed','5%contrast','60%randomized'}) ; 
subplot(2,2,2),barwitherr(squeeze(std(mean(msubs(:,:,boldt),3),0,1))'./sqrt(6),squeeze(mean(mean(msubs(:,:,boldt),3),1))') ; ylim([0,.05]) ; set(gca,'XTick',1:3,'XTickLabel',{'unperturbed','5%cont','60%rnd'}) ; title('BOLD mean (t=6:14s)') 
text(.8,squeeze(mean(mean(msubs(:,1,boldt),1),3))+ .0035,['%change=',num2str(squeeze(mean(mean(msubs(:,1,boldt),1),3))*100)]) ;
text(1.8,squeeze(mean(mean(msubs(:,2,boldt),1),3))+ .0035,['%change=',num2str(squeeze(mean(mean(msubs(:,2,boldt),1),3))*100)]) ;
text(2.8,squeeze(mean(mean(msubs(:,3,boldt),1),3))+ .0035,['%change=',num2str(squeeze(mean(mean(msubs(:,3,boldt),1),3))*100)]) ;
bolds(1) = squeeze(mean(mean(msubs(:,1,boldt),1),3)) ; bolds(2) = squeeze(mean(mean(msubs(:,2,boldt),1),3)) ; bolds(3) = squeeze(mean(mean(msubs(:,3,boldt),1),3)) ; 

msubs = squeeze(mean(aslallmtrials,3)) ; 
subplot(2,2,3),errorbar(squeeze(mean(msubs,1))',squeeze(std(msubs,0,1))'./sqrt(6),'LineWidth',2) ; vline(3,'k') ; vline(6,'k') ; title('ASL time series') ; xlim([1,8]) ;
set(gca,'XTick',1:8,'XTickLabel',-16:8:8*6) ; xlabel('time(s)') ; ylabel('task/rest') ; hline(1,'k') ; 
subplot(2,2,4),barwitherr(squeeze(std(mean(msubs(:,:,3:5),3),0,1))'./sqrt(6),squeeze(mean(mean(msubs(:,:,3:5),3),1))') ; set(gca,'XTick',1:3,'XTickLabel',{'unperturbed','5%cont','60%rnd'}) ; title('ASL mean (t=8:24s)') ; 
text(.8,squeeze(mean(mean(msubs(:,1,3:5),1),3))+ .195,['fractional=',num2str(squeeze(mean(mean(msubs(:,1,3:5),1),3))*100)]) ;
text(1.8,squeeze(mean(mean(msubs(:,2,3:5),1),3))+ .195,['fractional=',num2str(squeeze(mean(mean(msubs(:,2,3:5),1),3))*100)]) ;
text(2.8,squeeze(mean(mean(msubs(:,3,3:5),1),3))+ .195,['fractional=',num2str(squeeze(mean(mean(msubs(:,3,3:5),1),3))*100)]) ;
cbfs(1) = squeeze(mean(mean(msubs(:,1,3:5),1),3)) ;cbfs(2) = squeeze(mean(mean(msubs(:,2,3:5),1),3)) ;cbfs(3) = squeeze(mean(mean(msubs(:,3,3:5),1),3)) ;
%}

% plot the single subjects
% BOLD
figure,
for i=1:8 ; 
    subplot(3,3,i) ; errorbar(squeeze(mean(boldallmtrials(i,:,:,:),3))',squeeze(std(boldallmtrials(i,:,:,:),0,3))'/3,'LineWidth',1) ; 
    xlim([1,29]) ; ylim([-.01,.04]); set(gca,'XTick',1:2:29,'XTickLabel',-16:4:40) ; xlabel('time (seconds relative to stimulus onset)') ; ylabel('(task-rest)/rest') ;
    hline(0,'k') ; vline([9,21],'r') ; title(['subject ',num2str(i)]) ; 
end


subplot(2,2,1) ; 
errorbar(squeeze(mean(mean(boldallmtrials(:,1,:,:),3),1))',squeeze(std(mean(boldallmtrials(:,1,:,:),3),0,1))'/sqrt(8),'r','LineWidth',1) ; hold on ; 
errorbar(squeeze(mean(mean(boldallmtrials(:,2,:,:),3),1))',squeeze(std(mean(boldallmtrials(:,2,:,:),3),0,1))'/sqrt(8),'g','LineWidth',1) ; 
errorbar(squeeze(mean(mean(boldallmtrials(:,3,:,:),3),1))',squeeze(std(mean(boldallmtrials(:,3,:,:),3),0,1))'/sqrt(8),'b','LineWidth',1) ; 

xlim([1,29]) ; ylim([-.01,.04]) ;set(gca,'XTick',1:2:29,'XTickLabel',-16:4:40) ; xlabel('time (seconds relative to stimulus onset)') ; ylabel('(task-rest)/rest') ;
hline(0,'k') ; vline([9,21],'k') ; title('grand average BOLD') ; 
subplot(2,2,2) ; 
errorbar(squeeze(mean(mean(aslallmtrials(:,1,:,:),3),1))',squeeze(std(mean(aslallmtrials(:,1,:,:),3),0,1))'/sqrt(8),'r','LineWidth',1) ; hold on ; 
errorbar(squeeze(mean(mean(aslallmtrials(:,2,:,:),3),1))',squeeze(std(mean(aslallmtrials(:,2,:,:),3),0,1))'/sqrt(8),'g','LineWidth',1) ; 
errorbar(squeeze(mean(mean(aslallmtrials(:,3,:,:),3),1))',squeeze(std(mean(aslallmtrials(:,3,:,:),3),0,1))'/sqrt(8),'b','LineWidth',1) ; 
ylim([0.7,2]) ; xlim([1,8]) ; vline(3,'k') ; vline(6,'k') ;  
set(gca,'XTick',1:8,'XTickLabel',-16:8:8*6) ; xlabel('time (seconds relative to stimulus onset)') ; ylabel('task/rest') ; hline(1,'k') ;
title('grand average CBF') ; 
subplot(2,2,3) ; 
clear asls bolds 
bolds = boldallmtrials ; 
for i=1:size(aslallmtrials,1) ; for j=1:size(aslallmtrials,2) ; asls(i,j,:,:) = imresize(squeeze(aslallmtrials(i,j,:,:)),[9,29]) ; end ; end
bolds = squeeze(mean(bolds,3)) ; asls = squeeze(mean(asls,3)) ;
alpha = 0.25 ;  
beta = 1.5 ; 
M=0.1 ; 
cmro2 = ((1-bolds./M).^(1/beta)).*(asls.^(1-alpha/beta)) ; 
errorbar(squeeze(mean(cmro2(:,1,:),1))',squeeze(std(cmro2(:,1,:),0,1))'/sqrt(8),'r','LineWidth',2) ; hold on ; 
errorbar(squeeze(mean(cmro2(:,2,:),1))',squeeze(std(cmro2(:,2,:),0,1))'/sqrt(8),'g','LineWidth',2) ;
errorbar(squeeze(mean(cmro2(:,3,:),1))',squeeze(std(cmro2(:,3,:),0,1))'/sqrt(8),'b','LineWidth',2) ;
ylim([.75,1.5]) ; xlim([1,29]) ;  hline(1,'k') ;  vline([9,21],'r') ; title('modeled CMRO2') ; 
set(gca,'XTick',1:2:29,'XTickLabel',-16:4:40) ;hline(0,'k') ; vline([9,21],'k') ; ylabel('CMRO2 fractional change') ; xlabel('time (seconds relative to stimulus onset') ; 
legend({'unperturbed','5%contrast','60%random'}) ; 



mcount = 1 ; 
for M=.1:.01:.2
clear asls bolds 
subplot(3,4,mcount) ; mcount = mcount + 1 ; 
bolds = boldallmtrials ; 
for i=1:size(aslallmtrials,1) ; for j=1:size(aslallmtrials,2) ; asls(i,j,:,:) = imresize(squeeze(aslallmtrials(i,j,:,:)),[9,29]) ; end ; end
bolds = squeeze(mean(bolds,3)) ; asls = squeeze(mean(asls,3)) ;
alpha = 0.25 ;  
beta = 1.5 ; 
%M=0.1 ; 
cmro2 = ((1-bolds./M).^(1/beta)).*(asls.^(1-alpha/beta)) ; 
errorbar(squeeze(mean(cmro2,1))',squeeze(std(cmro2,0,1))'/sqrt(8)) ;ylim([.75,1.5]) ; xlim([1,29]) ;  hline(1,'k') ;  vline([9,21],'r') ; title(['M = ',num2str(M)]) ; 
end ; suptitle('cmro2 time course for various m-values') ; 
%%% get the cmro2
%{
alpha = 0.2 ;  
beta = 1.3 ; 
M=0.1 ; 
cmro2 = ((1-bolds./M).^(1/beta)).*(cbfs.^(1-alpha/beta)) ; 
mcount = 1 ; clear cmro2
btimes = boldt ; atimes = 3:5 ; 
for M=0.01:.01:0.2
masl = squeeze(mean(aslallmtrials,3)) ; mbold = squeeze(mean(boldallmtrials,3)) ; 
for i=1:6 ;
    boldsi = squeeze(mean(mbold(i,:,btimes),3)) ; 
    cbfsi = squeeze(mean(masl(i,:,atimes),3)) ;    
    cmro2(mcount,i,:) = ((1-boldsi./M).^(1/beta)).*(cbfsi.^(1-alpha/beta)) ; 
end
mcount = mcount + 1 ; 
end
figure,
icount =1 ;
for i=1:4:20 ; 
    subplot(3,3,icount) ; barwitherr(squeeze(std(cmro2(i,:,:),0,2))./sqrt(6),squeeze(mean(cmro2(i,:,:),2))) ; 
    %title(['M = ', num2str(0.01*i),', p(unpt > 5%cont) = ',num2str(anova1(real(squeeze(cmro2(i,:,:))),[],'off'))]); 
    title('M=0.1') ; 
    icount = icount + 1 ; set(gca,'XTickLabel',{'unpt','5%cont','60%rnd'}) ; ylabel('CMRO2 fractional change') ; 
end ; suptitle('modeled CMRO2 for 5 different M values') ; 


alpha = 0.2 ;  
beta = 1.3 ; 
M=0.1 ; 
%%%%% get the cmro2 time course
mcount = 1 ; clear cmro2t
btimes = boldt ; atimes = 3:5 ; 
for M=0.01:.01:0.2
for i=1:6 ;
    masl = squeeze(mean(aslallmtrials(i,:,:,:),3)) ; masl = imresize(masl,[3,29]) ; 
    mbold = squeeze(mean(boldallmtrials(i,:,:,:),3)) ; 
    for j=1:size(mbold,3)
        boldsi = squeeze(mean(mbold(i,:,j),3)) ; 
        cbfsi = squeeze(mean(masl(i,:,j),3)) ;    
        cmro2t(mcount,i,:,j) = ((1-boldsi./M).^(1/beta)).*(cbfsi.^(1-alpha/beta)) ; 
    end
end
mcount = mcount + 1 ; 
end

%%%% plot some single trial BOLD
msubs = squeeze(mean(boldallmtrials,3)) ; hold on ; 
errorbar(squeeze(mean(msubs(:,[1,2,3],:),1))',squeeze(std(msubs(:,[1,2,3],:),0,1))'./sqrt(6),'LineWidth',2) ; ylim([-.01,.055]) ; hline(0,'k') ; vline(9,'k') ; vline(21,'k') ; title('BOLD time series') ;
set(gca,'XTick',1:4:29,'XTickLabel',-16:8:40) ; xlabel('time(s)') ; ylabel('(task-rest)/rest') ; xlim([1,29]) ; legend({'high contrast','low contrast','RANDOMIZATION'}) ; 
subplot(2,2,2),barwitherr(squeeze(std(mean(msubs(:,:,13:18),3),0,1))'./sqrt(6),squeeze(mean(mean(msubs(:,:,13:18),3),1))') ; ylim([0,.05]) ; set(gca,'XTick',1:3,'XTickLabel',{'unperturbed','5%cont','60%rnd'}) ; title('BOLD mean (t=6:14s)') 
text(.8,squeeze(mean(mean(msubs(:,1,13:18),1),3))+ .0035,['%change=',num2str(squeeze(mean(mean(msubs(:,1,13:18),1),3))*100)]) ;
text(1.8,squeeze(mean(mean(msubs(:,2,13:18),1),3))+ .0035,['%change=',num2str(squeeze(mean(mean(msubs(:,2,13:18),1),3))*100)]) ;
text(2.8,squeeze(mean(mean(msubs(:,3,13:18),1),3))+ .0035,['%change=',num2str(squeeze(mean(mean(msubs(:,3,13:18),1),3))*100)]) ;
bolds(1) = squeeze(mean(mean(msubs(:,1,13:18),1),3)) ; bolds(2) = squeeze(mean(mean(msubs(:,2,13:18),1),3)) ; bolds(3) = squeeze(mean(mean(msubs(:,3,13:18),1),3)) ; 

%%%%%%%%% ASL
msubs = squeeze(mean(aslallmtrials,3)) ; 
subplot(1,1,1),errorbar(squeeze(mean(msubs(:,[1,2],:),1))',squeeze(std(msubs(:,[1,2],:),0,1))'./sqrt(6),'LineWidth',2) ; vline(2,'k') ; vline(5,'k') ; title('ASL time series') ; xlim([1,7]) ;
set(gca,'XTick',1:7,'XTickLabel',-8:8:8*6) ; xlabel('time(s)') ; ylabel('task/rest') ; 

a = squeeze(mcmr(:,:,1)-mcmr(:,:,2)) ;
b = squeeze(mcmr(:,:,1)-mcmr(:,:,3)) ;

a = real(reshape(a,[1,191*6])) ;
b = real(reshape(b,[1,191*6])) ; 
subplot(2,2,1) ; hist(a,20) ; title('unperturbed-low contrast cmro2') ; 
subplot(2,2,2) ; hist(b,20) ; title('unperturbed - randomization cmro2') ; 
suptitle(['M = 0.01-0.2, increments of 0.001, subjects = 6']) ;
%}

