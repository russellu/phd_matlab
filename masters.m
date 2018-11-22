clear all ; close all  ;
cd c:/shared/papsaves ; times = load('times') ; times = times.times ; freqs = load('freqs') ; freqs = freqs.freqs ; 

% fmri 
allstimepochs = load('allstimepochs') ; allstimepochs = allstimepochs.allstimepochs ; 
for i=1:22 
    for j=1:6 ;
        for k=1:45
            allstimepochs(i,j,k,:) = allstimepochs(i,j,k,:) - squeeze(allstimepochs(i,j,k,2)) ; 
        end
    end
end

%%% EEG
allersp = load('allersp.mat') ; allersp = allersp.allersp ; 
times = load('times') ; times = times.times ; freqs = load('freqs') ; freqs = freqs.freqs ; 
stims = [2,3,1,5,6] ;
%%% baseline correct
allersp = allersp - repmat(mean(allersp(:,:,:,:,:,times<0),6),[1,1,1,1,1,200]) ; 

% movie for single subjects: 3 stimulus types
subj = squeeze(allersp(1,:,:,:,:,:)) ; 
avg1 = zeros(60,200) ; avg2 = zeros(60,200) ; avg3 = zeros(60,200) ; 
f = figure; set(f,'Position',[100,100,1000,800]) ;
for i=1:135 ;    
    avg1 = avg1 + squeeze(mean(subj(1,1:4,i,:,:),2)) ; 
    avg2 = avg2 + squeeze(mean(subj(6,1:4,i,:,:),2)) ; 
    avg3 = avg3 + squeeze(mean(subj(2,1:4,i,:,:),2)) ; 
    subplot(2,3,2) ; 
    imagesc(times,freqs,avg1/i,[-4,4]) ; axis xy ;xlabel('time(s)') ; ylabel('frequency(hz)') ; vline([0,2],'k') ; title('unperturbed') ; 
    subplot(2,3,3) ; 
    imagesc(times,freqs,avg2/i,[-4,4]) ; axis xy ;vline([0,2],'k') ; title('60% random') ;   
    subplot(2,3,1) ; 
    imagesc(times,freqs,avg3/i,[-4,4]) ; axis xy ;vline([0,2],'k') ; title('5% contrast') ;       
    suptitle(['#trials = ',num2str(i),' SNR = ',num2str(sqrt(i))]) ;
    getframe(f) ; pause(1/i) ;     
end

%%% grand averages
% remove bad trials
mersp = squeeze(mean(mean(mean(allersp(:,:,1:4,:,freqs>5 & freqs<110,times>0 & times<2000),3),5),6)) ;
for i=1:22 ; for j=1:6 ; zsij = abs(zscore(squeeze(mersp(i,j,:)))) ; bads(i,j,:) = zsij>=2 ; end ; end
mcersp = squeeze(mean(allersp(:,:,1:4,:,:,:),3)) ; clear mtersp
for i=1:22 ; for j=1:6 ; mtersp(i,j,:,:) = squeeze(mean(mcersp(i,j,bads(i,j,:)==0,:,:),3)) ; end ; end

avg1 = zeros(60,200) ; avg2 = zeros(60,200) ; avg3 = zeros(60,200) ;  
f = figure; set(f,'Position',[100,100,1000,800]) ;
for i=1:22
    avg1 = avg1 + squeeze(mtersp(i,1,:,:)) ; avg2 = avg2 + squeeze(mtersp(i,2,:,:)) ; avg3 = avg3 + squeeze(mtersp(i,6,:,:)) ; 
    subplot(2,3,2) ; imagesc(times,freqs,avg1/i,[-2,2]) ; axis xy ; vline([0,2],'k') ; xlabel('time(s)') ; ylabel('frequency(hz)') ; title('unperturbed') ; 
    subplot(2,3,1) ; imagesc(times,freqs,avg2/i,[-2,2]) ; axis xy ; vline([0,2],'k') ; xlabel('time(s)') ; ylabel('frequency(hz)') ; title('5%contrast') ; 
    subplot(2,3,3) ; imagesc(times,freqs,avg3/i,[-2,2]) ; axis xy ; vline([0,2],'k') ; xlabel('time(s)') ; ylabel('frequency(hz)') ; title('60%random') ; 
    suptitle(['#subjects = ',num2str(i),' SNR = ',num2str(sqrt(i))]) ;
    getframe ; pause(0.25) ; 
end

cd c:/shared/allres/suhan2 ; ls 
ica = pop_loadset('ica_notch85.set') ; 
eps = pop_epoch(ica,{'S 11'},[-.85,2.85]) ; 
elabs = {eps.chanlocs.labels} ; 
for i=1:size(eps.data,1)
    [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(eps.data(i,:,:)),eps.pnts,[eps.xmin,eps.xmax],eps.srate,0,...
        'plotersp','off','plotitc','off','winsize',64,'freqs',[1,120],'nfreqs',60,'baseline',0) ;
end
mcersp = squeeze(mean(mean(ersp(:,freqs>50 & freqs<80,times>0 & times<2),2),3)) ; 
[sv,si] = sort(mcersp,'descend') ; 
s1 = squeeze(eps.data(34,:,:)) ; 
plot(s1(:,130)) ; xlabel('time(ms)') ; ylabel('microvolts') ; set(gca,'XTick',1:50:length(eps.times),'XTickLabel',round(eps.times(1:50:end))) ; vline([219,732],'k') ; 


for i=1:32 
    h = subplot(2,32,i) ; p = get(h, 'pos') ; 
    imagesc(squeeze(ersp(si(i),:,:)),[-3,3]) ; axis xy ; set(gca,'XTick',[],'YTick',[]) ; title(elabs{si(i)}) ;   
    if i>1
        p(1) = p_prev(1) + 0.02 ; set(h,'pos',p) ; 
    else set(h,'pos',p) ;
    end
    p_prev = p ; 
end
for i=33:64 
    h = subplot(2,32,i) ; p = get(h, 'pos') ; 
    imagesc(squeeze(ersp(si(i),:,:)),[-3,3]) ; axis xy ; set(gca,'XTick',[],'YTick',[]) ; title(elabs{si(i)}) ;   
    if i>33
        p(1) = p_prev(1) + 0.02 ; set(h,'pos',p) ; 
    else set(h,'pos',p) ;
    end
    p_prev = p ; 
end


%%%%%% display the FMRI ica components:
cd c:/shared/allfmris/sub_russell/ica ; ls
mix = load('melodic_mix') ; 
refbrain = load_untouch_nii('0a.nii.gz') ; refimg = refbrain.img ; 
comps = load_untouch_nii('melodic_IC.nii.gz') ; compimg = comps.img ;
% motion
subplot(2,1,1) ; plot(mix(:,1),'LineWidth',2) ; xlabel('time(TR=2s)') ; 
plotoverlayIntensity2D(squeeze(refimg(:,:,20)),squeeze(abs(compimg(:,:,20,1)))>4,squeeze(compimg(:,:,20,1)),270) ;

cd c:/shared/allfmris/sub_gabriella/ica ; ls
mix = load('melodic_mix') ; 
refbrain = load_untouch_nii('0a.nii.gz') ; refimg = refbrain.img ; 
comps = load_untouch_nii('melodic_IC.nii.gz') ; compimg = comps.img ;
% visual
subplot(2,1,1) ; plot(mix(1:120,9),'LineWidth',2) ; vline(3:8:120) ; ylim([-2,3]) ; xlabel('time(TR=2s)') ; 
plotoverlayIntensity2D(squeeze(refimg(:,:,17)),squeeze((compimg(:,:,17,9))>3),squeeze(compimg(:,:,17,9)),270) ;
% eyes
subplot(2,1,1) ; plot(mix(1:120,6),'LineWidth',2) ; vline(3:8:120) ; xlabel('time(TR=2s)') ; 
plotoverlayIntensity2D(squeeze(refimg(:,:,16)),squeeze(abs(compimg(:,:,16,6))>5),squeeze(compimg(:,:,16,6)),270) ;

%weird 1
compn = 34 ; slice = 23 ; 
subplot(2,1,1) ; plot(mix(1:120,compn),'LineWidth',2) ; vline(3:8:120) ; xlabel('time(TR=2s)') ; 
plotoverlayIntensity2D(squeeze(refimg(:,:,slice)),squeeze(abs(compimg(:,:,slice,compn))>3),squeeze(compimg(:,:,slice,compn)),270) ;


%%%%% swi with correlation map overlayed
cd c:/shared/corrvis ; ls 
swi = load_untouch_nii('swi.nii.gz') ; swimg = swi.img ;
corrswi = load_untouch_nii('swicorrs.nii.gz') ; corrimg = corrswi.img ; 

for slice = 40:2:90 ; 
    figure,
    plotoverlayIntensity2D(squeeze(swimg(:,:,slice)),abs(squeeze(corrimg(:,:,slice))),squeeze(corrimg(:,:,slice)),90) ; 
end

%%%% CBF visualization
cd c:/shared/asl2/nifti_russell2_asl_dicom/corrvis ; ls
blood = load_untouch_nii('blood_in_t1.nii.gz') ; bloodim = blood.img ; 
T1 = load_untouch_nii('T1.nii.gz' ) ; t1im = T1.img ; 
for slice = 30:2:70 ; 
    figure,
    plotoverlayIntensity2D(squeeze(t1im(:,:,slice)),mat2gray(squeeze(bloodim(:,:,slice))),squeeze(bloodim(:,:,slice)),270) ; 
end
%%%% mip
cd c:/shared/corrvis ; ls 
act = load_untouch_nii('MONG_01_RB_WIP_ToF_0.4mm_SENSE_4_1.nii') ; 

%%%%%% synchrony figure:
f = figure; set(f,'Position',[100,100,1000,600]) ;
a1 = zeros(50,50,5000) ; 
a2 = zeros(50,50,5000) ; 
a3 = zeros(50,50,5000) ; 
for i=1:50 ;
    for j=1:50 ;
        sincurve = sin((1:.25:5000)) ;
        a1(i,j,:) = sincurve(1:5000) ; 
        sincurve = sin((1:.25:5000)+rand*2) ;
        a2(i,j,:) = sincurve(1:5000) ; 
        sincurve = sin((1:.25:5000)+rand*100) ;
        a3(i,j,:) = sincurve(1:5000) ;           
    end
end
tcurve1 = [] ; tcurve2 = [] ; tcurve3 = [] ;
for i=1:5000  ;
    subplot(2,3,1) ; imagesc(a1(:,:,i),[-1,1]) ; title('perfect synchrony') ; 
    subplot(2,3,4) ; tcurve1(length(tcurve1)+1) = (sum(sum(a1(:,:,i)))) ; plot(tcurve1) ; ylim([-3000,3000]) ; 
    subplot(2,3,2) ; imagesc(a2(:,:,i),[-1,1]) ; title('partial synchrony') ; 
    subplot(2,3,5) ; tcurve2(length(tcurve2)+1) = (sum(sum(a2(:,:,i)))) ; plot(tcurve2) ; ylim([-3000,3000]) ; 
    subplot(2,3,3) ; imagesc(a3(:,:,i),[-1,1]) ; title('zero synchrony') ; 
    subplot(2,3,6) ; tcurve3(length(tcurve3)+1) = (sum(sum(a3(:,:,i)))) ; plot(tcurve3) ; ylim([-3000,3000]) ; 
    getframe ; pause(0.01) ; 
end

cd c:/shared/allt1s/ ; ls 
t1s = dir('*nii') ; 
for i=1:length(t1s) ; a = load_untouch_nii(t1s(i).name) ; allt1s{i} = a.img ; end

for i=1:22 
    h = subplot(1,22,i) ; p = get(h, 'pos') ; 
    imshow(uint8(mat2gray(rot90(squeeze(allt1s{i}(120,:,:))))*255)) ; colormap gray
    if i>1
        p(1) = p_prev(1) + 0.03 ; set(h,'pos',p) ; 
    else set(h,'pos',p) ;
    end
    p_prev = p ; 
end






