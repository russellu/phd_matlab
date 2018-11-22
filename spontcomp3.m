clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie'} ; 

comps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38]} ; % all subjects, right and left
fcomps = {[64,80,12],[41,31,40],[61,48,15],[6,84,79],[79,80,50],[75,77,41]} ; 

for sub=1:length(subs)

    %%%% BOLD FMRI processing:
    cd(['c:/shared/badger_mri/',subs{sub},'/nii/melodic']) ; 
    mix = load('melodic_mix') ; 
    weights = load_untouch_nii('melodic_IC.nii.gz') ; 
    allweights{sub} = weights.img ; 
    meansub = load_untouch_nii('mean.nii.gz') ; 
    allmeans{sub} = meansub.img ; 
    mixinds = 1:735:735*5 ; 
    for i=1:length(mixinds) ; boldinds{i} = mixinds(i):mixinds(i)+734 ; end
    boldinds{6} = 735*5+1:735*5+1+449 ; 
    for i=1:length(boldinds) ; segmix{i} = mix(boldinds{i},:) ; end
    
    %%%% EEG processing:
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'highfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ; gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ; allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    setNames = {allstims1(1).name,allstims2(2).name,gammas1(1).name,gammas2(2).name,movies(1).name,rests(1).name} ;   
    
    finalersp = zeros(3,64,50,200) ; 
    clear allcorrs ; 
    for setN = 1:6;    
        EEG = pop_loadset(setNames{setN}) ; 
        if setN == 1 ; alleegs{sub} = EEG ; end

        %%% lots of noise in the data. need to remove noisy epochs in analyzer,
        %%% and then re-run ICA. especially on valerie
        trigsr = {EEG.urevent.type} ; latsr = cell2mat({EEG.urevent.latency}) ; 
        r128s = find(strcmp(trigsr,'R128')) ;
        rlats = latsr(r128s) ; 
        firstlats(setN) = rlats(1) ; 
        lastlats(setN) = rlats(length(rlats)) ; 
        
        EEG = pop_select(EEG,'nopoint',[1,rlats(1) ; rlats(length(rlats)),size(EEG.data,2)]) ; 
        EEG = eeg_checkset(EEG) ; 

     % get the continuous power in all frequencies
        eegsecs = EEG.pnts./EEG.srate ; 
        eegTRs = round(eegsecs./0.693) ; 
        clear epow
        for c=1:length(comps{sub})  
            [epow(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(EEG.icaact(comps{sub}(c),:),EEG.pnts,[EEG.xmin,EEG.xmax],EEG.srate,0,...
                'plotersp','off','plotitc','off','baseline',NaN,'timesout',eegTRs,'nfreqs',50,'freqs',[1,100],'winsize',EEG.srate) ; 
        end   
        cmix = segmix{setN}(1:eegTRs,:) ;  
        cmix = eegfiltfft(cmix',1/.693,.02,1.5) ; cmix = cmix' ; 
        hrf = spm_hrf(0.693) ; 
        allpows{sub,setN} = epow ; 
        clear convpow
        for i=1:size(epow,1)
            for j=1:size(epow,2)
                conved = conv(squeeze(epow(i,j,:)),hrf,'full') ;
                convpow(i,j,:) = conved(1:eegTRs) ;             
            end
        end  
        for i=1:size(convpow,2)
           allcorrs(setN,i,:,:) = corr(squeeze(convpow(:,i,50:end-50))',squeeze(cmix(50:end-50,:))) ;        
        end
 
           if setN==3 || setN==4
                trigs = {'S  1','S  2','S  3'} ; 
                alltrigs = {EEG.urevent.type} ; 
                lats = cell2mat({EEG.urevent.latency}) ; 
                clear stimlats
                for t=1:length(trigs)
                    tr = find(strcmp(trigs{t},alltrigs)) ;
                    stimlats(t,:) = lats(tr) ; 
                end
                stimlats = (stimlats - firstlats(setN))./EEG.srate ; 
                stimtrs = round(stimlats./0.693) ; 
                fmricomps = segmix{setN} ; 
                fmricomps = eegfiltfft(fmricomps',1/.693,0.02,1.5) ; fmricomps = fmricomps'  ;
                % first, create the ideal time series, convolve, and correlate to
                % isolate components that were modulated by the stimulus presentation
                ideal = zeros(1,735) ; 
                stimduration = round(5./0.693) ; 
                res_stimtrs = reshape(stimtrs,[1,numel(stimtrs)]) ; 
                for i=1:length(res_stimtrs) ; ideal(res_stimtrs(i):res_stimtrs(i)+stimduration) = 1 ; end
                hrf = spm_hrf(0.693) ; 
                convedf = conv(ideal,hrf,'full') ; convedf = convedf(1:735) ;        
                corrs = corr(convedf(50:end-50)',fmricomps(50:end-50,:)) ; 
                subcorrs{sub}(setN-2,:) = corrs ; 
                %figure,
                %plot(mat2gray(convedf(50:end-50))) ; hold on ;
                %plot(mat2gray(squeeze(mean(mean(convpow(:,freqs>40 & freqs<70,50:end-50),1),2))),'r') ; 
                %{
                if setN==3
                subplottight(6,1,sub) ; 
                bar(ideal(300:600),'k') ; hold on ; plot(mat2gray(squeeze(mean(mean(epow(:,freqs>40 & freqs<65,300:600),1),2)))','r','LineWidth',2) ; 
                %subplot(2,1,2) ; 
                plot(mat2gray(squeeze(mean(mean(epow(:,freqs>8 & freqs<25,300:600),1),2)))','b','LineWidth',2) ; 
                end
                %}
           end

            

        clear xc ; 
        for f=1:length(freqs)
            xc(f,:)  = crosscorr(smooth(squeeze(mean(epow(1,freqs>=freqs(f)-2.5 & freqs<=freqs(f)+2.5,50:end-50),2)),10),squeeze(cmix(50:end-50,fcomps{sub}(1))),20) ; 
        end
        allxc(sub,setN,:,:) = xc ; 

    end
    

  %  subplot(2,3,sub) ; imagesc(xc) ; 
    allfinalersp(sub,:,:,:,:) = finalersp ; 
    allmix(sub,:,:) = cmix(:,fcomps{sub}) ; 
    allsubcorrs{sub} = allcorrs ; 
     
  %  figure,imagesc(squeeze(mean(epow,1)))
end

%{
for i=1:length(allsubcorrs)
   %subplot(2,3,i) ; 
   sc(i,:,:) = (squeeze(mean(mean(allsubcorrs{i}(:,:,:,fcomps{i}(1:2)),3),4))') ; 
end
plot(squeeze(mean(mean(sc(:,:,6),1),3)))
%}


%{
for i=1:length(subcorrs) ; 
    [sv,si] = (sort(squeeze(mean(subcorrs{i},1)),'descend')) ;
    allinds(i,:) = si(1:3) ; 
end
%}

%{

labels= {'retinotopy#1','retinotopy#2','gamm#1','gamma#2','movie','rest'} ; 

ftimes = (-20:20)*0.693 ; 
for i=1:6 ; 
    subplot(2,3,i) ; 
   imagesc(ftimes,freqs,squeeze(mean(allxc(:,i,:,:),1)),[-.25,.25]) ; xlabel('lag(s)') ; axis xy
   ylabel('freq(hz)') ; vline(0,'k') ; title(labels{i})  ; 
   %plot(squeeze(mean(mean(allxc(:,i,freqs>40 & freqs<100,:),1),3)),'r','LineWidth',2) ;  hold on ; 
   %plot(squeeze(mean(mean(allxc(:,i,freqs>8 & freqs<20,:),1),3)),'b','LineWidth',2) 
   %vline(20) ; hline(0,'k') ; 
   
end

ijcount = 1 ; 
for i=1:6 ; 
    for j=1:6
    subplot(6,6,ijcount) ;
    imagesc(ftimes,freqs,squeeze(mean(allxc(i,j,:,:),2)),[-.25,.25]) ; axis xy ; vline(0,'k') ; 
    title(['subject=',num2str(i)])  ; 
    ijcount = ijcount + 1; 
    if j==1 && i==6
    ylabel('freq(hz)') ; xlabel('lag(s)') ;
    end
    end
end


sublegend = {'S1','S2','S3','S4','S5','S6'} ; 
subplot(1,2,1) ;
plot(squeeze(mean(mean(allxc(:,5:6,freqs>8 & freqs<40,:),3),2))','LineWidth',2) ; legend(sublegend) ; vline(21,'k') ; 
set(gca,'XTick',1:5:41,'XTickLabel',ftimes(1:5:end)) ; xlabel('lag(s)') ; ylabel('correlation (r)') ; title('8-40Hz') ; 
subplot(1,2,2) ; 
plot(squeeze(mean(mean(allxc(:,5:6,freqs>40 & freqs<100,:),3),2))','LineWidth',2) ; legend(sublegend) ; vline(21,'k') ; 
set(gca,'XTick',1:5:41,'XTickLabel',ftimes(1:5:end)) ; xlabel('lag(s)') ; ylabel('correlation (r)') ;  title('40-100Hz') ; 

sublegend = {'S1','S2','S3','S4','S5','S6'} ; 
subplot(1,2,1) ;
plot(squeeze(mean(mean(allxc(:,3:4,freqs>8 & freqs<40,:),3),2))','LineWidth',2) ; legend(sublegend) ; vline(21,'k') ; 
set(gca,'XTick',1:5:41,'XTickLabel',ftimes(1:5:end)) ; xlabel('lag(s)') ; ylabel('correlation (r)') ; title('8-40Hz') ; 
subplot(1,2,2) ; 
plot(squeeze(mean(mean(allxc(:,3:4,freqs>40 & freqs<100,:),3),2))','LineWidth',2) ; legend(sublegend) ; vline(21,'k') ; 
set(gca,'XTick',1:5:41,'XTickLabel',ftimes(1:5:end)) ; xlabel('lag(s)') ; ylabel('correlation (r)') ;  title('40-100Hz') ; 

freqlabs = {'8-40Hz','40-100Hz'} ; 

subplot(1,2,1) ; 
basehigh = squeeze(mean(mean(mean(allxc(:,5:6,freqs>40,:),1),2),3)) ;
baselow = squeeze(mean(mean(mean(allxc(:,5:6,freqs<40 & freqs> 8,:),1),2),3)) ;
plot(squeeze(mean(mean(mean(allxc(:,5:6,freqs>8 & freqs<40,:),1),2),3)),'LineWidth',2) ; hold on ; 
plot(squeeze(mean(mean(mean(allxc(:,5:6,freqs>40,:),1),2),3)),'r','LineWidth',2) ;
hline(0,'k') ; xlim([0,42]) ; vline(21,'k') ; vline(find(baselow==min(baselow)),'b') ; vline(find(basehigh==max(basehigh)),'r') ; legend(freqlabs) ; 
set(gca,'XTick',1:5:41,'XTickLabel',ftimes(1:5:end)) ; xlabel('lag(s)') ; ylabel('correlation (r)') ; 
title('spontaneous (movie+rest)') ; 

subplot(1,2,2) ; 
taskhigh = squeeze(mean(mean(mean(allxc(:,3:4,freqs>40,:),1),2),3)) ;
tasklow = squeeze(mean(mean(mean(allxc(:,3:4,freqs<40 & freqs> 8,:),1),2),3)) ;
plot(squeeze(mean(mean(mean(allxc(:,3:4,freqs>8 & freqs<40,:),1),2),3)),'LineWidth',2) ; hold on ; 
plot(squeeze(mean(mean(mean(allxc(:,3:4,freqs>40,:),1),2),3)),'r','LineWidth',2) ; vline(find(tasklow==min(tasklow)),'b') ; vline(find(taskhigh==max(taskhigh)),'r') ; hline(0,'k') ; xlim([0,42]) ; vline(21,'k') ; 
legend(freqlabs) ; set(gca,'XTick',1:5:41,'XTickLabel',ftimes(1:5:end)) ; xlabel('lag(s)') ; ylabel('correlation (r)') ; 
title('event-related') ; 





restpow = allpows(:,6) ; clear reeg
for i=1:length(restpow) ; reeg(i,:,:) = mean(restpow{i},1) ; end ; 
mkern = squeeze(mean(allxc(:,6,:,:),1)) ; 
pred = zeros(6,size(reeg,3)) ; 
for i=1:6 ; 
   for j=51:size(reeg,3)
      pred(i,j) = sum(sum(mkern.*squeeze(reeg(i,:,j-20:j+20)))) ;  
       
   end
    
end

%%% use a random sampling approach to estimate the transfer function
inds = 50:400 ; 
for i=1:length(inds)
    
    
end

%}







subplot(1,2,1) ; imagesc(squeeze(mean(allxc(:,5,:,:),1)),[-0.35,.35]) ; axis xy ; 
subplot(1,2,2) ; imagesc(squeeze(mean(allxc(:,6,:,:),1)),[-0.35,.35]) ; axis xy ; 



