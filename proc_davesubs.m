cd c:/shared/raw_eeg/MONG_01_RB ; close all ; clear all ; 
ls 
mongs = dir('*vhdr') ; 
for i=1:length(mongs) ; 
   EEG = pop_loadbv('.',mongs(i).name) ;  
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   EEG = denoise_grad(EEG) ; 
   EEG = pop_resample(EEG,512) ; 
   eegs{i} = EEG ; 
   EEG = denoise_bcg(EEG) ;   
   if i==1 ; merged = EEG  ;else merged = pop_mergeset(EEG,merged) ; end     
end
filt = merged ; filt.data = eegfiltfft(filt.data,filt.srate,1,90) ; 
ica = pop_runica(filt,'runica') ; 

figure,for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(ica.icawinv(:,i)),ica.chanlocs) ; title(i) ; end
for i=1:length(eegs) ; 
    eegs{i}.icaact = icaact(eegs{i}.data,ica.icaweights*ica.icasphere,0) ;
    eegs{i}.icawinv = ica.icawinv ; eegs{i}.icasphere = ica.icasphere ; eegs{i}.icaweights = ica.icaweights ; eegs{i}.icachansind = ica.icachansind ; 
end 
for i=1:length(eegs) ; pop_saveset(eegs{i},['ica_',num2str(i)]) ; end 

mnames = {
'MONG_01_RB_FIX_BOX.vhdr','MONG_01_RB_FT_BOX.vhdr','MONG_01_RB_GR_BOX.vhdr','MONG_01_RB_GR_SIN.vhdr',
'MONG_01_RB_REST.vhdr','MONG_01_RB_RUSSELL_1.vhdr','MONG_01_RB_RUSSELL_2.vhdr'  } ; 

trigs = {'S  1','S  2','S  3','S  4','S  5','S  6'} ;
names = {'fovea','periphery','upper','lower','right','left'} ;
for t =1:length(trigs)
ep = pop_epoch(eegs{6},{trigs{t}},[-1,3]) ; clear ersp ; 
for j=1:64 ; 
    [ersp(j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(j,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,'plotersp','off','plotitc','off',...
        'freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',0) ; 
        
end
%figure,for i=1:64 ; subplot(5,13,i) ;imagesc(squeeze(ersp(i,:,:)),[-5,5]) ; title(i) ;end
figure,imagesc(squeeze(mean(ersp([40,15],:,:),1)),[-5,5])
allersp2(t,:,:) = squeeze(mean(ersp([40,15],:,:),1)) ;
end


% process the FMRI
% get the triggers for each volume separately from the EEG files
geegs{1} = eegs{6} ; geegs{2} = eegs{7} ; 
clear triginds volOnsets lats
for i=1:length(geegs) ; 
    volOnsets{i} = find(strcmp('R128',{geegs{i}.urevent.type})) ;    
    latsi = cell2mat({geegs{i}.urevent.latency}) ; 
    lats{i} = latsi(volOnsets{i}) ; 
    for t=1:length(trigs)
        trigindst = find(strcmp(trigs{t},{geegs{i}.urevent.type})) ; % trigger indices of stimuli
        triglats = latsi(trigindst) ; % latencies for stimuli
        for lt=1:length(triglats) % for all trials of that stimulus type
           latdiffs = abs(triglats(lt) - latsi(volOnsets{i})) ; % the difference between that trial's latency and all TR latencies
           pnlatdiffs = (triglats(lt) - latsi(volOnsets{i})) ; % pos and negative latdiffs
           triginds(i,t,lt) = find(latdiffs==min(latdiffs)) ; % find the minimum difference
           offsets(i,t,lt) = pnlatdiffs(find(latdiffs==min(latdiffs))) ; % save the offsets from that TR to the actual trigger (positive difference means 
        end
    end
end

% get the FMRI data ; 
TR = 0.9 ; task = round(10/TR) ; 
hrf = spm_hrf(TR) ; 
clear epochs
for i=1:size(triginds,1) ; 
    cd C:\shared\raw\MONG_01_RB\gamma
    f1 = load_untouch_nii(['reg_',num2str(i),'.nii.gz']) ; 
    fimg = f1.img ; 
    for j=1:size(triginds,2) ; disp(j) ; 
        for k=1:size(triginds,3)
            epochs(j,:,:,:,(i*10-10)+k,:) = fimg(:,:,:,triginds(i,j,k):triginds(i,j,k)+task) ; 
        end
    end
end

%%% get the EEG single trials
trigs = {'S  1','S  2','S  3','S  4','S  5','S  6'} ;
names = {'fovea','periphery','upper','lower','right','left'} ;
clear ersp ; 
for sc=1:2
for t =1:length(trigs)
ep = pop_epoch(eegs{sc+5},{trigs{t}},[-1,3]) ; 
for j=1:64 ; 
    for k=1:size(ep.icaact,3)
        [ersp(j,t,(sc*10-10)+k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(j,:,k)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,'plotersp','off','plotitc','off',...
            'freqs',[1,120],'nfreqs',60,'winsize',128,'baseline',NaN) ; 
    end
end
end
end
base = squeeze(mean(ersp(:,:,:,:,times<0),5)) ; base = repmat(base,[1,1,1,1,200]) ; 
corrected = ersp-base ; 

for comps=1:64 ; disp(comps) ; 
mcorrect = squeeze(mean(corrected(comps,:,:,:,:),1)) ; 
diffepochs = squeeze(mean(epochs(:,:,:,:,:,7:end),6)-epochs(:,:,:,:,:,1)) ; 
voxcorrs = zeros(6,size(diffepochs,2),size(diffepochs,3),size(diffepochs,4),60) ; 
for st= 1:size(diffepochs,1) ; % disp(st) ; 
for i=1:60 ;
    voxcorrs(st,:,:,:,i) = voxcorr(squeeze(diffepochs(st,:,:,:,:)),squeeze(mean(mcorrect(st,:,i,times>0 & times<2),4))) ; 
 % figure,
 %   imagesc(squeeze(mean(mcorrect(1,:,i,:),1)))
end
end
mvox = squeeze(mean(voxcorrs,1)) ; 
save_nii(make_nii(mvox),[num2str(comps),'.nii.gz']) ;

end
% not much on the spontaneous single trials...

cd c:\shared\raw\MONG_01_RB\gamma\ica2
a = load('melodic_mix') ; 
goods = [19,20,21,28,31,34,37,45,46,49,51,52,55,64,66,68,70,74,75] ;
all = 1:size(a,2) ; 
icount = 1 ; clear pfilts ; 
for i=1:5:100 
    pfilts(icount,:,:) = eegfiltfft(eegs{6}.icaact,eegs{6}.srate,i,i+5) ; 
    icount = icount + 1 ; 
end
incr = ceil(size(pfilts,3)/size(a,1)) ; 
pfilts = pfilts.^2 ; 

for i=1:size(pfilts,1) ; disp(i) ; 
    for j=1:size(pfilts,2) ;
        kcount = 1 ;
        for k=1:incr:size(pfilts,3)-incr
            resfilts(i,j,kcount) = squeeze(mean(pfilts(i,j,k:k+incr),3)) ; 
            kcount = kcount + 1 ; 
        end
    end
end

hrf = spm_hrf(0.9) ; 
for i=1:size(resfilts,1)
    for j=1:size(resfilts,2)
        conved(i,j,:) = conv(squeeze(resfilts(i,j,:)),hrf,'same') ; 
    end
end

for i=1:size(conved,1)
    for j=1:size(conved,2)
        for k=1:size(a,2)
            corrs(i,j,k) = corr2(squeeze(a(1:963,k)),squeeze(conved(i,j,:))) ;             
        end
    end
end

for i=1:86 ; b(:,i) = detrend(a(:,i)) ; end


%{
for i=1:6 ; 
    subplot(2,3,i) ; imagesc(times,freqs,squeeze(allersp(i,:,:)),[-5,5]) ;
    title(names{i}) ; vline([0,2],'k') ; 
end
suptitle('scan 2') ; 
subplot(1,2,1) ; 
plot(squeeze(mean(mean(allersp2(:,freqs>10&freqs<25,times>0&times<2),2),3)),squeeze(mean(mean(allersp2(:,freqs>60&freqs<80,times>0&times<2),2),3)),'o') 
title(['scan1 corr2 = ',num2str(corr2(squeeze(mean(mean(allersp2(:,freqs>10&freqs<25,times>0&times<2),2),3)),squeeze(mean(mean(allersp2(:,freqs>60&freqs<80,times>0&times<2),2),3))))])
xlabel('10-25Hz') ; ylabel('60-80Hz') ; 
subplot(1,2,2) ; 
plot(squeeze(mean(mean(allersp(:,freqs>10&freqs<25,times>0&times<2),2),3)),squeeze(mean(mean(allersp(:,freqs>60&freqs<80,times>0&times<2),2),3)),'o') 
title(['scan2 corr2 = ',num2str(corr2(squeeze(mean(mean(allersp(:,freqs>10&freqs<25,times>0&times<2),2),3)),squeeze(mean(mean(allersp(:,freqs>60&freqs<80,times>0&times<2),2),3))))])
xlabel('10-25Hz') ; ylabel('60-80Hz') ; 
%}

cd c:\shared\raw\MONG_01_RB\ica_4
%cd c:\shared\raw\MONG_01_RB\ica_1
a = load('melodic_mix') ; 
b = eegfiltfft(a',1/0.9,0.01,.5) ;
b = b(:,20:end-20) ; 
eeg = eegs{4} ;
actdata = eegs{4}.icaact ; 
clear corrmats ; clear ersp  
for co=1:64 ; 
[ersp(co,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(actdata(co,:)),eeg.pnts,[eeg.xmin,eeg.xmax],eeg.srate,0,...
    'baseline',NaN,'freqs',[1,120],'nfreqs',120,'timesout',size(a,1),'plotersp','off','plotitc','off','winsize',512) ; 
end
clear convs
for i=1:64 ;
    for j=1:size(ersp,2)
        convs(i,j,:) = conv(squeeze(ersp(i,j,:)),hrf) ; 
    end
end
convs = convs(:,:,20:size(b,2)+20-1) ; 
clear corrmat ; 
for i=1:size(convs,1) 
    for j=1:size(convs,2)
        for k=1:size(b,1)
            corrmat(i,j,k) = corr2(b(k,:)',squeeze(convs(i,j,:))) ; 
        end
    end
end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(corrmat(i,:,:)),[-1,1]) ; title(i) ; end

%%%% get the triggers
onsets = find(strcmp('S  1',{eegs{3}.urevent.type})) ; 
lats = {eegs{3}.urevent.latency} ; onsetlats = round((cell2mat({lats{onsets}})./eegs{3}.srate)/.9-round(20*.9)) ; 
img = zeros(1,size(convs,3)) ; for i=1:length(onsetlats)-1 ; img(onsetlats(i):onsetlats(i)+20/.9) = 1 ; end ;  img = double(~img) ; img(img==0) = 0.5 ; 


ec = 10 ; fr = 21 ; fc = 13 ; 
%figure, imagesc(img,[0,1]) ; colormap gray ; hold on  ; shiftamt = .5 ;
figure,plot(1.5-(mat2gray(squeeze(mean(convs(ec,fr,:),2)))),'r','LineWidth',3) ; 
hold on ; plot(1.5-(mat2gray(b(fc,:))),'o','LineWidth',1) ;  plot(1.5-(mat2gray(b(fc,:))),'LineWidth',2) ; legend({'eeg','fmri'}) ; 
xlabel('time(s)') ; set(gca,'XTick',1:25:length(times)-39,'XTickLabel',round(times(1:25:end-39))) ; vline(onsetlats,'k') ; set(gca,'YTick',[]) ; ylabel('arb. units') ;


plot((squeeze(mean(convs(ec,fr,:),2))),(b(fc,:))','.') ; title(num2str(corr2((squeeze(mean(convs(ec,fr,:),2))),(b(fc,:))'))) ; ylabel('bold') ; xlabel('eeg') ; 




% for spontaneous fluctuations

ec = 43 ; fr = 6 ; fc = 13 ; 
%figure, imagesc(img,[0,1]) ; colormap gray ; hold on  ; shiftamt = .5 ;
figure,plot((mat2gray(squeeze(mean(convs(ec,fr,:),2)))),'r','LineWidth',3) ; 
hold on ; plot((mat2gray(b(fc,:))),'o','LineWidth',1) ;  plot((mat2gray(b(fc,:))),'LineWidth',2) ; legend({'eeg','fmri'}) ; 
xlabel('time(s)') ; set(gca,'XTick',1:25:length(times)-39,'XTickLabel',round(times(1:25:end-39))) ; vline(onsetlats,'k') ; set(gca,'YTick',[]) ; ylabel('arb. units') ;xlim([1,321]) ;
