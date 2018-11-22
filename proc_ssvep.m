cd('c:/shared/ssvep') ; ls 
clear all ; close all ;

% inside the MR scanner files:
esets = dir('Test_Russell*vhdr') ;
for eset=1:length(esets)
EEG = pop_loadbv('.',esets(eset).name) ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
EEG = denoise_grad(EEG) ; 
EEG = pop_resample(EEG,500) ; 
%EEG = denoise_bcg(EEG) ; 

allsets{eset} = EEG ; 
end


osets = dir('out_ssvep_0*vhdr') ;
for oset=1:length(osets)
EEG = pop_loadbv('.',osets(oset).name) ;
EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
EEG = pop_resample(EEG,500) ; 

allosets{oset} = EEG ; 
end

merged = pop_mergeset(allosets{1},allsets{1}) ; 
merged = pop_mergeset(merged,allosets{2}) ; 
merged = pop_mergeset(merged,allosets{3}) ; 
merged = pop_mergeset(merged,allsets{2}) ; 
merged = pop_mergeset(merged,allosets{4}) ; 
merged = pop_mergeset(merged,allosets{5}) ; 
merged = pop_mergeset(merged,allsets{3}) ; 
merged = pop_mergeset(merged,pop_mergeset(allosets{6},allosets{7})) ; 

merged2 = pop_resample(merged,300) ; merged2.data = eegfiltfft(merged2.data,merged2.srate,1,150) ;
merged2 = pop_runica(merged2,'runica') ; 
for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(merged2.icawinv(:,i)),merged2.chanlocs) ; title(i) ; end

eps = pop_epoch(merged2,{'S  4'},[-.5,2.5]) ; 
plot(squeeze(eps.icaact(11,:,61:90)),'b') ; hold on ; plot(squeeze(eps.icaact(11,:,30:60)),'r') ;

insidemerged = pop_mergeset(pop_mergeset(allsets{1},allsets{2}),allsets{3}) ; 
imerge = pop_resample(insidemerged,300) ; imerge.data = eegfiltfft(imerge.data,imerge.srate,1,150) ; 
imerge = pop_runica(imerge,'runica') ; 
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(imerge.icawinv(:,i)),merged2.chanlocs) ; title(i) ; end

trigs = {'S  1','S  2','S  3','S  4'} ; clear ersp
for t=1:length(trigs)
    ieps = pop_epoch(merged2,{trigs{t}},[-.5,2.5]) ; 
    for s=1:size(ieps.icaact,1) ;
        [ersp(t,s,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ieps.icaact(s,:,31:60)),ieps.pnts,[ieps.xmin,ieps.xmax],ieps.srate,0,...
            'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',50) ; 
    end
end

channels = [11]  ;
%single trials
trigs = {'S  1','S  2','S  3','S  4'} ; clear tersp
for t=1:length(trigs)
    ieps = pop_epoch(merged2,{trigs{t}},[-.5,2.5]) ; 
    for s=1:length(channels)
        for tsk=1:size(ieps.icaact,3) ; 
            [tersp(t,s,tsk,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ieps.icaact(channels(s),:,tsk)),ieps.pnts,[ieps.xmin,ieps.xmax],ieps.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',50,'baseline',NaN) ; 
        end
    end
end











%{
for i=1:length(allsets)
   scanset_i = allsets{i} ; 
   scanset_i.icaact = icaact(scanset_i.data,merged2.icaweights*merged2.icasphere,0) ; 
   if i==1
      mergescans = scanset_i ; 
   else mergescans = pop_mergeset(mergescans,scanset_i) ; 
   end
end

for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(merged2.icawinv(:,i)),merged2.chanlocs) ; end
mergeps = pop_epoch(mergescans,{'S  1'},[-.5,2.5]) ; 
%}



%%% correlate icawinv with distance brain
cd C:\shared\raw\MONG_01_RB ; 
dists = load_untouch_nii('res_alldists.nii.gz') ; distim = dists.img ; 
winv = merged2.icawinv ; winv(32,:) = [] ; distim(:,:,:,32) = [] ; 
components = zeros(size(distim,1),size(distim,2),size(distim,3),64) ; 
for i=1:64 ; disp(i) ; 
   components(:,:,:,i) = voxcorr(distim,squeeze(winv(:,i))) ; 
end

[n1,x1] = hist(diff(allsets{1}.data(5,1:100000)),100) ; 
[n2,x2] = hist(diff(allosets{1}.data(5,1:100000)),100) ; 

plot(n1,'LineWidth',3) ; hold on; plot(n2,'r','LineWidth',3); xlim([x2(1),100]) ; legend({'inside bore','outside bore'}) ; title('timeseries DERIVATIVE histogram') ; 

s1 = allsets{1}.data(10,20000:100000) ;
s2 = allosets{1}.data(10,20000:100000) ;

j = wiener2(s1,[3,1],s2) ; 
plot(j) ; hold on ; plot(s1,'r') ;plot(s1-j,'k') ; 

plot(s1-j) ; 




