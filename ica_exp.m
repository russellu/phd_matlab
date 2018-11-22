cd('c:/shared/denoising_MR') ; 
clear all ; close all ; 
%{
EEGdenoise = pop_loadbv('.','Russell_test_2015-08-05_EEGoutsideMRIroom.vhdr') ; 
EEGnoise = pop_loadbv('.','Russell_test_2015-08-05_EEGnoEPI.vhdr') ; 
EEGdenoise = pop_chanedit(EEGdenoise,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
EEGnoise = pop_chanedit(EEGnoise,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

n1 = pop_resample(EEGnoise,500) ; 
dn1 = denoise_bcg(n1) ; 
n2 = pop_resample(EEGdenoise,500) ;
merged = pop_mergeset(dn1,n2) ; 

%}
cd C:\shared\sim_1\first_simultaneous_recordings
merged = pop_loadset('bothmerged.set') ; 

wv = merged.data(46,3000:100000) ; 
[c,l] = wavedec(wv,3,'db1');
[A,D] = dwt(wv,'sym4');

dec = mdwtdec('r',wv,2,'db2');
[XD,decDEN,THRESH] = mswden('den',dec,'sqtwolog','sln');


merged = denoise_bcg(merged) ;
merged.data = eegfiltfft(merged.data,merged.srate,1,100) ; 
%merged = pop_select(merged,'nochannel',32) ; 
merged = pop_runica(merged,'runica') ;
for i=1:64 ; subplot(5,13,i) ; topoplot(squeeze(merged.icawinv(:,i)),merged.chanlocs) ; title(i) ; end

%eps = pop_epoch(merged,{'S  1'},[-1,3]) ; 
eps = pop_epoch(merged,{'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 14','S 15','S 16','S 17','S 18','S 19','S 20',},[-1,3]) ; 
v = squeeze(eps.icaact(31,:,:)) ; zts = zscore(v,0,1) ; maxzts = zscore(max(abs(zts),[],1)) ; 
goods = find(maxzts<2) ; 


for i=1:size(eps.icaact,1)
    [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(eps.icaact(i,:,goods(1:150))),...
        eps.pnts,[eps.xmin,eps.xmax],eps.srate,0,...
        'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'baseline',0,...
        'winsize',150) ; 

end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-3,3]) ; title(i) ; end


imagesc(squeeze(mean(ersp([14,20,34],:,:))),[-3,3])

