clear all ; close all ; 
cd C:\shared\retino_eeg
discr = dir('russ_discrete*vhdr') ; 
for i=1:length(discr)
   EEG = pop_loadbv('.',discr(i).name) ; 
   EEG = pop_resample(EEG,256) ; 
   if i==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG) ; end
end
merged = pop_chanedit(merged,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
temp = merged.data(1:32,:) ; merged.data(1:32,:) = merged.data(33:64,:) ; merged.data(33:64,:) = temp ; 
mergefilt = merged ; mergefilt.data = eegfiltfft(merged.data,merged.srate,1,128) ; 
ica = pop_runica(mergefilt,'runica') ; ica.icaact = ica.icaweights*ica.icasphere*ica.data ; 
ica.data = ica.icaact ; 
allep = pop_epoch(ica,{'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 22','S 23','S 24','S 25'},[-.25,2.25]) ;
clear allersp ; 
for i=1:64 ;
        [allersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
end

for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(allersp(i,:,:)),[-4,4]) ; title(i) ; axis xy ; end



goodcs = [11,12,20,21] ; bads = zeros(1,64) ; bads(goodcs) = 1 ; bads = find(bads==0) ; 
winv = pinv(ica.icaweights*ica.icasphere) ; 
acts = ica.icaact ; acts(bads,:) = 0 ; denoised = winv*acts ; 
ica.data = denoised ; 

stims = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 22','S 23','S 24','S 25'} ; 
titles = {'full','left','right','top','bot','360','315','270','225','180','135','90','45','bot_right','bot_left','top_left','top_right'} ; 
freqs = 1:4:120 ; 
alldats = zeros(17,64,length(freqs),768,20) ; 
for i=1:length(stims) 
    stimep = pop_epoch(ica,{stims{i}},[-.5,2.5]) ; 
    reseps = reshape(stimep.data,[64,768*20]) ; 
    for f=1:length(freqs)
        filt = eegfiltfft(reseps,ica.srate,freqs(f)-1,freqs(f)+1) ; 
        alldats(i,:,f,:,:) = reshape(filt,size(stimep.data)) ; 
    end
end
alldats = abs(alldats) ; 
mdats = squeeze(mean(alldats(:,:,:,stimep.times<2000 & stimep.times>0,:),4) - mean(alldats(:,:,:,stimep.times<0 & stimep.times>-250,:),4)) ; 
figure,for i=1:30 ; subplot(5,6,i) ; topoplot(squeeze(mean(mdats(5,:,i,:),4)),ica.chanlocs,'maplimits',[-.15,.15]) ; end
for i=1:17 ; subplot(4,5,i) ; topoplot(squeeze(mean(mean(mdats(i,:,freqs>8 & freqs<26,:),3),4)),ica.chanlocs,'maplimits',[-.5,.5]) ; title(titles{i}) ; end

mgamma = squeeze(mean(mean(mdats(:,:,freqs>8& freqs<25,:),3),4)) ; cbar = [-.1,.1] ; 
subplot(2,3,1) ; topoplot(mgamma(1,:),ica.chanlocs,'maplimits',cbar) ; title('full field') ; colorbar ;
subplot(2,3,2) ; topoplot(sum(mgamma(2:3,:),1),ica.chanlocs,'maplimits',cbar) ; title('right + left') ; colorbar ;
subplot(2,3,3) ; topoplot(sum(mgamma(4:5,:),1),ica.chanlocs,'maplimits',cbar) ; title('top + bottom') ; colorbar ;
subplot(2,3,4) ; topoplot(sum(mgamma(14:end,:),1),ica.chanlocs,'maplimits',cbar) ; title('quad1+quad2+quad3+quad4') ; colorbar ;
subplot(2,3,5) ; topoplot(sum(mgamma(6:13,:),1),ica.chanlocs,'maplimits',cbar) ; title('oct1+oct2+oct3+oct4+oct5+oct6+oct7+oct8') ; colorbar ;

malpha = squeeze(mean(mean(mdats(:,:,5,:),3),4)) ; cbar = [-.1,.1] ; 
for i=1:17 ; subplot(4,5,i) ; topoplot(mgamma(i,:),ica.chanlocs) ; title(titles{i}) ; end
mgamma = squeeze(mean(mean(mdats(:,:,12:16,:),3),4)) ; cbar = [-.1,.1] ; 
figure,for i=1:17 ; subplot(4,5,i) ; topoplot(mgamma(i,:),ica.chanlocs) ; title(titles{i}) ; end

for i=1:17
    [f,topos(i,:,:)] = topoplot(mgamma(i,:),ica.chanlocs) ; 
    [f,topos2(i,:,:)] = topoplot(-malpha(i,:),ica.chanlocs) ; 
end

clear rgbtopos
rgbtopos(:,:,1) = uint8(mat2gray(squeeze(topos2(17,:,:)))*255) ; 
%rgbtopos(:,:,2) = uint8(mat2gray(squeeze(topos(16,:,:)))*255) ; 
rgbtopos(:,:,3) = uint8(mat2gray(squeeze(topos2(15,:,:)))*255) ; 
imshow(rgbtopos) ; axis xy ; 


%{
clear ersp ; 
for i=1:length(stims)
    epi = pop_epoch(ica,{stims{i}},[-.5,2.5]) ; disp(i) ; 
    for j=1:64
        for k=1:20 ;
            [ersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.data(j,:,k)),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off') ; 
        end
    end
end
bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,200]) ; 
stdbersp = squeeze(mean(std(mean(bersp,4),0,5),2)) ; 
[sv,si] = sort(stdbersp,2,'ascend') ; 
clear mersp ; 
for i=1:17 ; 
    mersp(i,:) = squeeze(mean(mean(mean(bersp(i,:,si(i,5:end),freqs>8 & freqs<25,times>0 & times<2),3),4),5)) ; 
end
subplot(1,2,1) ; topoplot(mersp(1,:),ica.chanlocs,'maplimits',[-8,8]) ; 
subplot(1,2,2) ; topoplot(mersp(17,:)+mersp(16,:)+mersp(15,:)+mersp(14,:),ica.chanlocs,'maplimits',[-8,8]) ; 
%}







