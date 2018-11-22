%cd('c:/Vision/Raw Files/Russell_2015-10-15/') ; ls  ; clear all  ; close all
cd('c:/Vision/Raw Files/felix/') ; ls  ; clear all  ; close all

% get the raw EEG data
sounds=dir('*vhdr') ;
for i=4%:max(size(sounds)) ; 
   EEG = pop_loadbv('.',sounds(i).name) ; 
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
  % gradcmat = corr(EEG.data') ;
    EEG = denoise_grad2(EEG) ; % remove gradient artifact using sliding window AAS
    bcgcmat = corr(EEG.data') ; 
   %eegs{i} = EEG ; 
   %if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end   
end







% correlate all channels with each other, and average the top 5 posterior
% channels to form a heartbeat template independent of the BCG
corrs = corr(EEG.data') ; 
[sv,si] = sort(corrs(48,:),'descend') ;
sumbcg = (mean(EEG.data(si(1:5),:),1)) ;

largewind = round(EEG.srate*1.5) ; % window size =>  1/2 second
smallwind = round(EEG.srate/3) ; 
% find the initial min peaks, large window size (1.5s)
icount = 1 ; clear mininds slowepochs
for w=1:largewind:length(sumbcg)-largewind    
    minind = find(sumbcg(w:w+largewind)==min(sumbcg(w:w+largewind)),1) ;
    if sumbcg(minind+w) ~= 0 
        mininds(icount) = minind+w ;  
        slowepochs(icount,:) = sumbcg(mininds(icount)-smallwind:mininds(icount)+smallwind) ; 
        icount = icount + 1 ; 
    end    
end
mslows = squeeze(mean(slowepochs)) ; 
% find the initial min peaks, small window size (0.5s)
icount = 1 ; clear mininds shortepochs
for w=1:smallwind:length(sumbcg)-smallwind    
    minind = find(sumbcg(w:w+smallwind)==min(sumbcg(w:w+smallwind)),1) ;
    if sumbcg(minind+w) ~= 0 
        mininds(icount) = minind+w ;  
        shortepochs(icount,:) = smooth(sumbcg(mininds(icount)-smallwind:mininds(icount)+smallwind),15) ; 
        icount = icount + 1 ; 
    end    
end
correpochs = corr(shortepochs',mslows') ; 
icount = 1 ; clear goods ; 
for i=1:2:length(correpochs)-2
    maxi = find(correpochs(i:i+2)==max(correpochs(i:i+2)),1) ; 
    goods(icount) = maxi+i-1 ; 
    icount = icount + 1 ; 
end
goods0 = goods ; 
tvals = mininds(goods) ; 
difftvals = diff(tvals) ; 
repeats = find(difftvals==0) ; goods(repeats+1) = [] ;
difftvals2 = diff(mininds(goods)) ;
goodcorrs = correpochs(goods) ; 
% set a realistic limit: less than half a second,
mindiff = round(EEG.srate/1.6) ; 
icount = 1 ; clear bads
for i=1:length(difftvals2)
    if difftvals2(i) < mindiff
        corrv1 = goodcorrs(i) ;  corrv2 = goodcorrs(i+1) ; 
        if corrv1 < corrv2
            bads(icount) = i ;
        else
            bads(icount) = i+1 ; 
        end
        icount = icount + 1 ; 
    end    
end
goods2 = goods ; 
goods2(bads) = [] ;
%plot(EEG.data(32,:)) ; vline(mininds(goods0),'g') ;vline(mininds(goods)) ; vline(mininds(goods2),'k') ; 

% get the epochs, and the epoch indices
clear gepochs gepochinds
for i=1:length(goods2)
    gepochs(i,:,:) = EEG.data(:,mininds(goods2(i))-round(EEG.srate/2):mininds(goods2(i))+round(EEG.srate/2)) ;  
    gepochinds(i,:) = mininds(goods2(i))-round(EEG.srate/2):mininds(goods2(i))+round(EEG.srate/2) ;  
end

fwind = 20; % width of the BCG window
% frequency space subtraction for each heartbeat...based on the fwind
% heartbeats from that index onwards (or behind, if near end)
edat2 = EEG.data ; 
for i=1:size(gepochs,1) ; 
    if i+fwind > size(gepochs,1)
        meani = squeeze(mean(gepochs(i-fwind:i,:,:),1)) ; 
        cfft = fft(meani,[],2) ;
        gepochi = fft(squeeze(gepochs(i,:,:)),[],2) ; 
        gepochsub = real(ifft(gepochi-cfft,[],2)) ; 
        edat2(:,gepochinds(i,:)) = gepochsub ; 
    else     
        meani = squeeze(mean(gepochs(i:i+fwind,:,:),1)) ; 
        cfft = fft(meani,[],2) ;
        gepochi = fft(squeeze(gepochs(i,:,:)),[],2) ; 
        gepochsub = real(ifft(gepochi-cfft,[],2)) ; 
        edat2(:,gepochinds(i,:)) = gepochsub ;        
    end    
end

EEG2 = EEG ; 
EEG2.data = edat2 ; 
EEG2.data = eegfiltfft(EEG2.data,EEG2.srate,1,128) ; 
ica = pop_runica(EEG2,'runica') ; 
tp(ica)


[s,f] = spectopo(ica.icaact,0,ica.srate,'plot','off') ; 

% a second pass?
%clear gepochs2 gepochinds2
%for i=1:length(goods2)
%    gepochs2(i,:,:) = edat2(:,mininds(goods2(i))-round(EEG.srate/2):mininds(goods2(i))+round(EEG.srate/2)) ;  
%    gepochinds2(i,:) = mininds(goods2(i))-round(EEG.srate/2):mininds(goods2(i))+round(EEG.srate/2) ;  
%end

