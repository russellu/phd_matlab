
cd('C:\shared\badger\jeremie') ; ls  ; clear all  ; close all
sounds=dir('*retino_*vhdr') ;
for i=2%:max(size(sounds)) ; 
   EEG = pop_loadbv('.',sounds(i).name) ; 
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
   EEG = denoise_grad3(EEG) ; 
   EEG = pop_resample(EEG,256) ;
   %EEG = denoise_bcg2(EEG) ; 
   %eegs{i} = EEG ; 
   %if i==1 ; merged = EEG ; else merged = pop_mergeset(EEG,merged) ; end   
end


filteeg = eegfiltfft(EEG.data,EEG.srate,5,50) ; 
meanpost = filteeg(32,:) ; 
for elec = 1:64 
clear diffepochs
wsize = round(EEG.srate*1) ; wsize2 = round(EEG.srate*.5) ; icount = 1 ;
clear maxinds mepochs; 
for i=wsize:wsize:size(EEG.data,2)-wsize*2
    maxinds(icount) = find(meanpost(i:i+wsize)==max(meanpost(i:i+wsize)))+i ; 
    mepochs(icount,:) = meanpost(maxinds(icount)-round(wsize2):maxinds(icount)+round(wsize2)) ; 
    icount = icount + 1 ; 
end
meanepoch = mean(mepochs,1) ; 

icount = 1 ; clear maxinds mepochs corrmepochs ; 
for i=wsize:wsize2:size(EEG.data,2)-wsize*2
    maxinds(icount) = find(meanpost(i:i+wsize2)==max(meanpost(i:i+wsize2)))+i ; 
    mepochs(icount,:) = EEG.data(elec,maxinds(icount)-round(wsize/4):maxinds(icount)+round(wsize/2)) ; 
    corrmepochs(icount,:) = meanpost(maxinds(icount)-round(wsize2):maxinds(icount)+round(wsize2)) ; 
    icount = icount + 1 ; 
end
mepochs((diff(maxinds)==0),:) = [] ; corrmepochs((diff(maxinds)==0),:) = [] ; maxinds((diff(maxinds)==0)) = []  ;
thresh = EEG.srate*.5 ; % 1/3 a second
%for s=1:10 
bads = zeros(size(maxinds)) ; 
for i=2:length(maxinds)
    tdiff = maxinds(i)-maxinds(i-1) ; 
    if tdiff < thresh
        correp1 = corr2(corrmepochs(i-1,:),meanepoch) ; 
        correp2 = corr2(corrmepochs(i,:),meanepoch) ; 
        if correp1 < correp2
            bads(i-1) = 1 ; 
        else 
            bads(i) = 1 ; 
        end
    end
end
maxinds(bads==1) = [] ; 
mepochs(bads==1,:) = [] ; corrmepochs(bads==1,:) = [] ;
%end
%meancorrs = corr(corrmepochs',meanepoch') ; 
%badzs = find(meancorrs<.5) ; 
%maxinds(badzs) = [] ; 
%mepochs(badzs,:) = [] ; 

cmat = corr(corrmepochs') ; 
sumc = mean(cmat) ; badc = find(sumc<0.75) ; 
%plot(filteeg(32,:)) ;  vline(maxinds) ; vline(maxinds(badc),'k') ;
% shift the poorly correlating epochs, try to find a better match.

shiftsz = 128 ; 
for i=1:length(badc)
    % shift and correlate
    indi = badc(i) ; 
    jcount = 1 ; clear corrs
    for j=-shiftsz:shiftsz
        shiftedepoch = meanpost(maxinds(indi)+j-round(wsize2):maxinds(indi)+j+round(wsize2)) ; 
        corrs(jcount) = corr2(shiftedepoch,meanepoch) ; 
        jcount = jcount + 1 ;  
    end
    maxc = find(corrs==max(corrs)) ; 
    newbads(i) = maxinds(indi) + maxc-shiftsz ; 
    maxinds(indi) = maxinds(indi) + maxc-shiftsz ; 
    mepochs(indi,:) = EEG.data(elec,maxinds(indi)-round(wsize/4):maxinds(indi)+round(wsize/2)) ; 
    corrmepochs(indi,:) = meanpost(maxinds(indi)-round(wsize2):maxinds(indi)+round(wsize2)) ; 
end
cmat = corr(corrmepochs') ; 
sumc = mean(cmat) ; badc = find(sumc<0.75) ; 
shiftsz = 128 ; 
for i=1:length(badc)
    % shift and correlate
    indi = badc(i) ; 
    jcount = 1 ; clear corrs
    for j=-shiftsz:shiftsz
        if maxinds(indi)+j-round(wsize2) >0
        shiftedepoch = meanpost(maxinds(indi)+j-round(wsize2):maxinds(indi)+j+round(wsize2)) ; 
        corrs(jcount) = corr2(shiftedepoch,meanepoch) ; 
        end
        jcount = jcount + 1 ;  
    end
    maxc = find(corrs==max(corrs)) ; 
    newbads(i) = maxinds(indi) + maxc-shiftsz ; 
    maxinds(indi) = maxinds(indi) + maxc-shiftsz ; 
    mepochs(indi,:) = EEG.data(elec,maxinds(indi)-round(wsize/4):maxinds(indi)+round(wsize/2)) ; 
    corrmepochs(indi,:) = meanpost(maxinds(indi)-round(wsize2):maxinds(indi)+round(wsize2)) ; 
end
plot(filteeg(32,:)) ;  vline(maxinds) ;
avamt = 20 ;
for ep=1:size(mepochs,1) ; 
    cmat = corr(mepochs') ; 
    [~,si] = sort(cmat(ep,:),'descend') ; 
    diffepochs(elec,ep,:) = mepochs(si(1),:) - mean(mepochs(si(2:avamt),:)) ; 
end
end


EEG2 = EEG ; 
for i=1:length(maxinds)
   EEG2.data(:,maxinds(i)-round(wsize/4):maxinds(i)+round(wsize/2)) = squeeze(diffepochs(:,i,:)) ; 
end



