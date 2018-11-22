cd('c:/shared/badger/alex') ; ls  ; clear all  ; close all
sounds=dir('*retino_gamma*vhdr') ;
for i=1
   EEG = pop_loadbv('.',sounds(i).name) ; 
   EEG = pop_chanedit(EEG,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
   
end
%{
[EEGc,neweeg5k,mgradord] = denoise_grad2(EEG) ; 
EEG2 = EEG ; 
gradtriginds = find(strcmp('R128',{EEG.urevent.type})) ; 
lats = {EEG.urevent.latency} ; 
gradlats = cell2mat(lats(gradtriginds)) ; gradlats(length(gradlats)) = [] ; 
diffgrad = diff(gradlats) ;
trlen = diffgrad(1)  ;
gradords = zeros(length(gradlats),size(EEG.data,1),trlen) ; 
% isolate each TR
for gr=1:length(gradlats)
    gradords(gr,:,:) = EEG.data(:,gradlats(gr):gradlats(gr)+trlen-1) ;      
end

% pad that shit, son
padamt = 50 ; % pad amount on either side
padgrads = zeros(size(gradords,1)+padamt*2,size(gradords,2),size(gradords,3)) ; 
padgrads(1:padamt,:,:) = flipdim(gradords(1:padamt,:,:),1) ; padgrads(end-padamt+1:end,:,:) = flipdim(gradords(end-padamt+1:end,:,:),1) ; 
padgrads(padamt+1:end-padamt,:,:) = gradords ; 

for e=1:64 ; disp(e) ; 
g1 = squeeze(padgrads(:,e,:)) ; 
meddiff = g1-imfilter(g1,fspecial('gaussian',[51,1],11)) ; 
subgrads(:,e,:) = meddiff(padamt+1:end-padamt,:) ; 
end

for gr=1:size(gradords,1)
    EEG2.data(:,gradlats(gr):gradlats(gr)+trlen-1) = squeeze(subgrads(gr,:,:)) ;           
end

eeg5k = EEG2 ; 
EEG2 = pop_resample(EEG2,256) ; 
gradtriginds = find(strcmp('R128',{EEG2.urevent.type})) ; 
lats = {EEG2.urevent.latency} ; 
gradlats = cell2mat(lats(gradtriginds)) ; 
diffgrad = diff(gradlats) ;
trlen = diffgrad(1)  ; 
EEG2.data(:,1:gradlats(1)+trlen/2) = 0 ; 
EEG2.data(:,gradlats(length(gradlats))+trlen/2:end) = 0 ; 

[s,~] = spectopo(EEG2.data,0,EEG2.srate,'plot','off') ;
[s2,f] = spectopo(EEGc.data,0,EEG2.srate,'plot','off') ;
subplot(1,2,1) ; imagesc(s) ; subplot(1,2,2) ; imagesc(s-s2); 
figure
plot(neweeg5k.data(10,:)) ; hold on ; plot(eeg5k.data(10,:),'r') ;
%}

ulimcount = 1 ; clear uepochs
for ulim = 5:25:100 
EEGn = EEG ; 
EEGn.data = EEGn.data - eegfiltfft(EEGn.data,EEGn.srate,0,ulim) ; 
gradtriginds = find(strcmp('R128',{EEG.urevent.type})) ; 
lats = {EEG.urevent.latency} ; 
gradlats = cell2mat(lats(gradtriginds)) ; gradlats(length(gradlats)) = [] ; 
diffgrad = diff(gradlats) ;trlen = diffgrad(1) ;
gradords = zeros(length(gradlats),size(EEG.data,1),trlen) ;
gradords2 = zeros(length(gradlats),size(EEG.data,1),trlen) ; 
for gr=1:length(gradlats) ; gradords2(gr,:,:) = EEG.data(:,gradlats(gr):gradlats(gr)+trlen-1) ; end
% isolate each TR
for gr=1:length(gradlats) ; gradords(gr,:,:) = EEGn.data(:,gradlats(gr):gradlats(gr)+trlen-1) ; end
cmat = corr(squeeze(gradords(:,1,:))') ; corrcount = 1 ; 
for ncorrs = 10:15:50
epochs = zeros(size(gradords,1),ncorrs,size(gradords,3)) ; 
for i=1:size(cmat,1) ; 
    [~,si] = sort(cmat(i,:),'descend') ; 
    epochs(i,:,:) = squeeze(gradords2(si(1:ncorrs),1,:)) ; 
end
stdepochs = squeeze(std(epochs,0,2)) ;
allstd(corrcount,:) = (mean(stdepochs,1)) ; corrcount = corrcount + 1 ; disp(ncorrs) ;  
uepochs(ulimcount,corrcount,:) = squeeze(mean(stdepochs,1)) ; 
end
ulimcount = ulimcount + 1 ;
end



EEGn = EEG ; EEG2 = EEG ; 
EEGn.data = EEGn.data - eegfiltfft(EEGn.data,EEGn.srate,0,200) ; 
gradtriginds = find(strcmp('R128',{EEG.urevent.type})) ; 
lats = {EEG.urevent.latency} ; 
gradlats = cell2mat(lats(gradtriginds)) ; gradlats(length(gradlats)) = [] ; 
diffgrad = diff(gradlats) ; trlen = diffgrad(1)  ;
gradords = zeros(length(gradlats),size(EEG.data,1),trlen) ;
gradords2 = zeros(length(gradlats),size(EEG.data,1),trlen) ; 
for gr=1:length(gradlats) ; gradords2(gr,:,:) = EEG.data(:,gradlats(gr):gradlats(gr)+trlen-1) ; end
for gr=1:length(gradlats) ; gradords(gr,:,:) = EEGn.data(:,gradlats(gr):gradlats(gr)+trlen-1) ; end
diffgradlats = diff(gradlats) ;

for i=1:size(gradords,2) ; disp(i)  ;
    cmat = corr(squeeze(gradords(:,i,:))') ;
    cmats(i,:,:) = cmat ; 
    for j=1:size(cmat,1)
        [cv,cinds] = sort(double(cmat(:,j)),'descend') ; 
        meangrad = squeeze(mean(gradords2(cinds(2:50),i,:),1)) ; 
        EEG2.data(i,gradlats(j):gradlats(j)+diffgradlats(2)-1) = EEG.data(i,gradlats(j):gradlats(j)+diffgradlats(2)-1)-meangrad' ; 
    end
end
%{
EEGres2 = EEG2 ; 
EEGres = eegfiltfft(EEG2.data,5000,78,80) ; 
EEGres2.data = EEGres ; 
ica = pop_runica(EEGres2,'runica') ; 

EEG23 = pop_resample(EEG2,256) ; 

gradtriginds = find(strcmp('R128',{EEG23.urevent.type})) ; 
lats = {EEG23.urevent.latency} ; 
gradlats = cell2mat(lats(gradtriginds)) ; 
diffgrad = diff(gradlats) ;
trlen = diffgrad(1)  ; 
EEG23.data(:,1:gradlats(1)+trlen/2) = 0 ; 
EEG23.data(:,gradlats(length(gradlats))+trlen/2:end) = 0 ; 
[s,f] = spectopo(EEG23.data,0,EEG23.srate,'plot','off');

s2 = s ; 

EEGuncorr = EEG23 ;  
plot(mat2gray(EEG2.data(10,5000*20:end-5000*20))) ; hold on ; plot(mat2gray(EEG.data(10,5000*20:end-5000*20)),'r')

ica = pop_runica(EEG23,'runica') ; 

[s,f] = spectopo(EEG2.data,0,EEG2.srate) ; 

[s2,f2] = spectopo(EEG.data,0,EEG2.srate) ; 

gradtriginds = find(strcmp('R128',{EEG.urevent.type})) ; 
lats = {EEG.urevent.latency} ; 
gradlats = cell2mat(lats(gradtriginds)) ; gradlats(length(gradlats)) = [] ; 
diffgrad = diff(gradlats) ; trlen = diffgrad(1)  ;
gradords = zeros(length(gradlats),size(EEG.data,1),trlen) ;
gradords2 = zeros(length(gradlats),size(EEG.data,1),trlen) ; 
for gr=1:length(gradlats) ; gradords2(gr,:,:) = EEG2.data(:,gradlats(gr):gradlats(gr)+trlen-1) ; end
for gr=1:length(gradlats) ; gradords(gr,:,:) = EEGn.data(:,gradlats(gr):gradlats(gr)+trlen-1) ; end
diffgradlats = diff(gradlats) ;
%}

% run ica on the high pass filtered gradient-subtracted data
EEGn = EEG2 ; 
trigs = {EEGn.urevent.type} ; 
st = find(strcmp('S 98',trigs)) ; en = find(strcmp('S 99',trigs)) ; 
lats = cell2mat({EEGn.urevent.latency}) ; 
EEGn = pop_select(EEGn,'point',[lats(st)+5000*2,lats(en)]) ; EEGs = pop_select(EEG2,'point',[lats(st)+5000*2,lats(en)]) ; 
filtn = EEGn ; filtn.data = filtn.data - eegfiltfft(EEGn.data,EEGn.srate,0,80) ; 
icafilt = pop_runica(filtn,'runica') ; 

EEGn.icaact = icaact(EEGn.data,icafilt.icaweights*icafilt.icasphere,0) ; 
EEGn.icachansind = icafilt.icachansind ; 
EEGn.icasphere = icafilt.icasphere ; 
EEGn.icasplinefile = icafilt.icasplinefile ; 
EEGn.icaweights = icafilt.icaweights ; 
EEGn.icawinv = icafilt.icawinv ;   


EEGns = pop_subcomp(EEGn,1) ; 
EEGsub = pop_resample(EEGns,256) ; EEGres = pop_resample(EEGs,256) ; 
[sp,fs] = spectopo(EEGsub.data,0,EEGsub.srate,'plot','off') ; 
[sp2,fs2] = spectopo(EEGres.data,0,EEGres.srate,'plot','off') ; 


EEG4 = pop_resample(EEGs,256) ; EEG4.data = EEGresn.data ; 

EEG4 = pop_runica(EEG4,'runica') ; 
[s,f] = spectopo(EEG4.icaact,0,EEG4.srate,'plot','off') ; 



res1 = pop_resample(EEG2,256) ; 
res2 = pop_resample(EEGs,256) ;
[s1,f1] = spectopo(res1.data,0,res1.srate,'plot','off') ; 
[s2,f2] = spectopo(res2.data,0,res2.srate,'plot','off') ; 

res2 = denoise_bcg(res2) ; 
res2filt = res2 ; res2filt.data = eegfiltfft(res2.data,res2.srate,1,128) ; 
ica = pop_runica(res2filt,'runica') ; 



ep = pop_epoch(ica,{'S  1','S  2'},[-2,7]) ; 
comps=1:64 ; clear ersp ; 
ep.icaact(:,:,[16,17]) = [] ; 
for c=1:length(comps) ; 
        [ersp(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(comps(c),:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'timesout',200) ;     
end
for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-7,7]) ; end











