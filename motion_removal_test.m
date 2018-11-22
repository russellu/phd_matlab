clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','karl','russell','tegan','sukhman','valerie'} ;
subns = [1,2,3,4,6,7,8,9] ; 
% EEG comps 
comps = {[5,12,14,16,21,19],...
    [6,14,16,17,25],...
    [6,10,11,13,16,28],...
    [6,9,23],...
    [12,15,18],...
    [8,15,20],...
    [7,13,21,31],...
    [6,15,21,22],...
    [8,15,18]} ;

for sub=8%:length(subns) ; 
cd(['c:/shared/badger_eeg2/',subs{subns(sub)}]) ; ls
mov = dir('bcgica*rest*set') ; 
EEG = pop_loadset(mov(1).name) ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
allica = load('allica') ; allica = allica.allica ; 
weights = allica{1} ; sphere = allica{2} ; EEG.icaact=weights*sphere*EEG.data ; 
ica = EEG ; ica.icaact = weights*sphere*EEG.data ; 
winv = pinv(weights*sphere) ;
rawts = mean(EEG.data([61,29],:),1) - mean(EEG.data([62,30],:),1) ;
freqs = 1:100 ; 
filts = zeros(size(ica.icaact,1),length(freqs),size(ica.icaact,2)) ; 
for j=1:length(freqs)
    filts(:,j,:) = eegfiltfft(ica.icaact,250,freqs(j)-1,freqs(j)+1) ; disp(j) ; 
end

filts = abs(filts) ;
resfilts = zeros(size(filts,1),size(filts,2),round(size(filts,3)/8)) ; 
for i=1:size(filts,1)
    resfilts(i,:,:) = imresize(squeeze(filts(i,:,:)),[size(filts,2),round(size(filts,3)/8)]) ; disp(i) ; 
end

meanres = (log(abs(squeeze(mean(resfilts(:,:,:),1))))) ;
meanres = meanres-min(meanres(:)) ; 
baserep = meanres - repmat(mean(meanres,2),[1,size(meanres,2)]) ; 
for i=1:size(baserep,1) ; baserep(i,zscore(baserep(i,:))<-2) = mean(baserep(i,:)) ; end
baserep = eegfiltfft(baserep,250/8,0.02,100) ; 
binrep = zeros(size(baserep)) ; negbinrep = zeros(size(baserep)) ; 
for i=1:100 
[z,mu,sig] = zscore(baserep(i,:)) ; 
bads = baserep(i,:) > mu+sig/2 ; 
goods = baserep(i,:) < mu  ; 
binrep(i,bads) = 1 ; negbinrep(i,goods) = 1 ; 
end
figure,subplot(1,2,1) ; imagesc(binrep) ; subplot(1,2,2) ; imagesc(baserep,[-2,2]) ; 
binrep = imdilate(binrep,strel(ones(3,9)))  ; 
binrep = 1-(negbinrep - binrep >=0) ;

newcomps = zeros(size(filts)) ; 
for i=1:64 ; disp(i) ; 
currentcomp = squeeze(filts(i,:,:)) ; 
[gx,gy] = meshgrid(1:size(binrep,2),1:size(binrep,1)) ; 
bads = find(binrep==1) ; goods = find(binrep==0) ; 
badx = gx(bads) ; bady = gy(bads) ; goodx = gx(goods) ; goody = gy(goods) ; 
goodv = currentcomp(goods) ; 
vq = griddata(goodx,goody,goodv,badx,bady,'cubic') ; 
replaced = currentcomp ; replaced(bads) = vq ; 
subplot(1,2,1),imagesc(baserep,[-2,2])
subplot(1,2,2),imagesc(replaced,[-2,2])
newcomps(i,:,:) = replaced ; 
end

imagesc(squeeze(mean(newcomps(comps{subns(sub)},:,:),1))) ; 

end




