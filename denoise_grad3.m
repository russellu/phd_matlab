function EEGret = denoise_grad3(EEG)
% function EEGret = denoise_grad3(EEG)

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

for i=1:size(gradords,2) ; disp(['processing channel ',num2str(i)])  ;
    cmat = corr(squeeze(gradords(:,i,:))') ;
    for j=1:size(cmat,1)
        [~,cinds] = sort(double(cmat(:,j)),'descend') ; 
        meangrad = squeeze(mean(gradords2(cinds(2:50),i,:),1)) ; 
        EEG2.data(i,gradlats(j):gradlats(j)+diffgradlats(2)-1) = EEG.data(i,gradlats(j):gradlats(j)+diffgradlats(2)-1)-meangrad' ; 
    end
end

trigs = {EEG2.urevent.type} ; 
st = find(strcmp('S 98',trigs)) ; en = find(strcmp('S 99',trigs)) ; 
lats = cell2mat({EEG2.urevent.latency}) ; 
allt = find(strcmp('R128',trigs)) ; 
l1 = lats(allt(1)) ; l2 = lats(allt(end)) ; 
%EEGret = pop_select(EEG2,'point',[lats(st)+5000*2,lats(en)]) ; 
EEGret = pop_select(EEG2,'point',[l1,l2]) ; 
%EEGret = EEG2 ; 
end