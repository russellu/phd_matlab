function [newzdat,eeg4,ss] = fbcg(EEG)
% function [newzdat,eeg4,ss] = fbcg(EEG)

% remove BCG artifacts
EEG.data(32,:) = rand(1,size(EEG.data,2))*0.0005 ; 
[s,f] = spectopo(EEG.data(:,5000:end-5000),0,250,'plot','off') ; EEG.data = eegfiltfft(EEG.data,250,1,250) ; 
%EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
filtEEG = EEG ; %filtEEG.data = eegfiltfft(EEG.data,250,0,250) ;
filtnograd = EEG ; filtnograd.data = eegfiltfft(EEG.data,250,3,25) ; 
ss = spectopo(filtEEG.data,0,250,'plot','off') ; 
filtnograd.data(32,:) = 0 ; 
ica = pop_runica(filtnograd(:,1:5:end),'runica') ; 
ica = ica_applyweights(filtnograd,ica) ; 

clear wskews ; 
wsize = 250*20 ; 
for i=1:5 ; jcount = 1 ; 
    for j=1:wsize:size(ica.icaact,2)-wsize
        wskews(i,jcount) = skewness(ica.icaact(i,j:j+wsize)) ; 
        jcount = jcount + 1 ; 
    end
end
skews = median((wskews),2) ; 
maxind = find(abs(skews(1:5))==max(abs(skews(1:5)))) ; 
polarity = 0 ; if skews(maxind)<0 ; polarity = -1 ; else polarity = 1 ; end
template = ica.icaact(maxind,:)*polarity ; figure,plot(template) ; 


wsize = 300 ; wincr = 50 ; 
icount = 1 ; clear allepochs maxinds
for i=wsize:wincr:length(template)-wsize*2
    epoch = template(i:i+wsize) ; 
    maxi = find(epoch==max(epoch))+i ; 
    maxinds(icount,:) = maxi-50:maxi+250 ; 
    allepochs(icount,:) = template(maxinds(icount,:)) - mean(template(maxinds(icount,:))) ; 
    icount = icount + 1 ; 
end
clear sdiffs
for i=1:size(allepochs,1)
    sdiffs(i,:) = sum((repmat(allepochs(i,1:150),[size(allepochs,1),1])-allepochs(:,1:150)).^2,2) ;
end
[sv,si] = sort(sdiffs,2,'ascend') ; 
allhrs = zeros(64,size(maxinds,1),size(maxinds,2)) ; 
for i=1:size(EEG.data,1)
    for j=1:size(maxinds,1)
        allhrs(i,j,:) = filtEEG.data(i,maxinds(j,:)) ; 
    end
end
% avg the indices 
mchanepochs = zeros(size(allhrs,1),size(allhrs,2),size(allhrs,3)) ; 
for i=1:size(allhrs,1) ;
    for j=1:size(allhrs,2)
    	mchanepochs(i,j,:) = squeeze(mean(allhrs(i,si(j,3:50),:),2)) ; 
    end      
end
zdat = EEG.data; 
for i=1:size(mchanepochs,1)
    for j=1:size(mchanepochs,2)
        zdat(i,maxinds(j,:)) = filtEEG.data(i,maxinds(j,:)) - squeeze(mchanepochs(i,j,:))' ; 
    end
end


eeg4 = EEG ; eeg4.data = zdat ; 
newzdat = zdat ; 

end
