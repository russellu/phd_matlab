clear all ; close all 
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','olga','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 

for sub=1:length(subs)
    cd(['c:/shared/allres/',subs{sub}]) ; 
    ls
    chandata = pop_loadset('newchandata.set') ; 
    ep = pop_epoch(chandata,{'S 11'},[-.5,2.5]) ; 
    freqs = 1:2:100 ;
    freqdata = zeros(length(freqs),size(ep.data,1),size(ep.data,2),size(ep.data,3)) ; 
    for i=1:64
        for j=1:length(freqs)
            freqdata(j,i,:,:) = eegfiltfft(squeeze(ep.data(i,:,:))',ep.srate,freqs(j)-2,freqs(j)+2)' ; 
        end
    end
    clear logfreqs fdiffs ; 
    logfreqs = abs(freqdata) ; 
    fdiffs = logfreqs - repmat(mean(logfreqs(:,:,1:50,:),3),[1,1,size(logfreqs,3),1]) ; 
    figure,subplot(2,2,1) ; 
    imagesc(squeeze(mean(mean(fdiffs(:,:,150:550,:),3),4))') ; 
    m = squeeze(mean(mean(fdiffs(:,:,150:550,:),3),4)) ; 
    subplot(2,2,2) ; topoplot(m(8,:),ep.chanlocs) ; 
    subplot(2,2,3) ; topoplot(mean(m(30:40,:)),ep.chanlocs) ; 
    suptitle(subs{sub}) ; 
    
    
end




