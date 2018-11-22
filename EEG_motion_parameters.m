clear all ; close all ; 

subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 

% xcorr bcg with all other channels

for sb=2:length(subs)
    
    cd(['E:\badger_eeg\',subs{sb}]);   
    sets = dir('retino*set'); 
    
    for set=1:length(sets)
        
        eeg = pop_loadset(sets(set).name); 
        raweeg = eeg; 
        trigs = {eeg.urevent.type};
        trtrigs = find(strcmpi('R128',trigs)); 
        lats = {eeg.urevent.latency}; 
        trlats = cell2mat(lats(trtrigs)); 
        
        cleaned = final_bcg(eeg); 
        
        
        
        freqs = 1:2:100; 
        clear res_filts bcg_res_filts    
        for f=1:length(freqs)
           filt_f = eegfiltfft(eeg.data,eeg.srate,freqs(f)-1,freqs(f)+1);  
           filt_f = filt_f(:,trlats(1):trlats(end)); 
           res_filts(f,:,:) = imresize(abs(filt_f),[64,length(trlats)]);  
           
           filt_f = eegfiltfft(cleaned.data,cleaned.srate,freqs(f)-1,freqs(f)+1);  
           filt_f = filt_f(:,trlats(1):trlats(end)); 
           bcg_res_filts(f,:,:) = imresize(abs(filt_f),[64,length(trlats)]);  
        end
        
        for i=1:50
           res_filts(i,:,:) = eegfiltfft(squeeze(res_filts(i,:,:)),1/0.693,0.005,2);  
           bcg_res_filts(i,:,:) = eegfiltfft(squeeze(bcg_res_filts(i,:,:)),1/0.693,0.005,2);  
        end
        
        for i=1:50
           for j=1:64
              xcorrs(sb,set,i,j,:) = xcorr(squeeze(res_filts(i,j,20:end-20)),squeeze(res_filts(i,32,20:end-20)),20,'coeff'); 
              bcg_xcorrs(sb,set,i,j,:) = xcorr(squeeze(bcg_res_filts(i,j,20:end-20)),squeeze(res_filts(i,32,20:end-20)),20,'coeff'); 
       
           end
        end
        
        
        
    end
    
    
    
end



mx = (squeeze(mean(mean(mean(xcorrs,1),2),4))); 
bcg_mx = (squeeze(mean(mean(mean(bcg_xcorrs,1),2),4))); 

times = (-20:20)*0.693; 
subplot(2,2,1) ; imagesc(times,1:2:100,mx,[-.35,.35]) ; axis xy ; colormap jet; 
subplot(2,2,2) ; imagesc(times,1:2:100,bcg_mx,[-.35,.35]) ; axis xy; colormap jet;
subplot(2,2,3); 
plot(1:2:100,squeeze(mean(mean(mean(xcorrs(:,:,:,[1:31,33:end],21),1),2),4))); hold on ; 
plot(1:2:100,squeeze(mean(mean(mean(bcg_xcorrs(:,:,:,[1:31,33:end],21),1),2),4)));
legend({'before subtraction', 'after subtraction'});
%{
sub_mx = squeeze(mean(mean(xcorrs(:,:,:,[1:31,33:end],21),2),4)); 


shadedErrorBar(1:2:100,mean(sub_mx(:,:),1),std(sub_mx(:,:),0,1)/sqrt(9)); 
xlabel('frequency(hz)'); ylabel('correlation'); title('T=0 correlation with BCG (all channels)');


%}