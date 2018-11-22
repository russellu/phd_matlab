cd(['c:/shared/badger_eeg/','valerie']) ; 
elecs = [59,45,31,46,60,9,20,18] ; 
restpow = load('post_restpow') ; restpow = restpow.restpow ; 


imagesc(squeeze(mean(restpow(elecs,35:end,:),1))) ; 

motion = squeeze(mean(mean(restpow(elecs,35:end,:),1),2)) ; 

mrestpow = squeeze(mean(restpow(elecs,:,:),1)) ; 

betas = motion\mrestpow' ; 

b = glmfit(motion',mrestpow(1,:)) ; 