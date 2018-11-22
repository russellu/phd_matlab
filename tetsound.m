cd('E:\Données EEG\ERP_Chirp'); 
eeg = pop_loadbv('.','C2_CI_ERP_Chirp.vhdr'); 
eeg = pop_resample(eeg,250);

filteeg = eeg;
filteeg.data = eeg.data - eegfiltfft(eeg.data,eeg.srate,59,61)- eegfiltfft(eeg.data,eeg.srate,0,1); 


allep = pop_epoch(filteeg,{'S  1'},[-.2,.6]);
%for i=1:64 ;subplot(5,13,i) ; plot(eps.times,squeeze(mean(eps.data(i,:,:),3)),'k') ; vline(0,'r'); end

[weight,sphere] = runica(filteeg.data,'maxsteps',128); 
icaeeg = filteeg ; icaeeg.data = weight*sphere*filteeg.data; 
allep = pop_epoch(icaeeg,{'S  1'},[-.2,.6]);

clear allersp itc
for i=1:64 
        [allersp(i,:,:),itc(i,:,:),powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
end

for i=1:64 ;subplot(5,13,i) ; plot(squeeze(mean(allep.data(i,:,:),3)),'k') ; vline(0,'r'); end