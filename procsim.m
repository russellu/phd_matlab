cd('c:/EEG traces')
cs = dir('c*') ; clear alleegs allersp
for i=1:max(size(cs)) ; 
    cd('c:/EEG traces') ; cd(cs(i).name) ; ls
    eeg = dlmread('Gamma -  Stimulus.txt') ; 
    alleegs(i,:) = eeg(:,2) ; 
    %[allersp(i,:,:),itc,powbase,times,freqs,erspboot,itcboot] = newtimef(alleegs(i,:),4096,[0,4096],1024,0,...
    %                                           'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'verbose','off','winsize',64,'baseline',NaN) ;    
    
end




