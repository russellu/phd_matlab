cd c:/shared/resmerged ; 
subs = dir('*') ; subs(1:2) = [] ; 

for s=26
    cd c:/shared/resmerged ; 
    cd(subs(s).name)  ;  
    EEG = pop_loadset('merged.set') ; 
    EEG.data = eegfiltfft(EEG.data,EEG.srate,1,128) ; 
    EEG.data = EEG.data - eegfiltfft(EEG.data,EEG.srate,59,61) ; 
    %temp = EEG.data(1:32,:) ; EEG.data(1:32,:) = EEG.data(33:64,:) ; EEG.data(33:64,:) = temp ; 
    gamma = eegfiltfft(EEG.data,EEG.srate,5,25) ; 
    gameeg = EEG ; gameeg.data = gamma ; 
    epochs = pop_epoch(EEG,{'S 11'},[-1,3]) ; 
    epochs.data = abs(epochs.data) ; 
    range = 1:10:1024 ;
    for i=1:length(range) ; disp(range(i)) ; 
       [h,grid_vals(:,:,i)] = topoplot(double(squeeze(mean(epochs.data(:,range(i),:),3))),EEG.chanlocs) ;  
    end
    clear ersp ;
    for i=1:64
          [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochs.data((i),:,:)),epochs.pnts,[epochs.xmin,epochs.xmax],epochs.srate,0,...
              'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'timesout',1000) ; 
    end
    
    clear grid_vals_l grid_vals_h
    for i=1:1000 ; disp(i) ; 
       [h,grid_vals_l(:,:,i)] = topoplot(double(squeeze(mean(ersp(:,4:12,i),2))),EEG.chanlocs) ;  
       [~,grid_vals_h(:,:,i)] = topoplot(double(squeeze(mean(ersp(:,25:35,i),2))),EEG.chanlocs) ;  

    end

    for i=1:size(grid_vals_l,3)
       subplot(1,2,1) ; gl = grid_vals_l(:,:,i) ; gl = flipud(gl) ; imagesc(gl,[-5,5]) ; 
        
       subplot(1,2,2) ; gl = grid_vals_l(:,:,i) ; gl = flipud(gl) ; imagesc(gl,[-5,5]) ; 
        
       getframe ; 
    end

end

