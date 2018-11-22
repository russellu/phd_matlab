function EEG2 = denoise_grad2(EEG) ;
%cd('C:\Vision\Raw Files\Russell_2015-10-15') ; 
%clear all ; close all ; 
%vhdrs = dir('*vhdr') ; 
%for i=2% :length(vhdrs) ; 
%    EEG = pop_loadbv('.',vhdrs(i).name) ; 
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
    mgradords = squeeze(mean(gradords,1)) ; % average across all TRs
    avgper = 25 ; % 25 period moving average 
    % perform the subtraction
    for gr=1:size(gradords,1)
        if gr+avgper < size(gradords,1) % if not near the end, take the future average
            EEG2.data(:,gradlats(gr):gradlats(gr)+trlen-1) = EEG.data(:,gradlats(gr):gradlats(gr)+trlen-1) - squeeze(mean(gradords(gr:gr+avgper,:,:),1)) ;           
        else % if near the end, take the previous average (bounds checking)
            EEG2.data(:,gradlats(gr):gradlats(gr)+trlen-1) = EEG.data(:,gradlats(gr):gradlats(gr)+trlen-1) - squeeze(mean(gradords(gr-avgper:gr,:,:),1)) ;           
        end
    end
    % resample return value (for some reason reampling after removing
    % points is fucked)
    EEG2 = pop_resample(EEG2,256) ; 
    gradtriginds = find(strcmp('R128',{EEG2.urevent.type})) ; 
    lats = {EEG2.urevent.latency} ; 
    gradlats = cell2mat(lats(gradtriginds)) ; 
    diffgrad = diff(gradlats) ;
    trlen = diffgrad(1)  ; 
    EEG2.data(:,1:gradlats(1)+trlen/2) = 0 ; 
    EEG2.data(:,gradlats(length(gradlats))+trlen/2:end) = 0 ; 

end


