function tp(EEG)
% plot topoplot
    for i=1:64 ; 
        subplot(5,13,i) ; 
        topoplot(EEG.icawinv(:,i),EEG.chanlocs) ; title(i) ;       
    end

end