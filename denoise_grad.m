function EEG2 = denoise_grad(EEG) 
% function EEG2 = denoise_grad(EEG) 
% REMOVE GRADIENT artifacts from an EEG recording, based on 'R128' as MR
% trigger
% uses average artifact subtraction method (AAS) 
% must be run individually on each set, due to the way it removes time
% points at the end of the denoising process (any time points before/after
% the first and last gradient triggers)

mrindices = find(strcmp({EEG.urevent.type},'R128')) ; 
lats = cell2mat({EEG.urevent.latency}) ;
gradlats = lats(mrindices) ; 

EEG2 = EEG ; 
for chan=1:64 ;
diffgradlats = diff(gradlats) ;
gradepochs = zeros(length(gradlats),diffgradlats(2)) ; 
for i=1:length(gradlats)
        if gradlats(i)+diffgradlats(1) < size(EEG.data,2) ; 
            gradepochs(i,:) = EEG.data(chan,gradlats(i):gradlats(i)+diffgradlats(1)-1) ; 
        end   
end
mepochs = squeeze(mean(gradepochs)) ; 
for i=1:length(gradlats)-1
     EEG2.data(chan,gradlats(i):gradlats(i)+diffgradlats(2)-1) =  EEG.data(chan,gradlats(i):gradlats(i)+diffgradlats(2)-1)-mepochs ; 
end
end

EEG2 = pop_select(EEG2,'point',[gradlats(2),gradlats(length(gradlats)-1)]) ; 

end




