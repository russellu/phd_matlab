cd c:/shared/alleegs ;
subs=dir('*') ;
for sub=1:size(subs,1)
    if sub > 2
        disp(subs(sub).name) ;
        cd(['c:/shared/alleegs/',subs(sub).name]) ; 
        resamps=dir('*set') ; 
        for resamp=1:size(resamps,1)
            disp(resamps(resamp).name) ; 
            currentSet = eeg_checkset(pop_loadset('',resamps(resamp).name)) ;
            currentSet = eeg_checkset(pop_resample(currentSet,300)) ;
            outName = ['newvis_',num2str(resamp),'.set'] ; 
            currentSet = eeg_checkset(pop_saveset(currentSet,outName)) ;
        end
    end
end


