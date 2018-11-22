subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','lisa','marc','marie','mathieu',...
    'maxime','mingham','patricia','po','russell','sunachakan','tah','thititip','vincent'} ; 


for s=2:length(subs)
    cd(['C:\shared\all_white_normals\fmris\sub_',subs{s}]) ; 
    cd trigs ; stimTimes = load('stimTimes1.mat') ; cd ..
    stimTimes = stimTimes.stimTimes ; 
    stimTimes = cell2mat(stimTimes) ; 
    times = stimTimes(1:2:end) ; trs = round(times/2) ; 
    fmris = dir('common_*') ; 
    clear allfmris 
    for c=1:length(fmris)
        meanfmri = load_untouch_nii(fmris(c).name) ; 
        fimg = meanfmri.img ; 
        epochs = zeros(size(fimg,1),size(fimg,2),size(fimg,3),30,10) ; 
        for i=1:length(trs)
            epochs(:,:,:,i,:) = fimg(:,:,:,trs(i)-2:trs(i)+7) ; 
        end

        resepochs = reshape(epochs,[size(epochs,1)*size(epochs,2)*size(epochs,3),size(epochs,4),size(epochs,5)]) ; 
        [h,p,ci,stats] = ttest(squeeze(mean(resepochs(:,:,5:6),3))',squeeze(resepochs(:,:,3))') ; 
        rest = reshape(stats.tstat,size(fimg(:,:,:,1))) ; 
        allfmris(:,:,:,c) = rest ; 
    end
    
    
    f1 = load_untouch_nii('f1.nii.gz') ; f1.img = mean(allfmris,4) ; save_untouch_nii(f1,'allttests.nii.gz') ; 

end
