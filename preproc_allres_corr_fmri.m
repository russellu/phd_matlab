clear all ; close all ;
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','lisa','marc','marie','mathieu',...
    'maxime','mingham','patricia','po','russell','sunachakan','tah','thititip','vincent'} ; 

comps = {[18,13],[13,10],[10,8],[13,7],[12,23],[6,15],[18,8],[27,30],[16,20],[18,9],[10,9],[29,15],[6,8],[10,12],...
    [14,36],[13,8],[34,19],[8,16],[21,20],[9,13],[23,16],[17,12],[18,7],[11,25]} ;

meanpos = load('C:\shared\all_white_normals\fmris\sub_vincent\meanpos') ; meanpos = meanpos.meanpos ; 
meanneg = load('C:\shared\all_white_normals\fmris\sub_vincent\meanneg') ; meanneg = meanneg.meanneg ; 

for s=1:length(subs)
    cd(['C:\shared\all_white_normals\fmris\sub_',subs{s}]) ; 

    fmris = dir('common_*') ; 
    cd melodic ; ic = load_untouch_nii('melodic_IC.nii.gz') ; icimg = ic.img ; 
    posroi = icimg(:,:,:,comps{s}(1)) ; negroi = icimg(:,:,:,comps{s}(2)) ; cd .. 
    negroi(:,1:round(size(negroi,2)/2),:) = 0 ; posroi(:,1:round(size(posroi,2)/2),:) = 0 ;  

    [sv,si] = sort(negroi(:),'descend') ; 
    zneg = zeros(size(negroi)) ; zneg(si(1:100)) = 1 ; 
    
    [sv,si] = sort(posroi(:),'descend') ; 
    zpos = zeros(size(posroi)) ; zpos(si(1:200)) = 1 ; 
    clear img ; figure,
    img(:,:,1) = uint8(mat2gray((sum(zpos,3)))*255) ; 
    img(:,:,3) = uint8(mat2gray((sum(zneg,3)))*255) ; 
    imagesc(img(20:end-20,40:end,:)) ; title(num2str(s)) ; 

    clear allfmris allposcorrs allnegcorrs negepochs posepochs
    for c=1:length(fmris)
        cd trigs ; 
        stimTimes = load(['stimTimes',num2str(c),'.mat']) ; cd ..
        stimTimes = stimTimes.stimTimes ; 
        stimTimes = cell2mat(stimTimes) ; 
        times = stimTimes(1:2:end) ; trs = round(times/2) ; 
        types = stimTimes(2:2:end) ; 
        
        meanfmri = load_untouch_nii(fmris(c).name) ; 
        fimg = meanfmri.img ; 
        resimg = reshape(fimg,[numel(fimg(:,:,:,1)),size(fimg,4)]) ; 
        
        posinds = find(zpos==1) ; neginds = find(zneg==1) ; 
        pos_ts = resimg(posinds,:) ; 
        neg_ts = resimg(neginds,:) ; 
        
        for i=1:6
            typesi = find(types==i) ; 
            trsi = trs(typesi) ; 
            for j=1:length(typesi)
               epochinds = trsi(j)-2:trsi(j)+8 ; 
               negepochs(c,i,j,:,:) = neg_ts(:,epochinds) ; 
               posepochs(c,i,j,:,:) = pos_ts(:,epochinds) ; 
            end
        end
        
    
       %{
        epochs = zeros(size(fimg,1),size(fimg,2),size(fimg,3),30,10) ; 
        for i=1:length(trs)
            epochs(:,:,:,i,:) = fimg(:,:,:,trs(i)-2:trs(i)+7) ; 
        end

        resepochs = reshape(epochs,[size(epochs,1)*size(epochs,2)*size(epochs,3),size(epochs,4),size(epochs,5)]) ; 
        [h,p,ci,stats] = ttest(squeeze(mean(resepochs(:,:,5:6),3))',squeeze(resepochs(:,:,3))') ; 
        rest = reshape(stats.tstat,size(fimg(:,:,:,1))) ; 
        allfmris(:,:,:,c) = rest ; 
        allposcorrs(:,:,:,c) = voxcorr(fimg(:,:,:,35:end-35),meanpos(35:end-35)) ; 
        allnegcorrs(:,:,:,c) = voxcorr(fimg(:,:,:,35:end-35),meanneg(35:end-35)) ; 
        %}
    end
    allpos(s,:,:,:,:,:) = posepochs ; 
    allneg(s,:,:,:,:,:) = negepochs ; 
    %f1 = load_untouch_nii('f1.nii.gz') ; f1.img = mean(allfmris,4) ; save_untouch_nii(f1,'allttests.nii.gz') ; 
    %f1.img = mean(allposcorrs,4) ; save_untouch_nii(f1,'meanpos.nii.gz') ; 
    %f1.img = mean(allnegcorrs,4) ; save_untouch_nii(f1,'meanneg.nii.gz') ; 
   
end


subneg = squeeze(mean(mean(mean(allneg,2),4),5)) ; 
subpos = squeeze(mean(mean(mean(allpos,2),4),5)) ; 
colors = {'b','g','r','c','m','y'}  ;

for i=1:5
    shadedErrorBar([],squeeze(mean(subpos(:,i,:),1)),squeeze(std(subpos(:,i,:),0,1))/sqrt(24),{colors{i}}) ; hold on ; 
end




for s=1:length(subs)
    cd(['C:\shared\all_white_normals\fmris\sub_',subs{s}]) ; 
    poscorrs = load_untouch_nii('meanpos.nii.gz') ; 
    posimg = poscorrs.img ; posimg(:,1:end-40,:) = 0 ; 
    voxsums(s) = sum(posimg(:)>.35) ; 
    subplot(3,8,s) ; imagesc(squeeze(sum(posimg>.3,3))) ; title(subs{s}) ; 
end

pchange = squeeze(mean(mean(mean(mean(allpos,2),5),4),3)) ; 
pchange = (mean(pchange(:,5:6),2)-mean(pchange(:,1:2),2)) ./ mean(pchange(:,1:2),2) ; 

for s=1:length(subs)
    cd(['C:\shared\all_white_normals\fmris\sub_',subs{s}]) ; 
    pchangei = pchange(s) ; 
    save('pchangei','pchangei') ; 
    pvoxi = voxsums(s) ; 
    save('pvoxi','pvoxi') ;
    allposi = squeeze(allpos(s,:,:,:,:,:)) ; 
    save('allposi','allposi') ; 
end






