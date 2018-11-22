clear all ; close all; 
subs = {'badger_alex','badger_dina','badger_genevieve','badger_jeremie','badger_russell','badger_sukhman','badger_tegan','badger_valerie'}; % 'badger_karl',
postelecs = [59,45,31,46,60,9,20,10];
niinames = {'retino_allstim*01*gz','retino_allstim*02*gz','retino_gamma*01*gz','retino_gamma*02*gz','retino_movie*gz','retino_rest*gz'};
saveniinames = {'retino_allstim_01','retino_allstim_02','retino_gamma_01','retino_gamma_02','retino_movie','retino_rest'};

for sb=1:length(subs)
    
    cd(['E:\eegs\',subs{sb}]);
    nonconved = load('nonconved'); nonconved = nonconved.nonconved; 
    pos_nonconved = load('pos_nonconved'); pos_nonconved = pos_nonconved.pos_nonconved; 
    pos_restrim = load('pos_restrim'); pos_restrim = pos_restrim.pos_restrim; 
    allrestrim = load('allrestrim'); allrestrim = allrestrim.allrestrim; 
    
    disp(subs{sb}); 
    
    clear corrs; 
    for st=1:length(niinames) 
        
        eeg_power = allrestrim{st}; 
        
        cd(['E:\fmris\',subs{sb},'\atlas_fmri']);
        fmriname = dir(['bp_clean_',niinames{st}]); 
        fmri = load_untouch_nii(fmriname(1).name); 
        for i=1:50; disp(i); 
           corrs(:,:,:,i,st) = voxcorr(fmri.img(:,:,:,50:size(eeg_power,3)-50),squeeze(mean(eeg_power(postelecs,i,50:end-50),1)));  
        end

        %{
        cd(['E:\fmris\',subs{sb},'\']);
        viscorrs = load_untouch_nii('atlas_gamma_mcorrs.nii.gz'); 
        viscorrs.img(isnan(viscorrs.img)) = 0 ;
        [sv,si] = sort(viscorrs.img(:),'descend'); 
        resfmri = reshape(fmri.img,[numel(fmri.img(:,:,:,1)),size(fmri.img,4)]);
        v1vox = resfmri(si(1:250),:); 
        tcorrs = zeros(1,size(v1vox,2)); 
        for t=25:size(v1vox,2)-25
           tcorrs(t) = mean(mean(corr(v1vox(:,t-10:t+10)')));         
        end
 
        mrestrim = squeeze(mean(eeg_power(postelecs,:,:),1)); 
        mnonconved = squeeze(mean(nonconved(postelecs,:,:),1));
        clear xcorrs; 
        for i=1:50
           %xcorrs(i,:) = xcorr(mrestrim(i,50:end-50),tcorrs(50:end-50),20,'coeff') ;
           % xcorrs(i,:) = crosscorr(mrestrim(i,50:end-50),tcorrs(50:end-50)) ;
           jcount = 1;
            for j=-20:20
                xcorrs(i,jcount) = corr2(squeeze(mrestrim(i,50:end-50)),tcorrs(50+j:end-50+j)); 
                jcount = jcount + 1; 
            end
        end
    %}
    end
    
    allcorrs(:,:,:,:,:,sb) = corrs; 

end


%cd e:/coupling_saved; 
%ref = load_untouch_nii('ref50.nii.gz'); 
%for i=1:6 ; ref.img = squeeze(mean(allcorrs(:,:,:,:,i,:),6)); save_untouch_nii(ref,['gsr_',saveniinames{i},'.nii.gz']); end

fcorrs = squeeze(mean(allcorrs(:,:,:,4:8,:,:),4)); 
taskinds = {[1,2],[3,4],[5],[6],1:6};
taskindlabs = {'block','event-related','continuous','rest','grand avg'};
cd e:/fmris/badger_russell ; 
f1 = load_untouch_nii('f1.nii.gz'); 
mask = f1.img > (mean(f1.img(:))*2); 
f1.img = f1.img.*mask; 

for i=1:length(taskinds)
   subplot(2,7,i); 
   
   corrmapi = squeeze(mean(mean(mean(fcorrs(:,:,10:20,taskinds{i},:),3),4),5)); 
   corrmapi(isnan(corrmapi)) = 0; corrmapi(corrmapi>0.08) = 0.08; corrmapi(corrmapi<-0.08) = -0.08; corrmapi(1,1) = -0.08; corrmapi(1,2) = 0.08; 
   meanmap = squeeze(mean(f1.img(:,:,11:19),3)); 
   plotoverlayIntensity2D(meanmap,mat2gray(abs(corrmapi)),corrmapi,270); title(taskindlabs{i}); 
   
   subplot(2,7,i+7); 
   corrmapi = squeeze(mean(mean(mean(fcorrs(25:40,:,:,taskinds{i},:),1),4),5)); 
   corrmapi(isnan(corrmapi)) = 0; corrmapi(corrmapi>0.08) = 0.08; corrmapi(corrmapi<-0.08) = -0.08; corrmapi(1,1) = -0.08; corrmapi(1,2) = 0.08; 
   meanmap = squeeze(mean(f1.img(28:36,:,:),1)); 
   plotoverlayIntensity2D(meanmap,mat2gray(abs(corrmapi)),corrmapi,90);
   
   
end
subplot(2,7,7) ; imagesc([-.08,.08]) ; h = colorbar ; title(h,'r'); colormap jet; 






