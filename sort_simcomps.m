clear all ; close all; 
subs = {'badger_alex','badger_dina','badger_genevieve','badger_jeremie','badger_russell','badger_sukhman','badger_tegan','badger_valerie'}; % 'badger_karl',
viscomps = {[14,26,15,],[8,1,8,19,21],[5,7,11,26,23],[5,20],[10,13,29],[6,11,13,19],[7,18,19,21],[7,12]}; % [9,13,18],
dmncomps = {[5,8],[9,11],[22,30],[9,17],[7,11],[15,18],[8,14],[23,24]}; % [15,25],
bcgcomps = {[1,2,4],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,3],[1,2,3],[1,2,3]}; %[1,2,3],

postelecs = [59,45,31,46,60,9,20,10];

for sb=1:length(subs)
cd(['E:\eegs\',subs{sb}])

sets = dir('retino*set');
%{
for st=1:length(sets)
    eeg = pop_loadset(sets(st).name);
    if st==1; merged = eeg;  else merged = pop_mergeset(merged,eeg); end 
end
%}
%filtmerged = merged; filtmerged.data = eegfiltfft(merged.data,merged.srate,4,60); 
%[weights,sphere] = runica(filtmerged.data,'maxsteps',128); 
%winv = pinv(weights*sphere); 
%figure,for i=1:64;  subplot(5,13,i) ; topoplot(winv(:,i),merged.chanlocs); end ; suptitle(subs{sb}); 
%poscomps{1} = weights; poscomps{2} = sphere; save('poscomps','poscomps'); 

poscomps = load('poscomps'); poscomps = poscomps.poscomps; 
weights = poscomps{1}; sphere = poscomps{2};  
winv = pinv(weights*sphere); 
setnames = {'retino_allstim*01*set','retino_allstim*02*set','retino_gamma*01*set','retino_gamma*02*set','retino_movie*set','retino_rest*set'};
niinames = {'retino_allstim*01*gz','retino_allstim*02*gz','retino_gamma*01*gz','retino_gamma*02*gz','retino_movie*gz','retino_rest*gz'};
saveniinames = {'retino_allstim','retino_allstim','retino_gamma','retino_gamma','retino_movie','retino_rest'};

settimes = {[-1,12],[-1,12],[-1,6],[-1,6]};
allstim_trigs = {'S 11','S 12','S 13','S 14','S 21','S 22','S 23','S 24','S 31','S 32','S 33','S 34','S 41','S 42','S 43','S 44','S 51','S 52','S 53','S 54','S 61','S 62','S 63','S 64','S 71','S 72','S 73','S 74','S 81','S 82','S 83','S 84'};
gamma_trigs = {'S  1','S  2','S  3'};
settrigs = {allstim_trigs,allstim_trigs,gamma_trigs,gamma_trigs}; 
epoch_stim = [10000,10000,5000,5000];
clear corrs; 
freqs = 1:2:100; 
clear allrestrim; 

for st=1:length(setnames)
    cd(['E:\eegs\',subs{sb}])
    eegname = dir(setnames{st});
    eeg = pop_loadset(eegname.name); 
    acts = weights*sphere*eeg.data; 
    bads = zeros(1,64);  bads(dmncomps{sb}) = 1; bads = find(bads==0); %bads(viscomps{sb}) = 1;
    acts(bads,:) = 0; 
    invacts = winv*acts; 
    
    lats = cell2mat({eeg.urevent.latency});
    types = {eeg.urevent.type};
    r128s = find(strcmpi('R128',types)); 
    xvec = acts(1,lats(r128s(1)):lats(r128s(end))); 
    
    allfilt = zeros(64,length(freqs),size(acts,2));       
    for f=1:length(freqs)
       eeg_f = eegfiltfft(invacts,eeg.srate,freqs(f)-1.5,freqs(f)+1.5);  
       allfilt(:,f,:) = eeg_f;   
    end
    
    trimfilt = abs(allfilt(:,:,lats(r128s(1)):lats(r128s(end)))); 
    clear restrim
    for i=1:64
       restrim(i,:,:) = imresize(squeeze(trimfilt(i,:,:)),[size(trimfilt,2),length(r128s)]); 
    end
    hrf = spm_hrf(0.693); 
    
    clear nonconved
    for i=1:64
        for j=1:50
            conved = conv(squeeze(restrim(i,j,:)),hrf,'full'); 
            conved = conved(1:size(restrim,3)); 
            nonconved(i,j,:) = squeeze(restrim(i,j,:)); 
            restrim(i,j,:) = conved; 
        end
    end
    
    allrestrim{st} = restrim; 
    
    %{
    cd(['E:\fmris\',subs{sb},'\atlas_fmri']);
    fmriname = dir(['bp_clean_',niinames{st}]); 
    fmri = load_untouch_nii(fmriname(1).name); 
    for i=1:50; disp(i); 
       corrs(:,:,:,i,st) = voxcorr(fmri.img(:,:,:,50:size(restrim,3)-50),squeeze(mean(restrim(postelecs,i,50:end-50),1)));  
    end
    
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
    
    mrestrim = squeeze(mean(restrim(postelecs,:,:),1)); 
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
    
    %{
    tempeeg = eeg;
    clear alldiffs; 
    for f=1:length(freqs)
        tempeeg.data = squeeze(allfilt(:,f,:)); 

        tempeeg_ep = pop_epoch(tempeeg,settrigs{st},settimes{st}); 
        elec_diffs = mean(abs(tempeeg_ep.data(:,tempeeg_ep.times>0 & tempeeg_ep.times<epoch_stim(st),:)),2) - mean(abs(tempeeg_ep.data(:,tempeeg_ep.times<0,:)),2); 
        alldiffs(:,f,:) = elec_diffs;           
    end    
    stdiffs(st,:,:) = squeeze(mean(alldiffs,3)) ;
    %}
end

pos_restrim = allrestrim; 
save('pos_restrim','pos_restrim'); 

%subxcorrs(sb,:,:) = xcorrs; 

%{
mkdir niicorrs;
freq_ts = load_untouch_nii('ref50.nii.gz'); 
cd niicorrs; 
for i=1:6 ; freq_ts.img = squeeze(corrs(:,:,:,st)); save_untouch_nii(freq_ts,[saveniinames{i},'.nii.gz']); end
freq_ts.img = squeeze(mean(corrs,5)); save_untouch_nii(freq_ts,'mean_states.nii.gz'); 
%}
%subdiffs(sb,:,:,:) = stdiffs; 


%{
gamma1 = pop_loadset('retino_gamma_01.set'); 
gamma2 = pop_loadset('retino_gamma_02.set'); 
gamma_merged = pop_mergeset(gamma1,gamma2); 
gamma_merged.data = weights*sphere*gamma_merged.data; 
triggers = {'S  1','S  2','S  3'};
allep = pop_epoch(gamma_merged,triggers,[-1,6]);
clear stersp1; 
for i=1:64
    [stersp1(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
            'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
end
%figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(stersp(i,:,:)),[-3,3]) ; axis xy;  end

allstim1name = dir('retino_allstim*set'); 

allstims1 = pop_loadset(allstim1name(1).name); 
allstims2 = pop_loadset(allstim1name(2).name); 
allstims_merged = pop_mergeset(allstims1,allstims2); 
allstims_merged.data = weights*sphere*allstims_merged.data; 
triggers = {'S 11','S 12','S 13','S 14','S 21','S 22','S 23','S 24','S 31','S 32','S 33','S 34','S 41','S 42','S 43','S 44','S 51','S 52','S 53','S 54','S 61','S 62','S 63','S 64','S 71','S 72','S 73','S 74','S 81','S 82','S 83','S 84'};
allep = pop_epoch(allstims_merged,triggers,[-1,12]);
clear stersp2;
for i=1:64
    [stersp2(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data(i,:,:)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
            'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
end
%figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(stersp(i,:,:)),[-2,2]) ; axis xy; colormap jet;  end ;suptitle(subs{sb}); 
figure,
for i=1:64 
    subplottight(6,22,i*2) ;
    imagesc(squeeze(stersp1(i,:,:)+stersp2(i,:,:)),[-3,3]) ; axis xy ; colormap jet ; set(gca,'XTick',[],'YTick',[]); 
    subplottight(6,22,i*2-1) ;
    topoplot(winv(:,i),eeg.chanlocs); 
    text(0,.75,num2str(i)); 
end
%}

end
