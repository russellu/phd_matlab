clear all ; close all; 
subs = {'badger_alex','badger_dina','badger_genevieve','badger_jeremie','badger_russell','badger_sukhman','badger_tegan','badger_valerie'}; %'badger_karl',

%dmn_indices = [14,2,8,1,2,14,19,31,35];
dmn_indices = [15,2,8,1,2,15,18,32,4]; 

for sb=1:length(subs)
cd(['E:\fmris\',subs{sb},'\mel'])
ics = load_untouch_nii('atlas_IC.nii.gz'); 
meanimg = load_untouch_nii('mean.nii.gz'); 
%{
figure,
subplot(2,2,1) ; imagesc(imrotate(squeeze(max(ics.img(:,:,:,dmn_indices(sb)),[],3)),270)); 
subplot(2,2,2) ; imagesc(imrotate(squeeze(max(ics.img(:,:,:,dmn_indices(sb)),[],2)),90)); 
subplot(2,2,3) ; imagesc(imrotate(flipud(squeeze(max(ics.img(:,:,:,dmn_indices(sb)),[],1))),90)); colormap jet; 
%}

dmn = ics.img(:,:,:,dmn_indices(sb)); 
[sv,si] = sort(dmn(:),'descend'); 
cd(['E:\fmris\',subs{sb},'\atlas_fmri'])
allstim = load_untouch_nii('bp_clean_retino_movie.nii.gz'); 
res_f = reshape(allstim.img,[numel(allstim.img(:,:,:,1)),size(allstim.img,4)]); 
dmn_ts = (squeeze(res_f(si(1:1000),:))'); 

clear corrs; 
for i=1:size(dmn_ts,1)-50
   corrs(i) = mean(mean(corr(dmn_ts(i:i+50,:))));  
end


corrs = zeros(1,size(dmn_ts,1)); 
for i=50:size(dmn_ts,1)-50
   corrs(i) = mean(mean(corr(dmn_ts(i-10:i+10,:))));  
end

%{
cd(['E:\eegs\',subs{sb}]);
fullweights = load('fullweights'); fullweights = fullweights.fullweights; 
winv = pinv(fullweights); 
resteeg = pop_loadset('den_retino_rest.set'); 

acts = fullweights*resteeg.data; 

lats = cell2mat({resteeg.urevent.latency});
types = {resteeg.urevent.type};
r128s = find(strcmpi('R128',types)); 
xvec = acts(1,lats(r128s(1)):lats(r128s(end))); 

freqs = 1:2:80; 
allfilts = zeros(64,length(freqs),450); 
convfilts = zeros(64,length(freqs),450); 
for f=1:length(freqs)
    filts = eegfiltfft(acts,resteeg.srate,freqs(f)-1.5,freqs(f)+1.5); 
    allfilts(:,f,:) = imresize(abs(filts(:,lats(r128s(1)):lats(r128s(end)))),[64,450],'bicubic'); 
end
%}
cd(['E:\fmris\',subs{sb},'\den_retino_movie']);
alpha = load_untouch_nii('alphabeta_atlas_hz.nii.gz'); 
hrfalpha = zeros(size(alpha.img)); 
hrf = spm_hrf(0.693); 
for i=1:size(alpha.img,1); disp(i); 
    for j=1:size(alpha.img,2)
        for k=1:size(alpha.img,3)
            convedijk = conv(squeeze(alpha.img(i,j,k,:)),hrf,'full'); 
            hrfalpha(i,j,k,:) = convedijk(1:size(alpha.img,4)); 
        end
    end
end
dmn_ts(isnan(dmn_ts)) = 0 ; hrfalpha(isnan(hrfalpha)) = 0; 
%vcorrs = voxcorr(hrfalpha(:,:,:,51:end-51),mean(dmn_ts(51:end-51,:),2)); 

vcorrs = voxcorr(hrfalpha(:,:,:,51:size(hrfalpha,4)-51),corrs(51:size(hrfalpha,4)-51)); 
dmncorrs = voxcorr(hrfalpha(:,:,:,51:size(hrfalpha,4)-51),mean(dmn_ts(51:size(hrfalpha,4)-51,:),2)); 

allvcorrs(sb,:,:,:) = vcorrs; 
alldmncorrs(sb,:,:,:) = dmncorrs; 

% compare the dynamic coupling to the resting state connectivity? 
% get dynamic correlation of hrf convolved alpha/beta with BOLD, compare
% that to correlation within BOLD DMN

%{
cd(['E:\fmris\',subs{sb},'\den_retino_rest']);
atlas_hzs = dir('atlas_hz*nii'); 
for ahz=1:50; disp(ahz); 
   freqs_i = load_untouch_nii(atlas_hzs(ahz).name);  
end
%}

%}
%{
figure,
sz = round(sqrt(size(ics.img,4)))+1; 
for i=1:size(ics.img,4)
    subplottight(sz,sz,i) ;
    imagesc(imrotate(squeeze(max(ics.img(:,:,:,i),[],3)),270)); set(gca,'XTickLabel',[],'YTickLabel',[]); 
    title(num2str(i)); 
end
%}

end


