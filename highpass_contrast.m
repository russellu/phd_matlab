cd E:\Russ_contrast

dats = dir('*dat'); % brainvision output,
for dat=1:length(dats)%  get all and merge
EEG = pop_loadbv('.',strrep(dats(dat).name,'.dat','.vhdr')); %load using header
EEG.data(32,:) = rand(1,length(EEG.data))*5; 
if dat==1; merged = EEG; else merged = pop_mergeset(EEG,merged); end 
eegs{dat} = EEG; 
end

mergefilt = eegfiltfft(merged.data,merged.srate,40,80); 
[weights,sphere] = runica(mergefilt,'maxsteps',128); 


newmerged = merged; 
newmerged.data = weights*sphere*merged.data; 
trigs = {'S 11','S 12','S 13','S 14','S 15','S 16','S 17','S 18'};
clear ersp
for i=1:length(trigs);  disp(i); 
    epica = pop_epoch(newmerged,{trigs{i}},[-2,8]); 
    for j=1:64
            [ersp(i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(j,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
    end
end

comps = [30,34];
mersp = squeeze(mean(mean(ersp(:,comps,:,times>0.5 & times < 7),2),4));

save('mersp','mersp'); 
save('freqs','freqs');


cd E:\contrast_fmris
fmris = dir('desp*gz');
for fmri=1:length(fmris)
   f = load_untouch_nii(fmris(fmri).name); 
   allfmris(fmri,:,:,:,:) = f.img;    
end

clear trigts percts
for i=1:length(eegs)
    eegi = eegs{i}; 
    lats = cell2mat({eegi.urevent.latency}); 
    atrigs = {eegi.urevent.type}; 
    s80 = find(strcmpi('S 80',atrigs)); 
    s80lat = lats(s80); 
    triglats = zeros(8,3); 
    for tr=1:length(trigs)
        inds = find(strcmpi(trigs{tr},atrigs)); 
        triglats(tr,:) = lats(inds); 
    end
    triglats = triglats - s80lat; 
    ts = zeros(1,round(eegi.srate*(555*0.9))); 
    
    for x=1:size(triglats,1)
        for y=1:size(triglats,2)
            ts(triglats(x,y):triglats(x,y)+eegi.srate*7) = 1; 
        end
    end
    ts = imresize(ts,[1,555]);
    conved = conv(ts,spm_hrf(0.9),'full'); 
    conved = conved(1:length(ts)); 
    
    fmri = squeeze(allfmris(i,:,:,:,:)); 
    resfmri = reshape(fmri,[64*64*36,555]); 
    filtfmri = eegfiltfft(resfmri,1/0.9,0.005,1); 
    resfilt = reshape(filtfmri,size(fmri));
    
    corrs = voxcorr(resfilt(:,:,:,30:end-30),conved(30:end-30)); corrs(isnan(corrs)) = 0; 
    [sv,si] = sort(corrs(:),'descend');
    ts = mean(filtfmri(si(1:100),:)); 
    pts = mean(resfmri(si(1:100),:)); 
    triglats = round(triglats/(eegi.srate*0.9)); 
    
    for x=1:size(triglats,1)
        for y=1:size(triglats,2)
            trigts(i,x,y,:) = ts(triglats(x,y)-5:triglats(x,y)+22);           
            percts(i,x,y,:) = pts(triglats(x,y)-5:triglats(x,y)+22);           
        end
    end

end


btrigts = (percts - repmat(mean(percts(:,:,:,4:7),4),[1,1,1,28])) ./repmat(mean(percts(:,:,:,4:7),4),[1,1,1,28]) ; 
figure,bar(squeeze(mean(mean(mean(trigts(:,:,:,15:20),1),3),4))) ; ylim([0.06,0.1]); 



mbtrigts = squeeze(mean(mean(percts,1),3)); 

percs = (mbtrigts - repmat(mean(mean(mbtrigts(:,4:7))),[8,28])) ./ repmat(mean(mean(mbtrigts(:,4:7))),[8,28]); 


contrasts = {'25%','35%','45%','55%','65%','75%','85%','95%'};
bar(squeeze(mean(percs(:,15:18),2))) ; xlabel('contrast') ; ylabel('%change'); 
set(gca,'XTickLabel',contrasts);  title('BOLD');

save('percs','percs'); 



