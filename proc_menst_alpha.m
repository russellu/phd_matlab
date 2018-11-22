clear all ; close all; 
cd E:\angelina_2; 
eegs = dir('*_luteal_0*vhdr'); 
for i=1:length(eegs)
eeg = pop_loadbv('.',eegs(i).name) ; 
eeg = pop_chanedit(eeg,'lookup','C:\mscripts2\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp') ;
eeg = pop_resample(eeg,256); 
if i==1; merged = eeg ; else merged=  pop_mergeset(merged,eeg); end
end

merged.data = merged.data - eegfiltfft(merged.data,merged.srate,59,61); 
filtmerged = eegfiltfft(merged.data,merged.srate,1,90); 
trigs = {'S  1'}; 

[weights,sphere]= runica(filtmerged,'maxsteps',128); 
winv = pinv(weights*sphere); 
newmerged = merged; newmerged.data = weights*sphere*merged.data; 
clear ersp
for i=1:length(trigs);  disp(i); 
    epica = pop_epoch(newmerged,{trigs{i}},[-15,20]); 
    for j=1:64
        for k=1:20
            [ersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(j,:,k)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',200) ; 
        end
    end
end
bersp = squeeze(ersp - repmat(mean(ersp(1,:,:,:,times<0),5),[1,1,1,1,200])); 

for i=1:64 ; subplot(5,13,i),imagesc(squeeze(mean(bersp(i,:,:,:),2)),[-8,8]); axis xy ; colormap jet; title(i);  end
comps = [2,3,9]; 

imagesc(times,freqs,squeeze(mean(mean(bersp(comps,:,:,:),2),1)),[-8,8]); axis xy; colormap jet; title('eyes closed vs eyes open'); 

for i=1:20 ; subplot(2,10,i) ; imagesc(times,freqs,squeeze(mean(bersp(comps,i,:,:),1)),[-10,10]); axis xy ; colormap jet; end



