clear all ; close all; 
cd E:\data_for_Russ\Lyes_no_adapt_30_oct

eeg = pop_loadbv('.','LYES_noadaptation.vhdr');
eeg = pop_chanedit(eeg,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
filt = eegfiltfft(eeg.data,eeg.srate,25,60); 
[weights,sphere] = runica(filt,'maxsteps',128); 
epmerged = eeg; epmerged.data = weights*sphere*eeg.data; 
winv = pinv(weights*sphere); 
figure,for i=1:64 ; subplot(5,13,i) ; topoplot(winv(:,i),eeg.chanlocs); end
stims = {'S  2','S  4'};
clear mstersp; 
for i=1:length(stims); disp(i); 
    allep = pop_epoch(epmerged,{stims{i}},[-2,20]); 
    for j=1:64
        for k=1:size(allep.data,3)
            [mstersp(i,j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(allep.data((j),:,k)),allep.pnts,[allep.xmin,allep.xmax],allep.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',80,'baseline',NaN,'verbose','off','timesout',1000) ; 
        end
    end
end
bersp = mstersp - repmat(mean(mstersp(:,:,:,:,times<0),5),[1,1,1,1,1000]);    

%{
allbersp(1,:,:,:) = squeeze(mean(mean(bersp(:,[22,24,27],:,:,:),2),3)); 


time_inds = find(times > 0 & times < 18);    
gtime_inds = find(times>0.5 & times<18);    

oblique = [35:55,125:145];
vertical = 80:100; 
horizontal = [1:10,170:180];
cardinal = [80:100,1:10,170:180]; 

s2_angles = [90:-1:0,359:-1:45];
s4_angles = [45:1:359,0:1:90];

tbersp = allbersp(:,:,:,time_inds); 
clear res_angles
res_angles(1,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
res_angles(2,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 
res_angles(3,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
res_angles(4,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 
res_angles(5,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
res_angles(6,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 
res_angles(7,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
res_angles(8,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 
res_angles(9,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
res_angles(10,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 

tbersp = allbersp(:,:,:,gtime_inds); 
res_angles = res_angles(:,(length(time_inds)-length(gtime_inds)+1):end); 

tbersp_amp = zeros(size(tbersp,1),size(tbersp,2),size(tbersp,3),360);    
for i=0:359
    for j=1:2
        bersp_ind = find(i==res_angles(j,:)); 
        tbersp_amp(:,j,:,i+1) = squeeze(mean(tbersp(:,j,:,bersp_ind),4)); 
    end
end

clear sub_db
sub_db(:,1,:,:) = squeeze(mean(tbersp_amp(:,1:2,:,:),2)); 
%}






time_inds = find(times > 0 & times < 18);    
gtime_inds = find(times>0.5 & times<18);    

s2_angles = [90:-1:0,359:-1:45];
s4_angles = [45:1:359,0:1:90];

tbersp = bersp(:,:,:,:,time_inds); 
clear res_angles
res_angles(1,:) = imresize(s2_angles,[1,length(tbersp)],'nearest'); 
res_angles(2,:) = imresize(s4_angles,[1,length(tbersp)],'nearest'); 

tbersp = bersp(:,:,:,:,gtime_inds); 
res_angles = res_angles(:,(length(time_inds)-length(gtime_inds)+1):end); 

tbersp_amp = zeros(size(tbersp,1),size(tbersp,2),size(tbersp,3),size(tbersp,4),360);    
for i=0:359
    for j=1:2
        bersp_ind = find(i==res_angles(j,:)); 
        tbersp_amp(j,:,:,:,i+1) = squeeze(mean(tbersp(j,:,:,:,bersp_ind),5)); 
    end
end


mamp = squeeze(mean(mean(tbersp_amp(:,[22,24,27],:,:,:),1),2)); 
mamp = mamp(:,:,1:180)/2 + mamp(:,:,181:end)/2;

for i=1:40
   mamp(i,:,:) = imfilter(squeeze(mamp(i,:,:)),fspecial('gaussian',[3,9],5),'circular');  
end

pretrials = 1:20;
posttrials = 21:40; 

subplot(2,2,1) ; imagesc(1:180,freqs,squeeze(mean(mamp(pretrials,:,:),1)),[-5,5]);  axis xy ; colormap jet; title('pre-adaptation'); xlabel('orientation (deg)'); ylabel('frequency(hz)'); 
subplot(2,2,2) ; imagesc(1:180,freqs,squeeze(mean(mamp(posttrials,:,:),1)),[-5,5]); axis xy ; colormap jet ; title('post-adaptation');  xlabel('orientation (deg)'); ylabel('frequency(hz)'); 
subplot(2,2,3); 
plot(squeeze(mean(mean(mamp(pretrials,20:30,:),1),2))); hold on; 
plot(squeeze(mean(mean(mamp(posttrials,20:30,:),1),2))); 

for i=1:180
   [h,p(i),~,~] = ttest(squeeze(mean(mamp(pretrials,20:30,i),2)),squeeze(mean(mamp(posttrials,20:30,i),2)));  
end
for i=1:180 ; if p(i) < 0.05 ; text(i,2,'*'); end; end ; title('*p<0.05'); xlabel('orientation(deg)'); ylabel('gamma power (dB)'); 












