
clear all ; close all ; 
cd c:/shared/Raw/subjects ; ls 
subs = dir('S6*vhdr') ; 
for ss=1:length(subs) ; 

sub = pop_loadbv('.',subs(ss).name) ; 
sub = pop_chanedit(sub,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 
%temp = sub.data(1:32,:) ; sub.data(1:32,:) = sub.data(33:64,:) ; sub.data(33:64,:) = temp ; 

sub = pop_resample(sub,256) ; 
sumdiff = find(zscore(sum(diff(sub.data,1,2).^2,2))>2) ; 
hrind = find(sumdiff==32) ; sumdiff(hrind) = [] ;
sub = pop_interp(sub,sumdiff,'spherical') ; 
sumdiff = find(zscore(sum(diff(sub.data,1,2).^2,2))>2) ; 
hrind = find(sumdiff==32) ; sumdiff(hrind) = [] ; 
sub = pop_interp(sub,sumdiff,'spherical') ; 
sub.data = sub.data-eegfiltfft(sub.data,sub.srate,59,61) ; 
sub.data = eegfiltfft(sub.data,sub.srate,0.5,128) ; 
sub = pop_runica(sub,'runica')  ;
sub = pop_saveset(sub,['ica_',strrep(subs(ss).name,'.vhdr','')]) ;
topoplot(sub.icawinv(:,1),sub.chanlocs);
temp = sub.icawinv(1:32,:) ; sub.icawinv(1:32,:) = sub.icawinv(33:64,:) ; sub.icawinv(33:64,:) = temp ; 
%{
trigs = {'S  2','S  8','S 32','S128','S  1'} ; 
for i=1:length(trigs) ; 
    epochsi = pop_epoch(sub,{trigs{i}},[-.5,3]) ; 
    for j=1:size(epochsi.data,1)
        [ersp(ss,i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochsi.data(j,:,:)),epochsi.pnts,[epochsi.xmin,epochsi.xmax],epochsi.srate,0,...
                'plotitc','off','plotersp','off','baseline',0,'freqs',[1,128],'nfreqs',64,'winsize',64) ; 
    end
end
%}
%{
clear ersp 
for i=1:length(trigs) ; 
    epochsi = pop_epoch(sub,{trigs{i}},[-.5,3]) ; 
    for j=1:size(epochsi.icaact,1)
        for k=1:size(epochsi.icaact,3)
            [ersp{i}(j,k,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epochsi.icaact(j,:,k)),epochsi.pnts,[epochsi.xmin,epochsi.xmax],epochsi.srate,0,...
                'plotitc','off','plotersp','off','baseline',NaN,'freqs',[1,128],'nfreqs',64,'winsize',64) ; 
        end
    end
end
save(['ersp_S',num2str(ss)],'ersp','-v7.3') ; 
%}
end




for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(mean(ersp(:,1,i,:,:),1))-squeeze(mean(ersp(:,2,i,:,:),1)),[-3,3]) ; end

for e=1:64
    for i=1:64
        for j=1:200
            [h,p,cl,ts] = ttest(squeeze(ersp(:,1,e,i,j)),squeeze(ersp(:,2,e,i,j))) ; 
            tmats(e,i,j) = ts.tstat ; 
        end
    end
end

for i=1:32 ; subplot(4,8,i) ; topoplot(double(squeeze(mean(tmats(:,i,times>.5 & times<1),3))),sub.chanlocs,'maplimits',[-4,4]) ; title(['hz = ',num2str(freqs(i))]) ; end
suptitle('t-test congruent minus non-congruent, power differences from 0-1 second following trial onset (S2 and S8)') ;

hz = find(freqs>12 & freqs<16) ; icount = 1 ;
for i=1:5:200-5 ;
   subplot(4,10,icount) ; topoplot(double(squeeze(mean(mean(tmats(:,hz,i:i+5),2),3))),sub.chanlocs,'maplimits',[-4,4]) ;  title(['time=',num2str(times(i)),'sec']) ; 
   icount = icount + 1 ;  
end
suptitle('t-test congruent minus non-congruent, in 12-16Hz as a function of time relative to trial onset') ; 







