clear all ; close all ; 
cd c:/shared ; 
sources  = load('MN_EEG_KERNEL_170129_2247_sources.mat') ; 
srcdata = sources.ImageGridAmp(:,5000:20000) ; 
%[s,f] = spectopo(sources.ImageGridAmp(:,5000:50000),0,250,'plot','off') ; 
gridlocs = sources.GridLoc ;
sub = load('C:\brainstorm_db\Protocol01\anat\alex\subjectimage_T1.mat') ; 

voxlocs = cs_convert(sub,'scs','voxel',gridlocs) ; 
voxlocs = round(voxlocs) ; 
clear sources ; 
filtsrc = eegfiltfft(srcdata,250,8,25) ; 
filtabs = imresize(abs(filtsrc),[size(filtsrc,1),500]) ; 


%sources = load('data_S__1_trial014.mat_sources.mat') ; 
%sdat = sources.ImageGridAmp ; 
%{
merged = pop_loadset('C:\shared\simdenoise\alex\merged.set') ; 
epi = pop_epoch(merged,{'S  1'},[-1,6]) ; 
ersp = zeros(size(sdat,1),60,200) ; 
for i=1:size(sdat,1) ; disp(i) ; 
[ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(sdat(i,:),epi.pnts+1,[epi.xmin,epi.xmax],epi.srate,0,...
    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off') ; 
end
bersp = ersp - repmat(mean(ersp(:,:,times<0),3),[1,1,200]) ; 
persp = permute(bersp,[3,2,1]) ; 
%}
