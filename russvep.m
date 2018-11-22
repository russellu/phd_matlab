cd C:\shared\russvep
EEG = pop_loadbv('.','r_ssvep_quad_01.vhdr') ; 
EEG = pop_chanedit(EEG,'lookup','C:\eeglab8_0_3_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ; 

EEG = pop_resample(EEG,256) ; 

filt = EEG ; filt.data = eegfiltfft(EEG.data,EEG.srate,1,128) ; 
filtica = pop_runica(filt,'runica') ; 




ep = pop_epoch(filtica,{'S 11'},[-3,120]) ; 

clear ersp ; 
for c=1:64 ;
   [ersp(c,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.icaact(c,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
       'plotersp','off','plotitc','off','timesout',1000,'freqs',[1,120],'nfreqs',120,'winsize',256) ;      
   
end
for i=1:64 ; ersp(i,:,:) = repmat(max(ersp(i,:,:),[],3),[1,1,1000])-ersp(i,:,:)  ; end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:))) ; title(i) ;end


[xg,yg] = meshgrid(-100:100,-100:100) ; 
[th,rh] = cart2pol(xg,yg) ; 
th = th*180 / pi +90 ; th = mod(th,360) ; 
timeinds = find(times>0 & times<15) ;  
anglestep = 360/length(timeinds) ; 
for i=1:64
    timevals = squeeze(ersp(i,16,timeinds)) ; 
    prf = zeros(size(th)) ; 
    for j=1:length(timevals) ;
        prf((th>(j-1)*anglestep-2 & th < j*anglestep*2)) = timevals(j) ;         
    end
    prfs(i,:,:) = prf ; 
end
figure,for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(prfs(i,:,:))) ; set(gca,'XTick',[],'YTick',[]) ; title(i); end



