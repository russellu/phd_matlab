clear all ; close all ; 
subs = {'alex','dina','genevieve','jeremie','russell','sukhman','tegan','valerie'}; 
scans = {'bp_reg_topup_mc_retino_allstims_01','bp_reg_topup_mc_retino_allstims_02','bp_reg_topup_mc_retino_gamma_01','bp_reg_topup_mc_retino_gamma_02','bp_reg_topup_mc_retino_movie','bp_reg_topup_mc_retino_rest'}; 
eegscans = {'retino*allstim*01*set','retino*allstim*02*set','retino*gamma*01*set','retino*gamma*02*set','retino*movie*set','retino*rest*set'}; 
comps = {[28,33,37],[18,46],[10,25,27],[11,32],[13,25],[28,39],[44,48,49],[35,46]};

for sb=1:length(subs)
    
    for sc=1:length(scans)
        %cd(['E:\rawbadger\badger_mri\',subs{sb},'\nii\warped']);   
        cd(['E:\badger_eeg\',subs{sb}]);  
        eegscan = dir(['bcg_',eegscans{sc}]); 
        disp(eegscan.name); 
        eeg = pop_loadset(eegscan.name); 
        cleaned{sc} = eeg; 
        %cleaned{sc} = final_bcg(eeg); 
        if sc==1 ; merged = cleaned{sc} ; else merged = pop_mergeset(cleaned{sc},merged); end

    end
    
    save('cleaned','cleaned'); 
    
    
    mergefilt = eegfiltfft(merged.data,merged.srate,40,80); 
    [weights,sphere] = runica(mergefilt(:,1:3:end),'maxsteps',128);
    winv = pinv(weights*sphere);  
    
    highweights{1} = weights; highweights{2} = sphere ; save('highweights','highweights'); 
    
    gammas = pop_mergeset(cleaned{3},cleaned{4}); 
    gammas.data = weights*sphere*gammas.data; 
    epica = pop_epoch(gammas,{'S  1','S  2'},[-1,6]); 
    for i=1:64
            [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(i,:,:)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                    'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off','timesout',200) ; 
    end
    figure,
    for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-6,6]) ; axis xy ; colormap jet; title(i); end

end






