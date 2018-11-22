subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 
comps = {[30,26,20,45],[32,23],[21,28],[47,46],[19,37],[47,38],[44,46]} ; % all subjects, right and left
%comps = {[4,5,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,24,26,27]} ; 

highfreq_allcomps = {[20,26,30,45],[23,32],[21,28,44],[47,46],[19,37],[38,47],[24,44,46]} ;
allfreq_allcomps = {[11,15,22,26,51],[6,14,16,17],[7,10,11,13,19],[7,11,22],[13,17,23],[7,17,22,25],[8,23,24]} ; 

for s=1:length(subs)
    cd(['c:/shared/badger_eeg/',subs{s}]) ;  
    ls 
    gammas=dir('allfreq*gamma*set') ; 
    for g=1:length(gammas)
        EEG = pop_loadset(gammas(g).name) ; 
        if g==1 ; merged = EEG ; else merged = pop_mergeset(merged,EEG,1) ; end
        
        
    end
    bads = zeros(1,64) ; bads(allfreq_allcomps{s}) = 1 ; bads = find(bads==0) ; 
    filtered = pop_subcomp(merged,bads,0) ; 
    s1 = pop_epoch(filtered,{'S  1','S  2'},[-5,7]) ;
    
    merges1 = pop_epoch(merged,{'S  1','S  2'},[-5,7]) ;
    % compute ersp
    %{
    for i=1:64 ;
       [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(merges1.icaact(i,:,:)),merges1.pnts,[merges1.xmin,merges1.xmax],merges1.srate,0,...
                                                            'plotersp','off','plotitc','off','verbose','off','freqs',[1,100],'nfreqs',50,'winsize',64,'baseline',0) ;      
    end
    for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(ersp(i,:,:)),[-5,5]) ; title(i) ; end
    %}
    

    % compute raw power
    clear freqpower freqs ; 
    for i=1:64
        jcount = 1 ; 
        for j=1:2:100
            freqpower(i,jcount,:,:) = eegfiltfft(squeeze(s1.data(i,:,:))',s1.srate,j-1,j+1) ; 
            freqs(jcount) = j ; 
            jcount = jcount + 1 ;  
        end
    end
    stimtinds = find(s1.times>0 & s1.times<5) ; 
    meanpow = squeeze(mean(mean(freqpower(:,:,:,stimtinds).^2,3),4)) ; 
    basepow = squeeze(mean(abs(freqpower(:,:,:,s1.times<0 & s1.times>-3000)),4)) ; stimpow = squeeze(mean(abs(freqpower(:,:,:,s1.times>0 & s1.times<5000)),4)) ; 
    powdiff = stimpow-basepow ; 
    basecorr = mean(powdiff,3) ; 
    figure,for i=1:50 ; subplot(5,10,i) ; topoplot(squeeze(basecorr(:,i)),s1.chanlocs) ; title(i) ; end
    
   % save('basecorr','basecorr') ; 

end
