clear all ; close all 
subs = {'alex','dina','genevieve','jeremie','russell','valerie','tegan'} ; 

comps = {[30,26],[32,23],[21,28],[47,46],[19,37],[47,38]} ;

for sub=1:length(subs) 
    cd(['c:/shared/badger_eeg/',subs{sub}]) ; 
    prefix = 'allfreq_' ; 
    gammas1=dir([prefix,'preproc*gamma*set']) ;gammas2=dir([prefix,'preproc*gamma*set']) ;
    allstims1=dir([prefix,'preproc*allstim*set']) ;allstims2=dir([prefix,'preproc*allstim*set']) ;
    movies=dir([prefix,'preproc*movie*set']) ; rests=dir([prefix,'preproc*rest*set']) ;
    setNames = {allstims1(1).name,allstims2(2).name,gammas1(1).name,gammas2(2).name,movies(1).name,rests(1).name} ;   
    
    for setN = 1:6 ;    
        EEG = pop_loadset(setNames{setN}) ; %EEG2 = pop_loadset(setNames{setN+1}) ; EEG = pop_mergeset(EEG,EEG2,1) ; 
        %{
        trigs = {'S  1','S  2','S  3'} ; 
        clear ersp
        for trig=1:3
            ep = pop_epoch(EEG,{trigs{trig}},[-2,7]) ; 
            
            for c=1:64 
                [ersp(trig,c,:,:),itc,powbase,times2,freqs2,~,~] = newtimef(ep.icaact(c,:,:),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                    'plotersp','off','plotitc','off','baseline',0,'timesout',100,'nfreqs',50,'freqs',[1,100],'winsize',64) ; 
            end
            
            %{
            for c=1:length(comps{sub})  
                for trial=1:size(ep.data,3) 
                    [sersp(trig,c,trial,:,:),itc,powbase,times2,freqs2,~,~] = newtimef(ep.icaact(comps{sub}(c),:,trial),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                        'plotersp','off','plotitc','off','baseline',NaN,'timesout',100,'nfreqs',50,'freqs',[1,100],'winsize',64) ; 
                end
            end
            %}      
        end
        %}
        
        [s,f] = spectopo(EEG.icaact(:,20*EEG.srate:end-20*EEG.srate),0,EEG.srate,'plot','off') ; 
        allspec(sub,setN,:,:) = s ; 
        %{
        for c=1:length(comps{sub})  
            for trial=1:size(ep.data,3) 
                [sersp(trig,c,trial,:,:),itc,powbase,times,freqs,~,~] = newtimef(ep.icaact(comps{sub}(c),:,trial),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                    'plotersp','off','plotitc','off','baseline',NaN,'timesout',100,'nfreqs',50,'freqs',[1,100],'winsize',64) ; 
            end
        end
        %}
        
    end      
    %bersp = sersp - repmat(mean(sersp(:,:,:,:,times2<0),5),[1,1,1,1,100]) ; 
    %for i=1:3 ; figure ; for j=1:32 ; subplottight(4,8,j) ; imagesc(squeeze(mean(sersp(i,:,j,:,:),2))) ; set(gca,'XTick',[],'YTick',[]) ; end ; end
    %figure ; for j=1:64 ; subplot(5,13,j) ; imagesc(squeeze(mean(ersp(1:2,j,:,:),1)),[-6,6]) ; title(j) ; end 
    %tp(EEG) ; 
    %figure,
    %plot(squeeze(mean(mean(sersp(1,:,:,:,times2>0 & times2<5),2),5))','r') ; hold on ; 
    %plot(squeeze(mean(mean(sersp(3,:,:,:,times2>0 & times2<5),2),5))','k') ; 
    %allsersp(sub,:,:,:,:,:) = sersp ; 
    %subeegs{sub} = EEG ; 
    tp(EEG) ;
end

%{
for i=1:6 ; subplot(2,3,i) ;
    plot(squeeze(mean(mean(allsersp(i,1,:,:,:,times2>0 & times2<5),3),6))','r') ; hold on ;
    plot(squeeze(mean(mean(allsersp(i,3,:,:,:,times2>0 & times2<5),3),6))','b') ;
    plot(squeeze(mean(mean(mean(allsersp(i,3,:,:,:,times2<0),3),6),4)),'k','LineWidth',3) ;
    set(gca,'XTick',1:5:length(freqs2),'XTickLabel',round(freqs2(1:5:end))) ; 
    xlabel('frequency(hz)') ; ylabel('log power') ; 
    title(['subject',num2str(i)]) ; 
end
%}














