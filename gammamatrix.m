clear all ; close all 
%whitebg('black')

trigcount = 1 ;
  for i=10:99
      if i<100
     	rtrigs{trigcount} = ['S',' ',num2str(i)] ;
      else 
          rtrigs{trigcount} = ['S',num2str(i)] ;
      end
      
      trigcount = trigcount + 1 ; 
  end
subjects = { 
    'russ9'
} ;

subjcount = 1 ; 
nparams = 90 ;
for subj=1:size(subjects,1) 
    cd(strcat('c:/raw_eeg/',subjects{subj})) ;
    disp(strcat('processing subject = ' , subjects{subj})) ;  
    fnames = dir('removed*vis*.set') ;     
    stimlats = cell(1,nparams) ; 
    epochEEGs = cell(size(fnames,1),nparams) ;
    allepochs = cell(1,nparams) ; 
    
    for fn=1:size(fnames,1) % for all set files that were loaded
        disp(strcat('processing subject = ' , subjects{subj})) ;  
        fname=fnames(fn).name ; disp(['current set = ',fname]) ;
        EEG = pop_loadset(fname); EEG = eeg_checkset(EEG);
        disp(['notch filtering 59.95:60.05Hz']) ; 
        f60 = eegfiltfft(EEG.data,EEG.srate,59.95,60.05) ; 
        EEG.data = EEG.data - f60 ; 
        ch=({EEG.chanlocs.labels});
        trigs = {EEG.urevent.type};
        latencies = cell2mat({EEG.urevent.latency});      
        for i=1:nparams
            stimlats{i} = latencies(ismember(trigs,rtrigs{i})) ;
            epochEEGs{fn,i} = pop_epoch(EEG,rtrigs(i),[-.4,1.25]) ; 
        end  
    end    
    for i=1:size(epochEEGs,2) % for all stimulus parameters 
        allepochs{i} = epochEEGs{1,i} ; 
        for j=2:size(epochEEGs,1)
            allepochs{i} = pop_mergeset(allepochs{i},epochEEGs{j,i}) ;
        end
    end
    for i=1:size(allepochs,2) % save all the merged trials as their own .set files
        pop_saveset(allepochs{i},'filename',[stimnames{i},'.set']) ;  
    end
    
    subjepochs{subj} = allepochs ; 
    
    %goods = load('goods.mat') ; 
    %goodmats{subj} = goods.goods ; 
    
end

% get the newtimef for each condition
for subj = 1:size(subjepochs,2)
    clear ersp
    for s=1:size(subjepochs{subj},2) 
        disp(['calculating ERSP, stim param = ',num2str(s), ' subject = ',subjects{subj}]) ; 
        for i=1:size(EEG.data,1)          
            [ersp(s,i,:,:),itc,powbase,times,freqs,erspboot,itcboot] = newtimef(squeeze(subjepochs{subj}{s}.data(i,:,:)),subjepochs{subj}{s}.pnts, [-400,1250], subjepochs{subj}{s}.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,150],'nfreqs',150,'verbose','off','baseline',0,'winsize',floor(subjepochs{subj}{s}.srate./2),'timesout',200);
        end
    end
    subjersp{subj} = ersp ; 
end

allpost = [22,51,23,52,24,55,56,57,26,27,28,54,58] ;

for i=1:90 ;
    subplot(9,10,i) ; 
    imagesc(squeeze(mean(ersp(i,allpost,:,:),2)),[-3,3]) ; 
    
end

timess = 27:200 ; 
elecs = [28,58,57,56,52,54,23] ;

clear gammamat
for f = 1:150-5
for i=1:9
    for j=1:10
     %   gammamat(i,j) = mean(mean(mean(ersp((i*10-10)+j,elecs,65:75,30:150),2),3),4) ; 
    %    betamat(i,j) = mean(mean(mean(ersp((i*10-10)+j,elecs,7:25,30:150),2),3),4) ; 
    allms(f,i,j) = mean(mean(mean(ersp((i*10-10)+j,elecs,f:f+5,30:150),2),3),4) ; 
    end
end
end
subplot(1,2,1) ; imagesc(imfilter(gammamat,fspecial('gaussian',3,1))) ;
subplot(1,2,2) ; imagesc(imfilter(betamat,fspecial('gaussian',3,1))) ; 

icount = 1 ;
for i=1:3:90 ; 
    subplot(4,8,icount) ;
    imagesc(squeeze(allms(i,:,:)),[-2,2]) ; 
    title(['hz = ',num2str(i),':',num2str(i+5)]) ; 
    icount = icount + 1 ;
end
suptitle('increasing x = increasing spatial frequency, increasing y = increasing temporal frequency');












