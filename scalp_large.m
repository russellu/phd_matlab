clear all ; close all  ;
cd c:/shared/resmerged ; 
subs=dir('*') ; subs(1:2) = [] ; 


comps = {[5,7,14,15,19,20,26,31,42],[3,5,7,8,10,11,13,16,23,26,31,38],[3,4,5,7,8,9,11,12,18],[6,7,11,12,13,14,16,18,19,26,29],[16,13,9],[3,4,6,17,22,23,26,35,42,43,44],[6,8,10,14,23,27,30],[2,4,5,7,13,16,17,18,23,28],[7,10,12,17,18,20,27],[10,11,13,21,22,23,24,32],[3,4,6,7,8,10,11,13,15,16,20,22,26],[7,9,11,15,16,31],[6,7,8,9,12,14,17,26,28,39],[1,4,6,8,9,12,13,16,19,20,31,33,40],[2,7,10,20,22,27,28],[7,8,9,10,12,13,15,22,27,29,32],[2,3,7,8,12,13,16,42],[12,13,16,17,19,21,22,31,41,45],[5,6,8,9,10,13,14,15,16,18,21],[3,4,5,6,7,13,21,27,28,48],[2,4,6,7,9,12,16,23,28,32,33,40,41],[3,5,8,15,19,11],[7,14,15,16,17,18,19,22,23,24,27,29,32],[1,7,8,14,16],[4,5,6,7,10,11,12,13,15,17,20,21,30,32,33,37],[2,3,5,6,12,13,14,18,22,25,28,30,35],[1,3,6,7,11,14,16,17,21],[9,14,15,22,28],[1,5,6,7,8,9,10,11,16,17,18],[1,2,36,7,8,11,13,14,15,16,20,23,24,27]} ; 


comps2 = {[20,7,26,15],[13,3,16,23],[12,7,3,11,9],[13,19,6,14],[13,11],[26,22,17,8,42],[14,27,23],[17,16,7,4],[20,7,12],[13,10,11,21],[7,16],[16,11,9,7],...
    [8,12,28,39,50,51],[8,16,19,33,40],[20,22,27],[13,7,12,27],[7,3,2,42],[22,13,12],[14,5,15,17,21],[4,6,7,27],[2,6,7,12,23],[19,15,8,3],[15,7,17,8],...
    [7,8,14],[4,15,32,37],[6,5,12,28],[7,1,17,3],[22,14,28],[10,9,16,18,54],[7,3,2]} ; 

for s=1:length(subs) ; 
    cd(['c:/shared/resmerged/',subs(s).name]) ;  
    ls 
    merged = pop_loadset('merged.set') ; 
    ica = load('saveica.mat') ; ica = ica.saveica ; 
    weights = ica{1} ; sphere = ica{2} ;
    winv = pinv(weights*sphere) ; 
    %{
    acts = weights*sphere*merged.data ; 
    newmerged = merged ; newmerged.data = acts ; 
    epi = pop_epoch(newmerged,{'S 11','S 14'},[-.8,2.8]) ; 
    [sv,si]  = sort((zscore(squeeze(max(std(epi.data(:,:,:),0,2),[],1)))),'ascend') ; 

    for i=1:64 ; disp(i) ; 
            [dersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.data(i,:,si(1:200))),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
    end

    gcf = figure ; for i=1:64 ; subplottight(10,14,i*2-1) ; topoplot(winv(:,i),merged.chanlocs) ;text(-0.55,0.4,num2str(i)) ; end ; 
%   for i=1:64 ; subplottight(10,14,i+70) ; plot(sp(i,:)) ; text(5,5,num2str(i)) ; xlim([0,129]) ; set(gca,'XTick',[],'YTick',[]) ; end
    for i=1:64 ; subplottight(10,14,i*2) ;imagesc(squeeze(dersp(i,:,:)),[-4,4]) ; axis xy ; set(gca,'XTick',[],'YTick',[]) ; end
    set(gcf, 'Position', get(0,'Screensize')); suptitle(subs(s).name) ; 
   %}

   
    
    goodcs = comps{s} ; 
    bads = zeros(1,64) ; bads(goodcs) = 1 ; bads = find(bads==0) ; 
    acts = weights*sphere*merged.data ; acts(bads,:) = 0 ; 
    invdat = winv*acts ; 
    
    newmerged = merged ; newmerged.data = invdat ; 
    epi = pop_epoch(newmerged,{'S 11','S 13','S 15','S 14'},[0.5,1.5]) ;  %,'S 12','S 13','S 14','S 15','S 16'
    basepi = pop_epoch(newmerged,{'S 11','S 13','S 15','S 14'},[-1,0]) ; 
    [sv,si]  = sort((zscore(squeeze(max(std(epi.data(:,:,:),0,2),[],1)))),'ascend') ; 

    freqs = (abs(fft(squeeze(epi.data),[],2))) ; 
    basefreqs = (abs(fft(squeeze(basepi.data),[],2))) ; 
    newfreqs = zeros(size(basefreqs)) ; for i=1:64 ; newfreqs(i,:,:) = imresize(squeeze(freqs(i,:,:)),[size(basefreqs,2),size(basefreqs,3)]) ; end
    a = squeeze(mean(freqs-basefreqs,3)) ;
    %figure,imagesc(a(:,1:128),[-.5,.5]) ; title(subs(s).name) ;
    allas(s,:,:) = a(:,1:128) ; 
    postersp = a ; save('postersp','postersp') ; 
  
    epi = pop_epoch(newmerged,{'S 11','S 13','S 15','S 14','S 12','S 16'},[-.8,2.8]) ;  %,'S 12','S 13','S 14','S 15','S 16'
    clear dersp
    for i=1:64 ; disp(i) ; 
      %  for j=1:135 ; 
            [dersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epi.data(i,:,si(1:450))),epi.pnts,[epi.xmin,epi.xmax],epi.srate,0,...
                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',0,'verbose','off') ; 
      %  end
    end
    elabs = {merged.chanlocs.labels} ; 
    gcf = figure ; for i=1:64 ; subplot(5,13,i) ; imagesc(squeeze(dersp(i,:,:)),[-3,3]) ; title(elabs{i}) ; end; suptitle(subs(s).name) 
    set(gcf, 'Position', get(0,'Screensize')); 
    save('dersp','dersp') ;
    
    
    
end