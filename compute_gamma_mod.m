clear all ; close all
cd C:\shared\allres
subs = {'alex','alexandra3','audrey','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','leila','lisa','marc',...
    'marie','mathieu','maxime','menglu','mingham','patricia','po','russell','suhan2','sunachakan','tah','tegan2','vincent'} ; 
trigs = {'S 11','S 12','S 13','S 14','S 15','S 16'} ; 

goodcs = {[5,8,9,15,19,20,25,28,29],[4,5,8,10,15,16,21,22,24,29,31,39,40],[2,3,4,5,6,7,8,10,14,15,16,21,44],[3,7,9,10,16,23,26,30],[14,15,17,25],...
    [3,5,6,15,18,19,25,29],[4,5,7,8,12,14,18,21,28],[3,5,6,7,10,13,14,23,40],[3,4,7,8,10,16,17,26,38],[6,10,12,13,14,16,22,42],[4,9,17,18,20,28,30,35],...
    [6,7,8,14,15,16,21,28,30,32],[4,5,8,9,12,18,23,26,35,38],[2,5,6,7,11,12,15,22,28,34,36],[2,3,4,6,7,10,19],[4,6,7,13,14,18,20,21,28,31,38],...
    [2,4,5,6,7,8,9,10,11,12,14,19,32],[4,6,8,10,11,13,20,23,31,35,37,48],[1,2,5,8,9,14,16,19,20,21,22,25,32],[2,7,8,11,26],[4,5,6,7,8,9,11,12,13,17,18,19,23,37],...
    [3,5,7,8,9,12,21,28],[5,6,7,9,10,13,15,16,21,23,25,28,31,32,33,41],[4,5,6,8,9,10,19,21,25,26],...
    [3,5,6,7,8,9,10,14,15,18,20,27,31,34],[3,4,5,6,7,8,10,11,12,13,15,16,27,29],[1,3,6,7,8,9,23,24,25,27,30],[9,10,21,25,29,31,33],...
    [4,5,7,8,9,11,14,15,17,21,29,42],[2,3,4,5,7,8,9,10,11,13,14,15,18,26,36]} ;

postgoodcs = {[8,9,25],[4,5,22,24],[14,15,21],[7,9,16],[14,15,25],[15,25,29],[8,12,21],[5,7,23,40],[3,7,8,17],[6,10,13,16],[4,35],[8,15,37],[8,12,18,26],...
    [2,6,15,22],[10,19],[4,13,20],[2,4,8],[4,6,10],[2,20,22],[2,7,11],[5,9,11],[5,9,12],[13,15,16,23],[6,8,9],[3,10,20,31,34],[6,7,12,29],[1,6,7],[10,25,31],...
    [8,14,17,21,42],[3,4,11]} ; 

postgoodcs2 = {[5,8,9,15,25,28],[4,5,8,10,15,22,24],[3,5,14,15,16,21],[7,9,16,23],[14,15,25],[15,18,19,25,29,31],[8,12,21],[5,7,23,40],[3,7,8,17],[6,10,13,16,42,48],[4,22,30,35],[7,8,15,32,37],[8,12,18,23,26,38],[2,6,7,11,15,22,36,43,53],[2,3,4,6,10,19],[4,13,20,21,28,31,38],[2,4,8,9,14,32],[4,6,8,10,23,35,48],[2,22,25],[2,7,8,11,26],[4,5,7,8,9,11],[5,7,8,9,12],[13,15,16,23,25,32],[6,8,9,10],[3,10,20,31,34],[6,7,10,11,12,15,29],[1,6,7,9,23],[9,10,25,31],[4,5,8,14,17,21,24,42],[3,4,11,18]} ;

badts = {...
    {[68,91,114],[29,44,92,131],[31,94,127,128],[64,128],[20,90],[21,75]},
    {[],[75],[],[56,127],[8,14,32],[111]},
    {[73,111],[100],[85],[31,116,118],[78,99],[27]},
    {[2],[2,80],[40,64,90],[54],[],[]},
    {[39],[],[7,126],[32],[73,113],[74,77]},
    {[],[11],[],[],[10],[120]},
    {[21],[],[128],[],[42,74],[42,135]},
    {[55,111],[52,75,100],[65,117],[],[47,75],[]},
    {[94],[106],[],[12,50,62,102],[],[]},
    {[47,94,112,117,128],[],[35,48,103,110,132],[26,61,135],[37,38,41],[14,22,106,111,129]},
    {[],[130],[80,117],[6,122],[],[]},
    {[88],[80,106],[],[26,59],[],[99]},
    {[34,90],[6,57,67,129],[21,57],[82,98],[35],[20]},
    {[],[],[],[44,98],[],[]},
    {[79,96],[],[],[58,94],[],[25]},
    {[58,119],[58],[18],[42],[38,57,60],[29,57]},
    {[9,13,47,100,104,122],[63,68],[15,19,44,46,47,54],[7,26,30,104],[21,62,66,69,125],[5,51,63]},
    {[41,50,54],[39,71],[66,73],[],[37],[17,20,53]},
    {[],[18,22],[10,22,110],[69],[41,103,115],[96]},
    {[11,127],[14,46,65,93],[98],[119,133],[28,57,83],[55,117]},
    {[],[],[],[],[],[11,93]},
    {[17,21,52,83],[54],[88],[38],[41],[129]},
    {[5,15,110],[50],[85,107,109],[31,77,106,107],[],[5,33]},
    {[],[84],[],[],[34],[]},
    {[],[],[],[],[],[9,41]},
    {[117,109],[1,16,23],[97],[33],[],[]},
    {[10,54,106],[3,24,106],[60],[107],[9,22,31,108],[]},
    {[79,80,89],[83],[],[17],[],[63]},
    {[56],[44,135],[76],[],[121],[]},
    {[35,44,108],[17,23,62,84],[14,98,130],[],[37],[14,27]}} ;

for sb=1:length(subs)
    cd(['c:/shared/allres/',subs{sb}]) ; 
    ls 
    fullica = load('fullica.mat') ; fullica = fullica.fullica ; 
    weights = fullica{1} ; sphere = fullica{2} ; winv = pinv(weights*sphere) ; 
    times = load('times') ; times = times.times;  
    freqs = load('freqs') ; freqs = freqs.freqs ; 
    eeg = pop_loadset('cleanfilt.set') ; 
    %ersp = load('ersp') ; ersp = ersp.ersp ; 
    bersp = load('bersp') ; bersp = bersp.bersp ;
    %figure,for i=1:6 ; subplot(2,3,i) ; imagesc(squeeze(mean(bersp(i,postgoodcs{sb},:,:),2)),[-4,4]) ; axis xy ; end ; suptitle(subs{sb}) ; 
    allbersp(sb,:,:,:) = squeeze(mean(bersp(:,postgoodcs2{sb},:,:),2)) ; 
    neweeg = eeg ; neweeg.data = weights*sphere*eeg.data ; 
 
    %{
    clear ersp ; 
    for tr=1:length(trigs)
    epica = pop_epoch(neweeg,{trigs{tr}},[-.6,2.5]) ;
    ts = zeros(1,135) ; ts(badts{sb}{tr}) = 1 ; goods = find(ts==0) ; 
    for i=1:length(postgoodcs{sb}) 
            [ersp(tr,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(postgoodcs{sb}(i),:,goods)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',50,'baseline',0,'verbose','off','timesout',200) ; 
    end
    end
    figure,for i=1:6 ; subplot(2,3,i) ; imagesc(squeeze(mean(ersp(i,:,:,:),2)),[-3,3]) ; axis xy ; end ; suptitle(subs{sb}) ; 
    %}
    %allersp(sb,:,:,:) = squeeze(mean(ersp,2)) ; 
        
       %{ 
    clear ersp ; 
    for tr=1:length(trigs)
    epica = pop_epoch(neweeg,{trigs{tr}},[-1,3]) ; 
    for i=1:length(postgoodcs{sb}) 
        for j=1:135
            [ersp(tr,i,j,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(epica.data(postgoodcs{sb}(i),:,j)),epica.pnts,[epica.xmin,epica.xmax],epica.srate,0,...
                                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',64,'baseline',NaN,'verbose','off','timesout',100) ; 
        end
    end
    end     
    bersp = ersp - repmat(mean(ersp(:,:,:,:,times<0),5),[1,1,1,1,100]) ; 
    mbersp = squeeze(mean(bersp,2)) ; 
    %for i=1:6 ; subplot(2,3,i) ; imagesc(squeeze(mean(mbersp(i,:,:,:),2)),[-4,4]) ;axis xy ; end
    for i=1:6 ; figure('units','normalized','outerposition',[0 0 1 1]) ; 
        for j=1:135 ; subplottight(10,14,j) ; imagesc(squeeze(mbersp(i,j,:,:)),[-18,18]) ; 
            axis xy ; set(gca,'XTick',[],'YTick',[]) ; text(5,50,num2str(j)) ; 
        end ; 
    end 
    %}
    si = load('si') ; si = si.si ;
    acts = weights*sphere*eeg.data ; 
    bads = zeros(1,64) ; bads(postgoodcs2{sb}) = 1 ; bads = find(bads==0) ; 
    acts(bads,:) = 0 ; 
    cleandat = winv*acts ; neweeg = eeg ; neweeg.data = cleandat ; 
    for tr=1:length(trigs)
        ep = pop_epoch(neweeg,{trigs{tr}},[-1,3]) ; clear ersp ; 
        fbase = abs(fft(ep.data(:,1:256,:),[],2)) ; 
        ftask1 = abs(fft(ep.data(:,257:257+255,:),[],2)) ; 
        fdiff1 = (ftask1)-(fbase) ; 
        ftask2 = abs(fft(ep.data(:,257+255:257+255+255,:),[],2)) ; 
        fdiff2 = (ftask2)-(fbase) ; 
        meandiff = (fdiff1+fdiff2) / 2 ; ts = zeros(1,135) ; ts(badts{sb}{tr}) = 1 ; goods = find(ts==0) ; 
        alldiffs(tr,:,:) = squeeze(mean(meandiff(:,1:128,si{tr}(1:120)),3)) ; 
        alltasks(tr,:,:) = squeeze(mean(ftask1/2+ftask2/2,3)) ;
        allbase(tr,:,:) = squeeze(mean(fbase,3)) ; 
        %subplot(2,3,tr) ; topoplot(squeeze(mean(mean(meandiff(:,60:80,:),2),3)),eeg.chanlocs,'maplimits',[-5,5])
    end
    figure,
    subplot(1,2,1) ; topoplot(squeeze(mean(mean(alldiffs([1,3,4,5],:,50:90),1),3)),eeg.chanlocs,'maplimits',[-3,3]) ; title(subs{sb}) ; 
    subplot(1,2,2) ; topoplot(squeeze(mean(mean(alldiffs([1,3,4,5],:,10:30),1),3)),eeg.chanlocs,'maplimits',[-50,50]) ; title(subs{sb}) ; 
    meangamma = squeeze(mean(mean(alldiffs([1,3,5],:,50:80),1),3)) ; 
    meanbeta = squeeze(mean(mean(alldiffs([1,3,5],:,10:30),1),3)) ; 
    allgamma(sb,:) = meangamma ; allbeta(sb,:) = meanbeta ;  
    save('alldiffs','alldiffs') ; save('alltasks','alltasks') ; save('allbase','allbase') ; 
    %{
    
    si = load('si') ; si = si.si ; 
    acts = weights*sphere*eeg.data ; 
    bads = zeros(1,64) ; bads(postgoodcs{sb}) = 1 ; bads = find(bads==0) ; 
    acts(bads,:) = 0 ; 
    cleandat = winv*acts ; neweeg = eeg ; neweeg.data = cleandat ; 
    cfreqs = 1:2:120 ; 
    for tr=1:length(trigs)
        allfreqs = zeros(125,64,60,640) ; 
        ep = pop_epoch(neweeg,{trigs{tr}},[-.5,2]) ; clear ersp ; 
        for i=1:64 
            [ersp(tr,i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ep.data(i,:,:)),ep.pnts,[ep.xmin,ep.xmax],ep.srate,0,...
                                'plotersp','off','plotitc','off','freqs',[1,120],'nfreqs',60,'winsize',50,'baseline',NaN,'verbose','off','timesout',200) ; 
        end
        ep.data = ep.data(:,:,1:135) ; 
        %csi = si{tr} ; %[~,csitr] = sort(csi,'ascend') ; 
        ts = zeros(1,135) ; ts(badts{sb}{tr}) = 1 ; goods = find(ts==0) ; 
        for elec=1:64
            for freq=1:length(cfreqs)
                allfreqs(:,elec,freq,:) = eegfiltfft(squeeze(ep.data(elec,:,goods(1:125)))',ep.srate,cfreqs(freq)-2,cfreqs(freq)+2) ; 
            end
        end
        allfreqs = abs(allfreqs) ; 
        basefreqs = allfreqs - repmat(mean(allfreqs(:,:,:,ep.times<0),4),[1,1,1,size(allfreqs,4)]) ; 
        stimfreqs(tr,:,:,:) = mean(basefreqs,1) ; 
    end
    figure,for i=1:6 ; subplot(2,4,i) ; imagesc(squeeze(mean(stimfreqs(i,:,:,ep.times>0 & ep.times<2000),4))) ; end ; suptitle(subs{sb}) ; 
    subplot(2,4,7) ; 
    topoplot(squeeze(mean(mean(mean(stimfreqs([1,5],:,32:40,ep.times>0 & ep.times<2000),3),4),1)),eeg.chanlocs,'maplimits',[-.008,.008])
    save('stimfreqs','stimfreqs') ; eptimes = ep.times ; save('eptimes','eptimes') ; 
    %}
end

