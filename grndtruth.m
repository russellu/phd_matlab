%cd('c:/Vision/Raw Files') ; ls 
cd('c:/shared/denoising_MR') ; 
clear all ; close all ; 
EEGdenoise = pop_loadbv('.','Russell_test_2015-08-05_EEGoutsideMRIroom.vhdr') ; 
EEGnoise = pop_loadbv('.','Russell_test_2015-08-05_EEGnoEPI.vhdr') ; 
EEGdenoise = pop_chanedit(EEGdenoise,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;
EEGnoise = pop_chanedit(EEGnoise,'lookup','C:\eeglab10_0_0_0b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp') ;

% get the power spectrum for both channels:
denoise500 = pop_resample(EEGdenoise,500) ; 
noise500 = pop_resample(EEGnoise,500) ; 
[td,fd] = spectopo(denoise500.data(46,:),0,denoise500.srate,'plot','off') ; 
[tn,fn] = spectopo(noise500.data(46,:),0,noise500.srate,'plot','off') ; 
%plot(td,'LineWidth',2) ; hold on ; plot(tn,'r','LineWidth',2) ; plot(tc,'k') ;set(gca,'XTick',(1:50:length(fn)),'XTickLabel',fn(1:50:length(fn))) ; legend({'oustide bore','inside bore','denoised inside bore'}) ; 

eeg2 = noise500 ;
%eeg2.data = eegfiltfft(eeg2.data,eeg2.srate,1,250) ; 
% get the peaks within a long sliding window(2 seconds) for the correlation
% function
clear maxinds hrepochs hrepochs2
peakdata = eeg2.data(32,:) ;  
winlength1 = eeg2.srate*2 ; maxinds = zeros(1,size(eeg2.data,2)) ; 
wprev = 25 ; % floor(eeg2.srate*.33) ; 
wnext = floor(eeg2.srate) ; 
for w=1:winlength1:size(eeg2.data,2)-winlength1
    windowi = peakdata(w:w+winlength1) ; 
    maxi = find(windowi==max(windowi)) ; 
    maxinds(w+maxi) = 1 ; 
end
peakinds = find(maxinds==1) ; 
for i=2:length(peakinds)-1
    hrepochs(i,:) = peakdata(peakinds(i)-wprev:peakinds(i)+wnext) ; 
end
mhr = squeeze(mean(hrepochs,1)) ; 

%%% get the second peaks (more numerous) and correlate
winlength2 = eeg2.srate ; maxinds = zeros(1,size(eeg2.data,2)) ; 
for w=1:round(winlength2):size(eeg2.data,2)-winlength2
    windowi = peakdata(w:w+winlength2) ; 
    maxi = find(windowi==max(windowi)) ; 
    maxinds(w+maxi) = 1 ; 
end
peakinds = find(maxinds==1) ; 
for i=1:length(peakinds)-1
    if peakinds(i)-wprev > 0 & peakinds(i) + wnext < length(peakdata)
        hrepochs2(i,:) = peakdata(peakinds(i)-wprev:peakinds(i)+wnext) ; 
    end
end
clear corrs ; 
for i=1:size(hrepochs2,1) ; corrs(i) = corr2(hrepochs2(i,:),mhr) ; end ; 
goods = corrs>0.65 ; goodepochs = hrepochs2(goods,:) ; 
goodpeakinds = peakinds(goods) ; 
% find the double peak indices and remove the lower one
zdiffs = zscore(diff(goodpeakinds)) ; 
tooclose = find(zdiffs<-2) ; 
for i=1:length(tooclose)
    if peakdata(goodpeakinds(tooclose(i))) < peakdata(goodpeakinds(tooclose(i))-1)       
        tooclose(i) = (tooclose(i)+1) ;  
    end
end
goodepochs(tooclose,:) = [] ; goodpeakinds(tooclose) = [] ; 
geps = goodepochs ; 
plot(eeg2.data(32,:)) ; vline(goodpeakinds)
imagesc(geps);
%for i=1:5 ; figure ; plot(geps(i,:)) ; hold on ; plot(mean(geps),'r') ; hold on ;  plot(mean(geps)-geps(i,:),'k') ; end
%corrs = corr(geps') ; 
%for i=1:size(corrs,1) ; [sortcorrs(i,:),sortcorrinds(i,:)] = sort(corrs(i,:),'descend') ; end 

%for elec = 1:64 ;
elec = 46 ; 
visdata = eeg2.data(elec,:) ; clear visepochs
templatedata = zeros(size(visdata)) ;
dist = 10 ; 
for i=1:length(goodpeakinds)
    visepochs(i,:) = visdata((goodpeakinds(i))-wprev:(goodpeakinds(i))+wnext) ; 
end
mgeps = mean(visepochs) ; 
meandiff = diff(mgeps) ; zcross = zeros(1,length(meandiff)) ; 
for i=1:length(meandiff)-1 ; if (meandiff(i) > 0 & meandiff(i+1)<0) | (meandiff(i) < 0 & meandiff(i+1) > 0 ); zcross(i) = 1 ; end ; end
%plot(meandiff) ; vline(find(zcross==1)) ;hline(0) ; hold on ; plot(mgeps) ; 
peakinds = find((zcross==1))+1 ; 
zdiffs = zscore(diff(peakinds)) ; peakinds(zdiffs<0)= [] ; 

%get the custom template for each heart beat
correpochs = corr(visepochs') ; 
for i=1:size(correpochs,1) ; 
    [sv,si] = sort(correpochs(i,:),'descend') ; 
    corrtemplates(i,:) = squeeze(mean(visepochs(si(2:400),:))) ; 
end

meancorrected = mgeps-smooth(mgeps,50)' ; 
%[sv,si] = sort(abs(meancorrected(peakinds)),'descend') ; 
zthresh = 0 ; 
peakinds2 = peakinds((zscore(abs(meancorrected(peakinds)))>zthresh)) ; % thresholded peaks
zinds = find(zscore(abs(meancorrected(peakinds)))>zthresh) ; zvalues = zscore(meancorrected(peakinds)) ; zvalues = zvalues(zinds) ; 
zinds = peakinds(zinds)  ;
%peakinds2 = 100:50:size(visepochs,2) ; 
%zinds = peakinds2 ; 
clear fullinds pfullinds ;
divconst = 2; dc = 14  ;
for pk=1:length(zinds) % get the surrounding window length for each peak (+/- distance to next/prev peak). 
   if pk==1
      fullinds{pk} = round(1:zinds(pk) + (zinds(pk+1)-zinds(pk))/divconst) ; 
      pfullinds{pk} = zinds(pk)-dc:zinds(pk)+dc ; 
   elseif pk==length(zinds)
       fullinds{pk} = round(zinds(pk)-(zinds(pk)-zinds(pk-1))/divconst : size(geps,2)) ; 
       pfullinds{pk} = round(zinds(pk)-dc : zinds(pk)+dc)  ; 
   else
      fullinds{pk} = round(zinds(pk)-(zinds(pk)-zinds(pk-1))/divconst : zinds(pk) + (zinds(pk+1)-zinds(pk))/divconst) ; 
      pfullinds{pk} = round(zinds(pk)-dc : zinds(pk) + dc) ; 
   end
end
%fullinds = fullinds(zinds) ; 
fullinds2 = fullinds ; 
%winds2 = winds((zscore(abs(meancorrected(peakinds)))>zthresh),:) ; % windows for thresholded peaks
%%% find a suitable scale and shift to align the peaks.

%for i=1:1 ; figure ; plot(smooth(visepochs(i,:),10)) ; vline(peakinds2) ; hold on ; plot(mgeps,'k') ;  plot(visepochs(i,:),'m') ;end
%%% get more epochs (split up the individual epochs into smaller sub-epochs)
%for i=1:size(visepochs,1) ; plot(visepochs(i,fullinds2{3})) ; hold on ; end
% fit the curve to the mean (Can use a more specialized mean to take into
% account subject motions)
clear allcurves allsegs sdiffs
for x=1:size(visepochs,1) ; 
    disp(x) ; %figure,
    for inds=1:length(pfullinds)
        curve1 = visepochs(x,pfullinds{inds}) ; 
        curve2 = corrtemplates(x,pfullinds{inds}) ; 
        shiftrange = -8:8 ; 
        transrange = -80:5:80 ; 
        multrange = [-1.2,-1,-.9,-.8,.8,.9,1,1.1,1.2] ; 
        clear allcurves sqrdiffs
        for i=1:length(shiftrange) % for all parameters ijk
            for j=1:length(transrange)
                for k=1:length(multrange)
                    if pfullinds{inds}+shiftrange(i) < length(mgeps) & min(pfullinds{inds}+shiftrange(i))>0
                        curve2 = corrtemplates(x,pfullinds{inds}+shiftrange(i))*multrange(k)+transrange(j) ;
                        sqrdiffs(i,j,k) = sum((curve1-curve2).^2) ;  
                        allcurves(i,j,k,:) = curve2 ; 
                    end
                end
            end
        end
        sqrdiffs(sqrdiffs==0) = NaN ; 
        indy = find(sqrdiffs==min(min(min(sqrdiffs)))) ; 
        [ix,iy,iz] = ind2sub(size(sqrdiffs),indy) ; 
        allsegs{x,inds} = squeeze(allcurves(ix,iy,iz,:)) ; 
        constants(x,inds,:) = [shiftrange(ix),transrange(iy),multrange(iz)] ; 
        %subplot(2,2,inds) ; plot(allsegs{x,inds}) ; hold on ; plot(visepochs(x,pfullinds{inds}),'r') ; 
    end
end


% put the constants into the indices (make a full time series)
clear allveci
for i=1:size(constants,1)
    %transveci = zeros(length(fullinds2),length(mgeps),3) ; 
    for j=1:length(fullinds2)
        %transveci(j,fullinds2{j},:) = repmat(squeeze(constants(i,j,:))',[length(fullinds2{j}),1]) ; 
        allveci(i,fullinds2{j},:) = repmat(squeeze(constants(i,j,:))',[length(fullinds2{j}),1]) ;
    end
end

% transform the average template, first blend the transformation:
clear smoothtrans ptemplates
for i=1:size(allveci,1)
    for j=1:size(allveci,3)
        smoothtrans(i,:,j) = smooth(squeeze(allveci(i,:,j)),10) ;
    end
end
ptemplates = zeros(size(smoothtrans,1),size(visepochs,2)) ; 
for i=1:size(smoothtrans,1)
    for j=1:size(smoothtrans,2)
        xij = round(smoothtrans(i,j,1)) ; yij = smoothtrans(i,j,2) ; zij = smoothtrans(i,j,3) ; 
        if min([j-xij,j+xij]) > 0 & j+xij <= length(mgeps)
            ptemplates(i,j) = corrtemplates(i,j+xij)*zij + yij ; 
        else ptemplates(i,j) = corrtemplates(i,j)*zij + yij ; 
        end
    end  
end
%{
% how to properly blend the template chunks?
clear templates
for i=1:size(allsegs,1) ;
    templates(i,:) = zeros(1,length(mgeps)) ; 
    for j=1:size(allsegs,2)
        templates(i,fullinds2{j}) = allsegs{i,j} ;
       % if j>1 % find the overlapping indices and average them
       %     bothinds = intersect(fullinds2{j},fullinds2{j-1}) ;
       %     indsi = 1:length(bothinds) ; 
       %     indsi1 = length(allsegs{i,j-1})-length(bothinds)+1 : length(allsegs{i,j-1}) ; 
       %     templates(i,bothinds) = (allsegs{i,j}(indsi) + allsegs{i,j-1}(indsi1)) ./ 2 ;
       % end
    end
end
%%% find all the discontinuities and display them:
for i=1:length(fullinds2) ; spikes(i) = fullinds2{i}(1) ; end
templates2 = templates ; 
for i=1:size(templates,1)
    %figure,
    %smoothtemplates(i,fullinds{1}) = 
   for j=1:length(spikes)
       if spikes(j)-20 > 0 & spikes(j) + 20 < length(mgeps) ; 
       cspike = templates(i,spikes(j)-20:spikes(j)+20) ; 
       smoothspike = smooth(cspike,5) ; 
     %  subplot(3,5,j) ; hold on ; 
     % plot(templates(i,spikes(j)-20:spikes(j)+20)) ;  
     % plot(visepochs(i,spikes(j)-20:spikes(j)+20),'r') ;  
     %  plot(smoothspike,'g') ; 
       templates2(i,spikes(j)-5:spikes(j)+5) = smoothspike(round(length(smoothspike)/2-5):round(length(smoothspike)/2+5)) ; 
       end
   end
end
%}

for i=1:size(ptemplates,1)
   figure, 
   subplot(1,2,1) ; hold on ; plot(visepochs(i,:)) ; plot(mgeps,'c') ; vline(peakinds2)
   plot(visepochs(i,:)-ptemplates(i,:),'k') ;
   plot(ptemplates(i,:),'k','LineWidth',2); 
   subplot(1,2,2) ; imagesc(squeeze(smoothtrans(i,:,:))) ; 
end

for i=1:size(visepochs,1) ; 
   catepochs((i-1)*size(visepochs,2)+1:(i-1)*size(visepochs,2)+size(visepochs,2)) = visepochs(i,:) ; 
   subepochs((i-1)*size(visepochs,2)+1:(i-1)*size(visepochs,2)+size(visepochs,2)) = visepochs(i,:)-ptemplates(i,:) ; 

    
end

[t1,f1] = spectopo(catepochs,0,500,'plot','off') ; 
[t2,f2] = spectopo((subepochs),0,500,'plot','off') ; 

plot(t1) ; hold on ; plot(t2,'r') ; plot(td,'k') ; legend({'inside bore','inside bore, denoised','outside bore'}) ; 

visdata2 = eeg2.data(elec,:) ;
for i=1:length(goodpeakinds)
    visdata2((goodpeakinds(i))-wprev:(goodpeakinds(i))+wnext) = visdata2((goodpeakinds(i))-wprev:(goodpeakinds(i))+wnext)-templates2(i,:) ; 
end
alldata(elec,:) = visdata2 ; 
%end
%plot(visdata) ; hold on ;plot(visdata2,'r') ; 
%save('alldata','alldata','-v7.3') ; 