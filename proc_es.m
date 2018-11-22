clear all ; close all ; 
cd('C:\shared\mkt') ; 
fid = fopen('ES.txt') ; 
data = textscan(fid,'%s %s %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
fclose(fid) ; 
d5 = data{6} ; d5 = d5(1:15:end) ; 
d5 = d5(length(d5)/2:end) ; 
%d5 = d5(end-50000:end) ; 

% make a moving average:
clear mvg p mean
for p = 1:500  
    for i=p+1:length(d5)
        if i==p+1
            mvg(i,p) = mean(d5(i-p+1:i)) ; 
        else
            mvg(i,p) = mvg(i-1,p) + d5(i)/(p) - d5(i-p)/(p) ;
        end  
    end
end

clear allprofits allntrades alltracks
m1s = 1:20:500 ; m2s = 1:20:500 ; 
m1count = 1 ; 
for m1=1:10:100 ; m2count = 1 ; disp(m1) ; 
    for m2=1:10:100 ; tpcount = 1 ; 
        for tp=0.5:4:30
            trade = false ; entry = 0 ; buy = false ; sell = false ; profit = 0 ; ntrades = 0 ; alltrades = [] ; sellentries = [] ; buyentries = [] ; sellexits = [] ; buyexits = [] ; 
            for i=201:length(d5)
                if trade == false
                    if d5(i) > mvg(i,m1) % buy mode
                        if d5(i) < mvg(i,m2) && d5(i-1) > mvg(i,m2) % cross down
                            trade = true ; entry = d5(i) ; buy = true ; ntrades = ntrades + 1  ; buyentries(length(buyentries)+1) = i ; 
                        end
                    elseif d5(i) < mvg(i,m1) % sell mode
                        if d5(i) > mvg(i,m2) && d5(i-1) < mvg(i,m2) % cross up
                            trade = true ; entry = d5(i) ; sell = true ; ntrades = ntrades + 1 ; sellentries(length(sellentries)+1) = i ; 
                        end
                    end
                elseif trade == true
                    if buy
                        if d5(i) - entry > tp || d5(i) - entry < -tp
                            buy = false ; profit = profit + d5(i)-entry ; trade = false ; alltrades(length(alltrades)+1) = d5(i) - entry ; buyexits(length(buyexits)+1) = i ; 
                        end
                    elseif sell
                        if entry - d5(i) > tp || entry - d5(i) < -tp
                            sell = false ; profit = profit + entry-d5(i) ; trade = false ; alltrades(length(alltrades)+1) = entry - d5(i) ; sellexits(length(sellexits)+1) = i ; 
                        end                   
                    end
                end    
            end ; allprofits(m1count,m2count,tpcount) = profit ; allntrades(m1count,m2count,tpcount) = ntrades ; alltracks{m1count,m2count,tpcount} = alltrades ;
            allbuyentries{m1count,m2count,tpcount} = buyentries ; allsellentries{m1count,m2count,tpcount} = sellentries ; allbuyexits{m1count,m2count,tpcount} = buyexits ; allsellexits{m1count,m2count,tpcount} = sellexits ; 
            tpcount = tpcount + 1 ; 
        end ; m2count = m2count + 1 ; 
    end ; m1count = m1count + 1 ; 
end

avgn = allprofits./allntrades ; 
figure,for i=1:size(avgn,3) ; subplot(ceil(sqrt(size(avgn,3))),ceil(sqrt(size(avgn,3))),i) ; imagesc(squeeze(avgn(:,:,i)),[.5,3]) ; title(i) ; end
figure,for i=1:size(avgn,3) ; subplot(ceil(sqrt(size(avgn,3))),ceil(sqrt(size(avgn,3))),i) ; imagesc(squeeze(allntrades(:,:,i)),[0,1000]) ; title(i) ; end
for i=1:size(alltracks,1) ; for j=1:size(alltracks,2) ; for k=1:size(alltracks,3) ; if ~isempty(alltracks{i,j,k}) ; corrs(i,j,k) = corr2(cumsum(alltracks{i,j,k}),1:length(alltracks{i,j,k}-1)) ; end ; end ; end ; end

figure,
for i=1:size(avgn,3) ; subplot(ceil(sqrt(size(avgn,3))),ceil(sqrt(size(avgn,3))),i) ; 
    imagesc(squeeze(allntrades(:,:,i)).*(squeeze(avgn(:,:,i))>.5).*(squeeze(avgn(:,:,i))),[1,3000]) ; title(i) ; 
end
figure,plot(((cumsum(alltracks{17,12,13}-.25))),'LineWidth',2) ; %hold on ;plot(mat2gray(d5),'r') ; 

x=19;y=23;z=15;
csum = cumsum(alltracks{x,y,z}-.25) ; diffs = alltracks{x,y,z} ; 
maxi = 0 ; for i=1:length(csum) ; if csum(i) > maxi ; maxi = csum(i) ; end ; drawdowns(i) = maxi-csum(i) ; end
maxdd = max(drawdowns) ; maxloss = 50*maxdd ; 
clear allequities
equity = 50000 ; tracks = alltracks{x,y,z} ; maxequity = 0 ; minequity = 1000000 ; 
for i=1:length(tracks)
    allowed = equity*.5 ; % percentage of equity willing to be risked
    maxtick = (allowed/maxdd)/4 ; % value of a single tick with this risk
    units = floor(maxtick/12.5) ; % round off to number of units 
    equity = equity + (tracks(i)-.25)*units*12.5 ; % compute new equity
    if equity > maxequity ; maxequity = equity ; end ; if equity < minequity ; minequity = equity ; end
    allequities(i) = equity ; 
end
plot(allequities) ; xlabel('trade #') ; ylabel('equity (start = $50,000)') ; 




