clear all ; close all ; 
cd('C:\shared\mkt\Forex') ; 
fid = fopen('EURUSD.txt') ; 
data = textscan(fid,'%s %s %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
fclose(fid) ; 
d5 = data{6} ; d5 = d5(1:5:end) ; 
d5 = d5(length(d5)/2:end) ; 
%d5 = d5(end-50000:end) ; 

% make a moving average:
clear mvg p mean
for p = 1:700  
    for i=p+1:length(d5)
        if i==p+1
            mvg(i,p) = mean(d5(i-p+1:i)) ; 
        else
            mvg(i,p) = mvg(i-1,p) + d5(i)/(p) - d5(i-p)/(p) ;
        end  
    end
end

diffmvg = diff(mvg) ; 

clear allprofits allntrades alltracks
m1s = 1:20:500 ; m2s = 1:20:500 ; 
m1count = 1 ; 
for m1=1:5:700 ; m2count = 1 ; disp(m1) ; 
    for m2=1:5:700 ; tpcount = 1 ; 
        for tp=0.0001:0.0001:0.005
            trade = false ; entry = 0 ; buy = false ; sell = false ; profit = 0 ; ntrades = 0 ; alltrades = [] ; sellentries = [] ; buyentries = [] ; sellexits = [] ; buyexits = [] ; 
            for i=701:length(d5)-1
                if trade == false
                    if d5(i) > mvg(i,m1) % buy mode
                        if d5(i) < mvg(i,m2) && d5(i-1) > mvg(i,m2) % cross down
                            trade = true ; entry = d5(i) ; buy = true ; ntrades = ntrades + 1  ;
                        end
                    elseif d5(i) < mvg(i,m1) % sell mode
                        if d5(i) > mvg(i,m2) && d5(i-1) < mvg(i,m2) % cross up
                            trade = true ; entry = d5(i) ; sell = true ; ntrades = ntrades + 1 ; 
                        end
                    end
                elseif trade == true
                    if buy
                        if d5(i)-entry > tp || d5(i)-entry <-tp
                            buy = false ; profit = profit + d5(i)-entry ; trade = false ; alltrades(length(alltrades)+1) = d5(i) - entry ; 
                        end
                    elseif sell
                        if  entry-d5(i) > tp || entry-d5(i) <-tp
                            sell = false ; profit = profit + entry-d5(i) ; trade = false ; alltrades(length(alltrades)+1) = entry - d5(i) ; 
                        end                   
                    end
                end    
            end ; 
            allprofits(m1count,m2count,tpcount) = profit ; allntrades(m1count,m2count,tpcount) = ntrades ; alltracks{m1count,m2count,tpcount} = alltrades ;
            allbuyentries{m1count,m2count,tpcount} = buyentries ; allsellentries{m1count,m2count,tpcount} = sellentries ; allbuyexits{m1count,m2count,tpcount} = buyexits ; allsellexits{m1count,m2count,tpcount} = sellexits ; 
            tpcount = tpcount + 1 ; 
        end ;m2count = m2count + 1 ;
    end ; m1count = m1count + 1 ; 
end ;
for i=1:size(allprofits,3) ; subplot(5,10,i); imagesc(squeeze(allprofits(:,:,i)./allntrades(:,:,i)),[-.0006,.0006]) ; title(i) ; end

figure,plot(cumsum(alltracks{2,31,4}))

