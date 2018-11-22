clear all ; close all ; 
cd('C:\shared\mkt') ; 
fid = fopen('ES.txt') ; 
data = textscan(fid,'%s %s %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
fclose(fid) ; offsetcount = 1 ; 
for offset=1:5 ; 
d5 = data{6} ; d5 = d5(offset:10:end) ; d5 = d5(250000:end) ; 
%d5 = d5(end-50000:end) ; 

% make a moving average:
clear mvg p mean
for p = 1:1000  
    for i=p+1:length(d5)
        if i==p+1
            mvg(i,p) = mean(d5(i-p+1:i)) ; 
        else
            mvg(i,p) = mvg(i-1,p) + d5(i)/(p) - d5(i-p)/(p) ;
        end  
    end
end
disp(offset) ; 
clear allprofits allntrades alltracks
m1s = 1:10:100 ; m2s = 1:10:100 ; tps = 0.5:3:30 ; 
m1count = 1 ; 
for m1=1:100:1000 ; m2count = 1 ; %disp(m1) ; 
    for m2=1:100:1000 ; tpcount = 1 ; 
        for tp=0.5:5:20
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

offsetprofits(:,:,:,offsetcount) = allprofits ; offsetntrades(:,:,:,offsetcount) = allntrades ;
offsetcount = offsetcount + 1 ;
end