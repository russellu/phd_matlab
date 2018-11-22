clear all ; close all ; 
cd('C:\shared\mkt') ; 
fid = fopen('ES.txt') ; 
data = textscan(fid,'%s %s %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
fclose(fid) ; 
d5 = data{6} ; d5 = d5(15:30:end) ; 

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

clear allprofits allntrades ; 
for tlim=1:15 ; disp(tlim) ; plimcount = 1 ; 
    for plim=1:.25:5 ;
        m1s = 1:5:100 ; m2s = 1:5:100 ; tps = 0.5:2:30 ; 
        m1 = m1s(16) ; m2 = m2s(13) ; tp = tps(11) ; 
        trade = false ; entry = 0 ; buy = false ; sell = false ; profit = 0 ; ntrades = 0 ; waitforbuy = false ; waitforsell = false ; realentry = 0 ; tcounter = 0 ; alltrades = [] ; 
        for i=300:length(d5)
            if trade == false
                if d5(i) > mvg(i,m1) % buy mode
                    if d5(i) < mvg(i,m2) && d5(i-1) > mvg(i,m2) % cross down
                        trade = true ; entry = d5(i) ; waitforbuy = true ; tcounter = 0 ; 
                    end
                elseif d5(i) < mvg(i,m1) % sell mode
                    if d5(i) > mvg(i,m2) && d5(i-1) < mvg(i,m2) % cross up
                        trade = true ; entry = d5(i) ; waitforsell = true ; tcounter = 0 ; 
                    end
                end
            elseif trade == true
                if waitforbuy % wait for it to go higher
                    if d5(i) - entry <= -plim
                        ntrades = ntrades + 1  ; %buyentries(length(buyentries)+1) = i ; 
                        realentry = d5(i) ; waitforbuy = false ; buy = true ; 
                    else tcounter = tcounter + 1 ; 
                    end
                elseif waitforsell % wait for it to go lower
                    if entry - d5(i) <= -plim
                        ntrades = ntrades + 1 ; %sellentries(length(sellentries)+1) = i ; 
                        realentry = d5(i) ; waitforsell = false ; sell = true ; 
                    else tcounter = tcounter + 1 ; 
                    end
                end
                if buy
                    if d5(i) - entry > tp || d5(i) - entry < -tp
                        buy = false ; profit = profit + d5(i)-realentry ; trade = false ; alltrades(length(alltrades)+1) = d5(i) - realentry ; 
                    end
                elseif sell
                    if entry - d5(i) > tp || entry - d5(i) < -tp
                        sell = false ; profit = profit + realentry-d5(i) ; trade = false ; alltrades(length(alltrades)+1) = realentry - d5(i) ; 
                    end                   
                end
                if tcounter > tlim ; trade = false ; waitforbuy = false ; waitforsell = false ; end 
            end    
        end 
            allprofits(tlim,plimcount) = profit ; allntrades(tlim,plimcount) = ntrades ; plimcount = plimcount + 1 ; tracks{tlim,plimcount} = alltrades ; 
    end
end
imagesc(allprofits./allntrades)

