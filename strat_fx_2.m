clear all ; close all ; 
currs = {'EURCAD','EURCHF','EURGBP','EURJPY','EURUSD','GBPCHF','GBPJPY','GBPUSD','NZDUSD','USDCAD','USDCHF','USDJPY'} ; 

for ccy=1:length(currs) ; 
    curr = currs{ccy} ; 
    cd('C:\shared\mkt\Forex') ; 
    fid = fopen([curr,'.txt']) ; 
    data = textscan(fid,'%s %s %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
    fclose(fid) ; 
    tshiftcount = 1 ;
    clear allprofits alltracks allntrades
    for tshift = 1:2:15
        orig = data{6} ; orig = orig(tshift:15:end) ; mtimes = data{1} ; mtimes = mtimes(tshift:15:end) ; 
        orig = orig(length(orig)/2:end) ; mtimes = mtimes(length(mtimes)/2:end) ; 
        nincrs = 20 ; 
        d5incr = length(orig)./nincrs ; 
        for d5i=1:nincrs ; 
            startindex = (d5i-1)*d5incr + 1 ; endindex = startindex+d5incr-1 ; 
            d5 = orig(startindex:endindex) ; 

            mvgs1 = 1:25:700 ; 
            mvgs2 = 1:25:700 ; 

            % make a moving average:
            clear mvg p mean
            for p = 1:length(mvgs1)  
                for i=mvgs1(p)+1:length(d5)
                    if i==mvgs1(p)+1
                        mvg(i,p) = mean(d5(i-mvgs1(p)+1:i)) ; 
                    else
                        mvg(i,p) = mvg(i-1,p) + d5(i)/(mvgs1(p)) - d5(i-mvgs1(p))/(mvgs1(p)) ;
                    end  
                end
            end

            for m1=1:length(mvgs1) ; disp(['tshift=',num2str(tshift),'d5i=',num2str(d5i),' m1=',num2str(m1)]) ;
                for m2=1:length(mvgs2) ; 
                    for m3=1:length(mvgs2) ; 
                        trade = false ; entry = 0 ; buy = false ; sell = false ; profit = 0 ; ntrades = 0 ; alltrades = [] ; sellentries = [] ; buyentries = [] ; sellexits = [] ; buyexits = [] ; 
                        for i=700:length(d5)
                            if trade == false
                                if d5(i) > mvg(i,(m1)) % buy mode
                                    if d5(i) < mvg(i,(m2)) && d5(i-1) > mvg(i,(m2)) % cross down
                                        trade = true ; entry = d5(i) ; buy = true ; ntrades = ntrades + 1  ;
                                    end
                                elseif d5(i) < mvg(i,(m1)) % sell mode
                                    if d5(i) > mvg(i,(m2)) && d5(i-1) < mvg(i,(m2)) % cross up
                                        trade = true ; entry = d5(i) ; sell = true ; ntrades = ntrades + 1 ;  
                                    end
                                end
                            elseif trade == true
                                if buy
                                    if d5(i) > mvg(i,(m3)) && d5(i-1) < mvg(i,(m3))
                                        buy = false ; profit = profit + d5(i)-entry ; trade = false ; alltrades(length(alltrades)+1) = d5(i) - entry ; 
                                    end
                                elseif sell
                                    if d5(i) < mvg(i,(m3)) && d5(i-1) > mvg(i,(m3))
                                        sell = false ; profit = profit + entry-d5(i) ; trade = false ; alltrades(length(alltrades)+1) = entry - d5(i) ; 
                                    end                   
                                end
                            end    
                        end  
                        allprofits(m1,m2,m3,d5i,tshiftcount) = profit ; 
                        allntrades(m1,m2,m3,d5i,tshiftcount) = ntrades ; 
                        alltracks{m1,m2,m3,d5i,tshiftcount} = alltrades ;
                    end
                end  
            end
        end  
        tshiftcount = tshiftcount + 1 ; 
    end
    alldata{1} = allprofits ; alldata{2} = allntrades ; alldata{3} = alltracks ; 
    save(['notp_alldata_',curr],'alldata') ; 
end
