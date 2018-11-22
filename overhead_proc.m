clear all ; close all ; 
currs =   {'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDSGD','AUDUSD','CADCHF','CADHKD','CADJPY','CHFJPY','CHFSGD','EURAUD','EURCAD','EURCHF','EURDKK','EURGBP','EURHKD','EURJPY','EURNZD','EURPLN','EURSEK','EURSGD','EURTRY','EURUSD','GBPAUD','GBPCAD','GBPCHF','GBPJPY','GBPNZD','GBPUSD','HKDJPY','NZDCAD','NZDCHF','NZDJPY','NZDUSD','SGDJPY','USDCAD','USDCHF','USDCNH','USDDKK','USDHKD','USDJPY','USDMXN','USDNOK','USDSGD','USDTRY','USDZAR','XAGUSD','XAUUSD','ZARJPY'} ; 
spreads = [ 2.5,     2.5,     0.95,    2.33,    5,       0.9,     1.95,    16.0     1.05,    1.5,     9,       1.9,     1.85,    1.0      4.0      0.85,    18.0,    0.6,     3.6,     20.0,    19.0,    5.25,    9.0,     0.25,    3.35,    3.5,     2.2,     1.6,     5.0,     0.96,    28.0,    3.0,     2.1,     2.1,     1.05,    2.5,     0.9,     1.05,    2.5,     4.9,     3.0,     0.35,    50.0,     22.0,    2.0,     6.0,     100.0,   3.05,    28.0,    0.7   ]; 

cd C:\shared\indices_4 ; ls 
for ccy=6:length(currs) ; 
    ccy_current = dir([currs{ccy},'*.csv']) ; 
    if (~isempty(strfind(ccy_current.name,'JPY')) && isempty(strfind(ccy_current.name,'HKD')) ) || ~isempty(strfind(ccy_current.name,'XAG')) || ~isempty(strfind(ccy_current.name,'XAU')) 
        lim1 = -0.4 ; lim2 = 0.4 ; mfactor = 0.01 ; 
    else 
        lim1 = -0.004 ; lim2 = 0.004 ; mfactor = 0.0001 ; 
    end
    totalcost = (spreads(ccy) + 0.35 + 0.9)*mfactor ; 
    fid = fopen(ccy_current.name) ; 
    dkdata = textscan(fid,'%s %s %f %f %f %f %f','delimiter',' ') ; 
    fclose(fid) ; 
    d5 = dkdata{6} ; d5 = d5(1:30:end) ; 
    
    mvgs1 = [5,15,30,50,90,150,220,320,450,600,900] ; 
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
    pipthreshs = [1,5,10,20,40,80]*mfactor ; 
    clear allprofits alltracks allntrades alltrackprofits
    for m1=1:length(mvgs1) ; disp(['curr=',currs{ccy},', m1=',num2str(m1),', totalcost=',num2str(totalcost)]) ;
        for m2=1:length(mvgs1) ; 
            for m3=1:length(mvgs1) ; 
                for pipthresh=1:length(pipthreshs)
                trade = false ; entry = 0 ; buy = false ; sell = false ; profit = 0 ; ntrades = 0 ; 
                alltrades = [] ; sellentries = [] ; buyentries = [] ; sellexits = [] ; buyexits = [] ;trackprofits = [] ; 
                for i=700:length(d5)
                    if mod(i,100)== 0 ; trackprofits(length(trackprofits)+1) = profit ; end
                    if trade == false
                        if d5(i)-pipthreshs(pipthresh) > mvg(i,(m1)) % buy mode
                            if d5(i) < mvg(i,(m2)) && d5(i-1) > mvg(i,(m2)) % cross down
                                trade = true ; entry = d5(i) ; buy = true ; ntrades = ntrades + 1  ;
                            end
                        elseif d5(i)+ pipthreshs(pipthresh) < mvg(i,(m1)) % sell mode
                            if d5(i) > mvg(i,(m2)) && d5(i-1) < mvg(i,(m2)) % cross up
                                trade = true ; entry = d5(i) ; sell = true ; ntrades = ntrades + 1 ;  
                            end
                        end
                    elseif trade == true
                        if buy
                            if d5(i) > mvg(i,(m3)) && d5(i-1) < mvg(i,(m3))
                                buy = false ; profit = profit + d5(i)-entry - totalcost ; trade = false ; alltrades(length(alltrades)+1) = d5(i) - entry - totalcost ; 
                            end
                        elseif sell
                            if d5(i) < mvg(i,(m3)) && d5(i-1) > mvg(i,(m3))
                                sell = false ; profit = profit + entry-d5(i) - totalcost ; trade = false ; alltrades(length(alltrades)+1) = entry - d5(i) -totalcost ; 
                            end                   
                        end
                    end    
                end  
                allprofits(m1,m2,m3,pipthresh) = profit ; 
                allntrades(m1,m2,m3,pipthresh) = ntrades ; 
                alltracks{m1,m2,m3,pipthresh} = alltrades ;
                alltrackprofits(m1,m2,m3,pipthresh,:) = trackprofits ; 
                end
            end
        end  
    end
    
    divp = allprofits./allntrades ; divp(isnan(divp)) = 0 ; 
    for i=1:size(alltrackprofits,1);for j=1:size(alltrackprofits,2);for k=1:size(alltrackprofits,3);for el=1:size(alltrackprofits,4);corrs(i,j,k,el)=corr2(squeeze(alltrackprofits(i,j,k,el,:))',1:size(alltrackprofits,5));end;end;end;end
    corrs(isnan(corrs)) = 0 ; 
    [sv,si] = sort(corrs(:),'descend') ; 
    [ix,iy,iz,it] = ind2sub(size(corrs),si) ; 
    clear hps
    for i=1:50 ;
        hps(i,:) = squeeze(alltrackprofits(ix(i),iy(i),iz(i),it(i),:)) ; 
        % figure,plot(squeeze(alltrackprofits(ix(i),iy(i),iz(i),it(i),:))) ; hold on ;  
    end
    results{1} = allprofits ; results{2} = allntrades ; results{3} = alltracks ; results{4} = alltrackprofits ; 
    save(['results_',currs{ccy}],'results') ; 
end