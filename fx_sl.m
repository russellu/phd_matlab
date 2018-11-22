clear all ; close all ; 
currs =   {'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDSGD','AUDUSD','CADCHF','CADHKD','CADJPY','CHFJPY','CHFSGD','EURAUD','EURCAD','EURCHF','EURDKK','EURGBP','EURHKD','EURJPY','EURNZD','EURPLN','EURSEK','EURSGD','EURTRY','EURUSD','GBPAUD','GBPCAD','GBPCHF','GBPJPY','GBPNZD','GBPUSD','HKDJPY','NZDCAD','NZDCHF','NZDJPY','NZDUSD','SGDJPY','USDCAD','USDCHF','USDCNH','USDDKK','USDHKD','USDJPY','USDMXN','USDNOK','USDSGD','USDTRY','USDZAR','XAGUSD','XAUUSD','ZARJPY'} ; 
spreads = [ 2.5,     2.5,     0.95,    2.33,    5,       0.9,     1.95,    16.0     1.05,    1.5,     9,       1.9,     1.85,    1.0      4.0      0.85,    18.0,    0.6,     3.6,     20.0,    19.0,    5.25,    9.0,     0.25,    3.35,    3.5,     2.2,     1.6,     5.0,     0.96,    28.0,    3.0,     2.1,     2.1,     1.05,    2.5,     0.9,     1.05,    2.5,     4.9,     3.0,     0.35,    50.0,     22.0,    2.0,     6.0,     100.0,   3.05,    28.0,    0.7   ]; 
[spreadsort,spreadsortind] = sort(spreads,'ascend') ; 
currs = currs(spreadsortind) ; spreads = spreads(spreadsortind) ; 
cd c:/shared/back_test; ls 
for ccy=1:30
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
   
    nincrs = 10 ; 
    stepsz = floor(length(d5)/nincrs) ; clear inds
    for i=1:nincrs ; inds(i,:) = (i-1)*stepsz+1:i*stepsz ; end ; 
    
    mvgs1 = 1:20:600 ; 
    pipthreshs = [0]*mfactor ; 
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
    clear allprofits alltracks allntrades
    for ind=1:size(inds,1)
        for m1=1:length(mvgs1) ; disp(['curr=',currs{ccy},', m1=',num2str(m1),', totalcost=',num2str(totalcost),' ind = ',num2str(ind)]) ;
            for m2=1:length(mvgs1) ; 
                    trade = false ; entry = 0 ; buy = false ; sell = false ; profit = 0 ; ntrades = 0 ; upl = 0 ; 
                    alltrades = [] ; sellentries = [] ; buyentries = [] ; sellexits = [] ; buyexits = [] ;trackprofits = [] ; 
                    for i=inds(ind,1):inds(ind,end)
                        if mod(i,100)== 0 ; trackprofits(length(trackprofits)+1) = profit ; end
                        if trade == false
                                if d5(i) < mvg(i,(m1)) && d5(i-1) > mvg(i,(m1)) % cross down
                                    trade = true ; entry = d5(i) ; buy = true ; ntrades = ntrades + 1  ;
                                end
                                if d5(i) > mvg(i,(m1)) && d5(i-1) < mvg(i,(m1)) % cross up
                                    trade = true ; entry = d5(i) ; sell = true ; ntrades = ntrades + 1 ;  
                                end
                        elseif trade == true
                            
                            if buy
                                upl = d5(i)-entry-totalcost ; 
                                if d5(i) > mvg(i,(m2)) && d5(i-1) < mvg(i,(m2)) || upl < -100*mfactor
                                    buy = false ; profit = profit + d5(i)-entry-totalcost ; trade = false ; alltrades(length(alltrades)+1) = d5(i) - entry - totalcost ; 
                                    upl = 0 ; 
                                end
                            elseif sell
                                upl = entry-d5(i)-totalcost ; 
                                if d5(i) < mvg(i,(m2)) && d5(i-1) > mvg(i,(m2)) || upl < -100*mfactor
                                    sell = false ; profit = profit + entry-d5(i)-totalcost ; trade = false ; alltrades(length(alltrades)+1) = entry - d5(i) -totalcost ; 
                                    upl = 0 ; 
                                end
                            end
                        end
                    end
                    allprofits(m1,m2,ind) = profit+upl ; 
                    allntrades(m1,m2,ind) = ntrades ; 
                    alltracks{m1,m2,ind} = alltrades ; 
            end
        end
    end
    savestuffs{1} = allprofits ; savestuffs{2} = allntrades ; savestuffs{3} = alltracks ; 
    save(['chunk_',currs{ccy}],'savestuffs') ; 
    mps = allprofits./allntrades ; 
    cvar = mean(mps,3)./std(mps,0,3) ;
    
    for i=1:size(alltracks,1)
        for j=1:size(alltracks,2)
            for k=1:size(alltracks,3)
                if ~isempty(alltracks{i,j,k})
                    corrs(i,j,k) = corr2(cumsum(alltracks{i,j,k}),1:length(cumsum(alltracks{i,j,k}))) ; 
                end
            end
        end
    end
    figure,
    subplot(1,2,1) ; imagesc(squeeze(mean(mps,3)),[-mfactor*10,mfactor*10]) ; title(num2str(sum(mps(:)>0)./numel(mps))) ; 
    subplot(1,2,2) ; imagesc(squeeze(mean(allprofits,3)),[-mfactor*400,mfactor*400]) ; suptitle(currs{ccy}) ; 
end

%for i=1:size(alltracks,3) ; ps(i,:) = imresize(cumsum(squeeze(alltracks{8,4,i})),[1,50]) ; end ; plot(mean(ps)) ; 


