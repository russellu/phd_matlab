clear all ; close all; 
currs =   {'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDUSD','CADCHF','CADJPY','CHFJPY','EURAUD','EURCAD','EURCHF',...
    'EURGBP','EURJPY','EURUSD','GBPCHF','GBPJPY','GBPUSD','NZDCHF','NZDJPY','NZDUSD','SGDJPY','USDCAD','USDCHF','USDJPY','USDSGD'} ; 
spreads = [ 2.5,     2.5,     0.95,    2.33,        0.9,     1.95,      1.05,    1.5 ,   1.9,     1.85,    1.0  ...
0.85,    0.6,  0.25, 2.2,     1.6,  0.96,  2.1,     2.1,     1.05,    2.5,     0.9,     1.05,    0.35,    2.0]; 
pipscales = [.0001,.0001,.01,.0001,.0001,.0001,.01,.01,.0001,.0001,.0001,.0001,.01,.0001,.0001,.01,.0001,.0001,.01,.0001,.01,.0001,.0001,.01,.0001]; 


cd C:\shared\dukas_mindata; 
for ccy=1:length(currs)
    namei = dir([currs{ccy},'*','Ask*']);
    disp(namei.name); 
    fid = fopen(namei.name) ; 
    dkdata = textscan(fid,'%s %s %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
    fclose(fid) ; 
    d5 = dkdata{5} ;

    williams = zeros(size(d5)); 
    period = 100; 
    for i=period+1:length(d5)
        maxi = max(d5(i-period:i)); 
        mini = min(d5(i-period:i)); 
        williams(i) = ((maxi-d5(i)) ./ (maxi-mini))*-100 ;
    end
    
    clear allprofits allntrades alltrades 
    tps = (1:3:30)*pipscales(ccy);  stops = (5:10:200)*pipscales(ccy); stops = -stops; 
    allptrackers = zeros(length(tps),length(stops),1880); 
    for tpcount=1:length(tps) ; disp(['tpcount = ',num2str(tpcount), ' currency = ', currs{ccy}]); 
        for stopcount=1:length(stops)
            ptrackercount = 1; 
            tp = tps(tpcount); stop = stops(stopcount); trades = []; 
            ntrades = 0; buy = false; sell = false; entry = 0; profit = 0; 
            for i=2000:1880000
                if mod(i,1000)==0
                    allptrackers(tpcount,stopcount,ptrackercount) = profit ; ptrackercount = ptrackercount + 1; 
                end
                if buy==false && sell==false
                    if williams(i) > -5 && williams(i-1) < -5
                        buy = true; entry = d5(i); ntrades = ntrades + 1; 
                    elseif williams(i) < -95 && williams(i-1) > -95
                        sell = true; entry = d5(i) ; ntrades = ntrades + 1; 
                    end
                elseif buy==true
                    if entry-d5(i) > tp || entry-d5(i) < stop
                        buy = false; profit = profit + entry-d5(i); trades(length(trades)+1) = entry-d5(i) - 0.35*pipscales(ccy) - spreads(ccy)*pipscales(ccy); 
                    end
                elseif sell==true
                    if d5(i)-entry > tp || d5(i)-entry < stop
                        sell = false; profit = profit + entry-d5(i); trades(length(trades)+1) = d5(i)-entry - 0.35*pipscales(ccy) - spreads(ccy)*pipscales(ccy);  
                    end
                end
            end
            allprofits(tpcount,stopcount) = profit; allntrades(tpcount,stopcount) = ntrades; alltrades{tpcount,stopcount} = trades; 
        end
    end
    currprofits(ccy,:,:) = allprofits./pipscales(ccy); currntrades(ccy,:,:) = allntrades; currptrackers(ccy,:,:,:) = allptrackers./pipscales(ccy);
end

for i=1:25 ; for j=1:10 ; for k=1:20 ; currcorrs(i,j,k) = corr2(squeeze(currptrackers(i,j,k,1:1879))',1:1879) ; end ; end ; end 
[sv,si] = sort(squeeze(mean(mean(currcorrs,2),3)),'descend'); ncurrs = 20; 
figure, subplot(2,2,1); imagesc(squeeze(mean(currprofits(si(1:ncurrs),:,:),1)./mean(currntrades(si(1:ncurrs),:,:),1)),[0,3]); colormap jet; 
mcurrps = squeeze(mean(currptrackers(si(1:ncurrs),:,:,:),1)); 
for i=1:10 ; for j=1:20 ; corrs(i,j) = corr2(1:1879,squeeze(mcurrps(i,j,1:1879))'); end ; end
subplot(2,2,2) ; imagesc(corrs,[0,1]) ; colormap jet; 
subplot(2,2,3) ; imagesc(squeeze(mean(currntrades(si(1:ncurrs),:,:))),[0,1500]); 
