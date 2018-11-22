clear all ; close all; 
currs =   {'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDUSD','CADCHF','CADJPY','CHFJPY','EURAUD','EURCAD','EURCHF',...
    'EURGBP','EURJPY','EURUSD','GBPCHF','GBPJPY','GBPUSD','NZDCHF','NZDJPY','NZDUSD','SGDJPY','USDCAD','USDCHF','USDJPY','USDSGD'} ; 
spreads = [ 2.5,     2.5,     0.95,    2.33,        0.9,     1.95,      1.05,    1.5 ,   1.9,     1.85,    1.0  ...
0.85,    0.6,  0.25, 2.2,     1.6,  0.96,  2.1,     2.1,     1.05,    2.5,     0.9,     1.05,    0.35,    2.0]; 
pipscales = [.0001,.0001,.01,.0001,.0001,.0001,.01,.01,.0001,.0001,.0001,.0001,.01,.0001,.0001,.01,.0001,.0001,.01,.0001,.01,.0001,.0001,.01,.0001]; 

cd C:\shared\dukas_mindata; 
for ccy=1:length(currs) 
namei = dir([currs{ccy},'*','Ask*']);
fid = fopen(namei.name) ; 
dkdata = textscan(fid,'%s %s %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
fclose(fid) ; 
d5 = dkdata{5} ;

rsilength = 35;
rsi_up = zeros(size(d5)); 
rsi_down = zeros(size(d5)); 
for i=rsilength+1:length(d5)
    seg = d5(i-rsilength:i); 
    diffs = diff(seg); 
    npos = sum(diffs>0); 
    nneg = sum(diffs<0); 
    rsi_up(i) = npos;
    rsi_down(i) = nneg; 
    
end
smthup = getmvg(rsi_up,5); 
smthdown = getmvg(rsi_down,5); 
rsi = 100-100./(1+smthup./smthdown); 

clear allprofits allntrades alltrades; 
tps = (1:2:20)*pipscales(ccy);  stops = (5:10:200)*pipscales(ccy); stops = -stops; 
for tpcount=1:length(tps) ;  disp(['tpcount = ',num2str(tpcount), ' currency = ', currs{ccy}]);
    for stopcount=1:length(stops); 
        tp = tps(tpcount); stop = stops(stopcount); trades = []; 
        ntrades = 0; buy = false; sell = false; entry = 0; profit = 0; 
        for i=2000:length(d5)
            if buy==false && sell==false
            if  rsi(i)>70
                buy = true; entry = d5(i); ntrades = ntrades + 1; 
            elseif rsi(i) <30
                sell = true; entry = d5(i) ; ntrades = ntrades + 1; 
            end
            elseif buy==true
                if entry-d5(i) > tp || entry-d5(i) < stop
                    buy = false; profit = profit + entry-d5(i); trades(length(trades)+1) = entry-d5(i); 
                end
            elseif sell==true
                if d5(i)-entry > tp || d5(i)-entry < stop
                    sell = false; profit = profit + entry-d5(i); trades(length(trades)+1) = d5(i)-entry;  
                end
            end
        end
        allprofits(tpcount,stopcount) = profit; allntrades(tpcount,stopcount) = ntrades; alltrades{tpcount,stopcount} = trades; 
    end
end
currprofits(ccy,:,:) = allprofits./pipscales(ccy); 
currntrades(ccy,:,:) = allntrades; 
end
figure,imagesc(squeeze(mean(currprofits./currntrades)));






