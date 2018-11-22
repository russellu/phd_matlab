clear all ; close all; 
currs =   {'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDUSD','CADCHF','CADJPY','CHFJPY','EURAUD','EURCAD','EURCHF',...
    'EURGBP','EURJPY','EURUSD','GBPCHF','GBPJPY','GBPUSD','NZDCHF','NZDJPY','NZDUSD','SGDJPY','USDCAD','USDCHF','USDJPY','USDSGD'} ; 
spreads = [ 2.5,     2.5,     0.95,    2.33,        0.9,     1.95,      1.05,    1.5 ,   1.9,     1.85,    1.0  ...
0.85,    0.6,  0.25, 2.2,     1.6,  0.96,  2.1,     2.1,     1.05,    2.5,     0.9,     1.05,    0.35,    2.0]; 
pipscales = [.0001,.0001,.01,.0001,.0001,.0001,.01,.01,.0001,.0001,.0001,.0001,.01,.0001,.0001,.01,.0001,.0001,.01,.0001,.01,.0001,.0001,.01,.0001]; 

cd C:\shared\dukas_mindata; 
for ccy =1:length(currs); 
namei = dir([currs{ccy},'*','Ask*']);
fid = fopen(namei.name) ; 
dkdata = textscan(fid,'%s %s %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
fclose(fid) ; 
d5 = dkdata{5} ;
stochperiods = 100; 
for stochperiod=1:length(stochperiods)
    period = stochperiods(stochperiod) ; 
    stoch = zeros(size(d5)); 
    for i=period+1:length(d5)
        stoch(i) = 100*(d5(i)-min(d5(i-period:i)))./(max(d5(i-period:i))-min(d5(i-period:i))); 
    end
    stoch(isnan(stoch)) = 0.5; 
    stochmvg = getmvg(stoch,25); 

    %{
    mvg1 = getmvg(d5,50); 
    bbands = zeros(2,length(d5)); 
    for i=150+1:length(d5)
        stdi = std(d5(i-50:i)); 
        bbands(1,i) = mvg1(i)+stdi*3.25;
    end
    %}
    clear allprofits allntrades alltrades; 
    tps = (1:2:20)*pipscales(ccy);  stops = (5:10:200)*pipscales(ccy); stops = -stops; 
    for tpcount=1:length(tps) ;  disp(['tpcount = ',num2str(tpcount), ' currency = ', currs{ccy}]);
        for stopcount=1:length(stops); 
            tp = tps(tpcount); stop = stops(stopcount); trades = []; 
            ntrades = 0; buy = false; sell = false; entry = 0; profit = 0; 
            for i=2000:length(d5)
                if buy==false && sell==false
                if  stochmvg(i-1) > 95 && stochmvg(i) < 95% (d5(i-1)>bbands(2,i) && d5(i)< bbands(2,i)) &&
                    sell = true; entry = d5(i); ntrades = ntrades + 1; 
                elseif stochmvg(i-1) <5 && stochmvg(i) > 5% (d5(i-1)<bbands(1,i) && d5(i) > bbands(1,i)) && 
                    buy = true; entry = d5(i) ; ntrades = ntrades + 1; 
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

    currprofits(ccy,stochperiod,:,:) = allprofits./pipscales(ccy); 
    currntrades(ccy,stochperiod,:,:) = allntrades; 
end
end
mprofits = currprofits./currntrades; 
imagesc(squeeze(mean(mprofits,1)))


%{
corrs = zeros(size(allntrades)) ; meanprofits = zeros(size(allntrades)); 
for i=1:size(alltrades,1)
    for j=1:size(alltrades,2)
        if ~isempty(alltrades{i,j})
            corrs(i,j) = corr2(cumsum(alltrades{i,j}),1:length(cumsum(alltrades{i,j}))); 
            meanprofits(i,j) = mean(alltrades{i,j}); 
        end
    end
end
meanprofits = meanprofits./pipscales(ccy); 
allcorrs(ccy,:,:) = corrs; allmeanprofits(ccy,:,:) = meanprofits;  
%}









