clear all ; close all ; 
timeframes = {'minute_0015','minute_0030','minute_0060','minute_0240','minute_1440'};

currs = {'AUSIDXAUD','BRENTCMDUSD','BTCUSD','CHEIDXCHF','COPPERCMDUSD','DEUIDXEUR','ESPIDXEUR','EUSIDXEUR','FRAIDXEUR','GASCMDUSD',...
         'GBRIDXGBP','HKGIDXHKD','JPNIDXJPY','LIGHTCMDUSD','USA30IDXUSD','USA500IDXUSD','USATECHIDXUSD','XAGUSD','XAUUSD',... 
         'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDSGD','AUDUSD','CADCHF','CADHKD','CADJPY','CHFJPY','CHFSGD','EURAUD','EURCAD','EURCHF','EURDKK',...
         'EURGBP','EURHKD','EURJPY','EURNOK','EURNZD','EURPLN','EURRUB','EURSEK','EURSGD','EURTRY','EURUSD','GBPAUD','GBPCAD','GBPCHF','GBPJPY','GBPNZD',...
         'GBPUSD','HKDJPY','NZDCAD','NZDCHF','NZDJPY','NZDUSD','SGDJPY','TRYJPY','USDCAD','USDCHF','USDCNH','USDDKK','USDHKD','USDJPY','USDMXN','USDNOK',...
         'USDPLN','USDRUB','USDSEK','USDSGD','USDTRY','USDZAR','ZARJPY'} ; 

pip_values = [77.01, 10000, 10000, 104.22, 1000, 122.96, 122.95, 122.95, 122.95, 1000,...
              140.8, 12.74,  0.94,  10000,  100,  100,    100,    100,    100,...
              0.78,  1.04,   0.94,  0.73,   0.76, 1,      1.04,   0.13,   0.94,   0.94,  0.76,  0.77,  0.78,  1.04,  0.17,...
              1.41,  0.13,   0.94,  0.13,   0.73, 0.29,   0.02,   0.12,   0.76,   0.25,  1,     0.77,  0.78,  1.04,  0.94,  0.73,...
              1,     0.0094, 0.78,  1.04,   0.94, 1.00,   0.94,   0.94,   0.78,   1.04,  0.16,  0.16,  0.13,  0.94,  0.055, 0.13,...
              0.29,  0.017,  0.12,  0.76,   0.25, 0.084,  0.94    ];
          
scales = [
          0.01,   0.01,   0.1,    0.01,   0.001,  0.01,   0.01,   0.01,   0.01,  0.001,...
          0.01,   0.01,   0.01,   0.01,   0.01,   0.01,   0.01,   0.01,   0.01, ...
          0.0001, 0.0001, 0.01,   0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.01,   0.01,...
          0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.01,   0.0001, 0.0001,...
          0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.01,...
          0.0001, 0.0001, 0.01,   0.0001, 0.0001, 0.01,   0.0001, 0.01,   0.01,   0.0001,...
          0.0001, 0.0001, 0.0001, 0.0001, 0.01,   0.0001, 0.0001, 0.0001, 0.0001, 0.0001,...
          0.0001, 0.0001, 0.0001, 0.01];

      for tf=1%:length(timeframes)
for ccy=1:length(currs)
    if(ccy>19)
        cd(['C:\shared\dukas\forex\',timeframes{tf}]) 
    else
        cd(['C:\shared\dukas_other\',timeframes{tf}]) 
    end
    namei = dir([currs{ccy},'*','Bid*']); disp(namei.name); 
    fid = fopen(namei.name) ; 
    dkbid = textscan(fid,'%s %f %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
    fclose(fid) ; 
    % ask
    namei = dir([currs{ccy},'*','Ask*']);
    fid = fopen(namei.name) ; 
    dkask = textscan(fid,'%s %f %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
    fclose(fid) ; 

    lengths = [cellfun(@length,dkask),cellfun(@length,dkbid)]; 
    shortl = min(lengths); 
    dkbid{1} = dkbid{1}(1:shortl); 
    dkbid{2} = dkbid{2}(1:shortl);
    dkbid{3} = dkbid{3}(1:shortl); 
    dkbid{4} = dkbid{4}(1:shortl); 
    dkbid{5} = dkbid{5}(1:shortl);
    dkbid{6} = dkbid{6}(1:shortl);
    dkbid{7} = dkbid{7}(1:shortl);       
    dkask{1} = dkask{1}(1:shortl); 
    dkask{2} = dkask{2}(1:shortl);
    dkask{3} = dkask{3}(1:shortl);
    dkask{4} = dkask{4}(1:shortl);
    dkask{5} = dkask{5}(1:shortl);
    dkask{6} = dkask{6}(1:shortl);
    dkask{7} = dkask{7}(1:shortl);

    bid_datetimes = dkbid{1}; 
    bid_datetimes = datenum(bid_datetimes,'yyyy.mm.dd HH:MM:SS'); 

    ask_datetimes = dkask{1}; 
    ask_datetimes = datenum(ask_datetimes,'yyyy.mm.dd HH:MM:SS'); 

    opens = single((dkbid{2}+dkask{2})/2) ; 
    highs = single((dkbid{3}+dkask{3})/2) ;  
    lows = single((dkbid{4}+dkask{4})/2) ; 
    closes = single((dkbid{5}+dkask{5})/2) ;  
    vols = single((dkbid{6}+dkask{6})/2); 
    spreads = single(dkask{5} - dkbid{5}); 
    timestamps = double(ask_datetimes); 

    bbandperiods = [40,80,150,250]; 
    bbandstds = [2.5,3]; 
    
    
    
    for p1=1:length(bbandperiods)
        for p2=1:length(bbandstds)
            bbands = getbbands(closes,bbandperiods(p1),bbandstds(p2)); %
            
            with_bbands = get_bbandcrossovers(closes,bbands); 
            against_bbands = with_bbands*-1; 

            close_ask = dkask{5}; 
            close_bid = dkbid{5}; 
            entry_ts = against_bbands;
            exit_ts = against_bbands;
            clear allprofit allntrades alltradetimes; 

            in_buy = false; in_sell = false; entry_price = 0; profit = 0; ntrades = 0; saved_trades = []; entry_index = 0; tradetime = 0; 
            for i=60:length(closes)        
                if ~in_buy && ~in_sell
                    if entry_ts(i) == 1
                        in_buy = true; entry_price = close_ask(i); ntrades = ntrades + 1; entry_index = i; 
                    elseif entry_ts(i) == -1
                        in_sell = true; entry_price = close_bid(i); ntrades = ntrades + 1; entry_index = i; 
                    end           
                elseif in_buy
                    if exit_ts(i) == -1
                        profit = profit + close_bid(i) - entry_price; in_buy = false; tradetime = tradetime + i-entry_index; 
                        saved_trades(length(saved_trades)+1) = (close_bid(i) - entry_price) ./ scales(ccy);
                    end
                elseif in_sell
                    if exit_ts(i) == 1
                        profit = profit + entry_price - close_ask(i); in_sell = false; tradetime = tradetime + i-entry_index; 
                        saved_trades(length(saved_trades)+1) = (entry_price - close_ask(i)) ./ scales(ccy);
                    end
                end        
            end            
            allprofit(p1,p2) = profit; allntrades(p1,p2) = ntrades; alltradetimes(p1,p2) = tradetime/ntrades; 
            if ~isempty(saved_trades); allcumsums(p1,p2,:) = imresize(cumsum(saved_trades),[1,500]); else; allcumsums(p1,p2,:) = zeros(1,500); end
        end
    end
    currprofits(ccy,:,:) = allprofit./scales(ccy);
    currntrades(ccy,:,:) = allntrades; 
    currcumsums(ccy,:,:,:) = allcumsums; 
    currtradetimes(ccy,:,:) = alltradetimes; 
    disp(currs{ccy}); 
          
end      

for i=1:54
    for j=1:size(currprofits,2)
        for k=1:size(currprofits,3)
            corrs(i,j,k) = corr2(1:500,squeeze(currcumsums(i,j,k,:))'); 
        end
    end
end
figure,
imagesc(squeeze(mean(corrs,1)))

name = ['bband2_bband2','_',timeframes{tf}]; 
fxstruct.currprofits = currprofits;
fxstruct.currntrades = currntrades;
fxstruct.currcumsums = currcumsums;
fxstruct.currtradetimes = currtradetimes; 
cd c:/shared/fxstructs; 
save(name,'fxstruct'); 

      end
%{
name = ['bband2_bband2','_',timeframes{1}]; 
fxstruct.currprofits = currprofits;
fxstruct.currntrades = currntrades;
fxstruct.currcumsums = currcumsums;
fxstruct.currtradetimes = currtradetimes; 
cd c:/shared/fxstructs; 
save(name,'fxstruct'); 
%}


