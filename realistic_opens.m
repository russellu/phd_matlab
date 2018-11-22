clear all ; close all ; 
timeframes = {'minute_0015'};
            
currs = {
         'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDSGD','AUDUSD','CADCHF','CADHKD','CADJPY','CHFJPY','CHFSGD','EURAUD','EURCAD','EURCHF','EURDKK',...
         'EURGBP','EURHKD','EURJPY','EURNOK','EURNZD','EURPLN','EURSEK','EURSGD','EURTRY','EURUSD','GBPAUD','GBPCAD','GBPCHF','GBPJPY','GBPNZD',...
         'GBPUSD','HKDJPY','NZDCAD','NZDCHF','NZDJPY','NZDUSD','SGDJPY','TRYJPY','USDCAD','USDCHF','USDCNH','USDDKK','USDHKD','USDJPY','USDMXN','USDNOK',...
         'USDPLN','USDSEK','USDSGD','USDTRY','USDZAR','ZARJPY'} ; 

pip_values = [
              0.78,  1.04,   0.94,  0.73,   0.76, 1,      1.04,   0.13,   0.94,   0.94,  0.76,  0.77,  0.78,  1.04,  0.17,...
              1.41,  0.13,   0.94,  0.13,   0.73, 0.29,   0.12,   0.76,   0.25,  1,     0.77,  0.78,  1.04,  0.94,  0.73,...
              1,     0.0094, 0.78,  1.04,   0.94, 1.00,   0.94,   0.94,   0.78,   1.04,  0.16,  0.16,  0.13,  0.94,  0.055, 0.13,...
              0.29,  0.12,  0.76,   0.25, 0.084,  0.94    ];
          
scales = [         
          0.0001, 0.0001, 0.01,   0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.01,   0.01,   0.0001, 0.0001, 0.0001, 0.0001, 0.0001,...
          0.0001, 0.0001, 0.01,   0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.01,   0.0001,...
          0.0001, 0.01,   0.0001, 0.0001, 0.01,   0.0001, 0.01,   0.01,   0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.01,   0.0001, 0.0001,...
          0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.01];

for ccy=1:length(currs)

    cd(['e:\dukas\fulldata\',timeframes{1}]) 

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

    threshs = [50,100,150]; 
    periods = [100,200,500]; 
    
    clear allprofit allntrades alltradetimes; 
    for t1=1:length(periods)
        for t2=1:length(threshs)
       
    %bbands = getbbands(opens,periods(t1),stds(t2)); 
    %indi = zeros(1,length(bbands)); 
    %for i=2:length(bbands) ; if opens(i) > bbands(1,i)&& opens(i-1) < bbands(1,i-1) ; ; indi(i) = -1; elseif opens(i) < bbands(2,i)&& opens(i-1) > bbands(2,i-1) ; indi(i) = 1 ; end ; end % 
    %bbanddiff = bbands(1,:) - bbands(2,:); 
    %bbanddiff = bbanddiff/scales(ccy); 
            
    %aroon = getaroon(opens,periods(t1)); 
    %indi = zeros(1,length(aroon)); 
    %for i=2:length(aroon) ; if aroon(i) > threshs(t2)  ; indi(i) = -1 ; elseif aroon(i) < -threshs(t2) ; indi(i) = 1; end ; end %&& aroon(i-1) < 98 && aroon(i-1) > -98 
    
    cci = getcci(lows,highs,opens,periods(t1));
    indi = zeros(1,length(cci)) ; 
    for i=2:length(cci) ; if cci(i) > threshs(t2) ; indi(i) = 1 ; elseif cci(i) < -threshs(t2) ; indi(i) = -1; end ; end
    
    close_ask = dkask{5}; 
    close_bid = dkbid{5}; 
    
    in_buy = false; in_sell = false; entry_price = 0; profit = 0; ntrades = 0; saved_trades = []; entry_index = 0; tradetime = 0; nightprofit = 0; 
    for i=1200:length(opens)        
        if ~in_buy && ~in_sell 
            if indi(i) == 1 
                in_buy = true; entry_price = close_ask(i); ntrades = ntrades + 1; entry_index = i; 
            elseif indi(i) == -1 
                in_sell = true; entry_price = close_bid(i); ntrades = ntrades + 1; entry_index = i; 
            end           
        elseif in_buy
            upl = (close_bid(i) - entry_price)/scales(ccy) ;
            if indi(i) == -1 || upl > 125 
                profit = profit + close_bid(i) - entry_price ; in_buy = false; tradetime = tradetime + i-entry_index; 
                saved_trades(length(saved_trades)+1) = (close_bid(i) - entry_price) ./ scales(ccy) ; 
            end
        elseif in_sell
            upl = (entry_price - close_ask(i))/scales(ccy);
            if indi(i) == 1 || upl > 125 
                profit = profit + entry_price - close_ask(i); in_sell = false; tradetime = tradetime + i-entry_index;
                saved_trades(length(saved_trades)+1) = (entry_price - close_ask(i)) ./ scales(ccy) ; 
            end
        end    

    end       
    
    allprofit(t1,t2) = profit; allntrades(t1,t2) = ntrades; alltradetimes(t1,t2) = tradetime/ntrades; 
    if ~isempty(saved_trades); allcumsums(t1,t2,:) = imresize(cumsum(saved_trades),[1,500]); else; allcumsums(t1,t2,:) = zeros(1,500); end
    
        end
    end
        
    currprofits(ccy,:,:) = allprofit./scales(ccy);
    currntrades(ccy,:,:) = allntrades; 
    currcumsums(ccy,:,:,:) = allcumsums; 
    currtradetimes(ccy,:,:) = alltradetimes; 
end      

for c=1:52 ;for i=1:size(currprofits,2); for j=1:size(currprofits,3); corrs(c,i,j) = corr2(squeeze(currcumsums(c,i,j,:))',1:500) ; end ; end; end
corrs(isnan(corrs)) = 0 ;
[sv,si] = sort(squeeze(mean(mean(corrs,2),3)),'descend'); 

count=1;
for i=1:size(currprofits,2)
    for j=1:size(currprofits,3)
        %[sv,si] = sort(corrs(:,i,j),'descend'); 
        subplottight(size(currprofits,2),size(currprofits,3),count); plot(squeeze(mean(currcumsums(si(1:25),i,j,:),1))); set(gca,'XTickLabel',[],'YTickLabel',[]); 
        count = count + 1 ; 
        mcorrs(i,j) = corr2(1:500,squeeze(mean(currcumsums(si(1:25),i,j,:),1))'); 
    end
end

