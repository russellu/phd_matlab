clear all ; close all ; 
timeframes = {'minute_0030'};
%{
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
%}
            
currs = {
         'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDSGD','AUDUSD','CADCHF','CADHKD','CADJPY','CHFJPY','CHFSGD','EURAUD','EURCAD','EURCHF','EURDKK',...
         'EURGBP','EURHKD','EURJPY','EURNOK','EURNZD','EURPLN','EURRUB','EURSEK','EURSGD','EURTRY','EURUSD','GBPAUD','GBPCAD','GBPCHF','GBPJPY','GBPNZD',...
         'GBPUSD','HKDJPY','NZDCAD','NZDCHF','NZDJPY','NZDUSD','SGDJPY','TRYJPY','USDCAD','USDCHF','USDCNH','USDDKK','USDHKD','USDJPY','USDMXN','USDNOK',...
         'USDPLN','USDRUB','USDSEK','USDSGD','USDTRY','USDZAR','ZARJPY'} ; 

pip_values = [
              0.78,  1.04,   0.94,  0.73,   0.76, 1,      1.04,   0.13,   0.94,   0.94,  0.76,  0.77,  0.78,  1.04,  0.17,...
              1.41,  0.13,   0.94,  0.13,   0.73, 0.29,   0.02,   0.12,   0.76,   0.25,  1,     0.77,  0.78,  1.04,  0.94,  0.73,...
              1,     0.0094, 0.78,  1.04,   0.94, 1.00,   0.94,   0.94,   0.78,   1.04,  0.16,  0.16,  0.13,  0.94,  0.055, 0.13,...
              0.29,  0.017,  0.12,  0.76,   0.25, 0.084,  0.94    ];
          
scales = [         
          0.0001, 0.0001, 0.01,   0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.01,   0.01,   0.0001, 0.0001, 0.0001, 0.0001, 0.0001,...
          0.0001, 0.0001, 0.01,   0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.01,   0.0001,...
          0.0001, 0.01,   0.0001, 0.0001, 0.01,   0.0001, 0.01,   0.01,   0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.01,   0.0001, 0.0001,...
          0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.01];

for ccy=1:length(currs)
   % if(ccy>19)
        cd(['e:\dukas\forex\',timeframes{1}]) 
    %else
    %    cd(['e:\dukas\dukas_other\',timeframes{1}]) 
    %end
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

    period = 40;     
    aroon = getaroon(closes, period); % -100 to 100
    cci = getcci(lows, highs, closes, period); % ~ -400 to 400 
    rsi = getrsi(opens,closes,period); % 0 to 100
    stoch = getstoch(closes,period); % 0 to 100
    smvg = getmvg(closes,period); %
    bbands = getbbands(closes,period,2); %
    macd = getmacd(closes,period,period/2)/scales(ccy); % ~ -50 to 50 
    
    filt_period = 350;     
    filt_aroon = getaroon(closes, filt_period); % -100 to 100
    filt_cci = getcci(lows, highs, closes, filt_period); % ~ -400 to 400 
    filt_rsi = getrsi(opens,closes,filt_period); % 0 to 100
    filt_stoch = getstoch(closes,filt_period); % 0 to 100
    filt_smvg = getmvg(closes,filt_period); %
    filt_bbands = getbbands(closes,filt_period,2); %
    filt_macd = getmacd(closes,filt_period,filt_period/2)/scales(ccy); % ~ -50 to 50 

    with_aroon = zeros(1,length(closes)); with_aroon(aroon<-50) = -1; with_aroon(aroon > 50) = 1; 
    against_aroon = with_aroon*-1; 
    
    with_cci = zeros(1,length(closes)); with_cci(cci<-50) = -1; with_cci(cci>50) = 1; 
    against_cci = with_cci*-1; 
    
    with_rsi = zeros(1,length(closes)); with_rsi(rsi<40) = -1; with_rsi(rsi>60) = 1;
    against_rsi = with_rsi*-1; 
    
    with_stoch = zeros(1,length(closes)); with_stoch(stoch<15) = -1; with_stoch(stoch>85) = 1; 
    against_stoch = with_stoch*-1; 
    
    with_smvg = get_crossovers(closes,smvg); 
    against_smvg = with_smvg*-1; 
    
    with_bbands = get_bbandcrossovers(closes,bbands); 
    against_bbands = with_bbands*-1; 
    
    with_macd = zeros(1,length(macd)); with_macd(macd<-30) = 1; with_macd(macd>30) = -1;
    against_macd = with_macd*-1; 
       
    
    filt_with_aroon = zeros(1,length(closes)); filt_with_aroon(aroon<-10) = -1; filt_with_aroon(aroon > 10) = 1; 
    filt_against_aroon = filt_with_aroon*-1; 
    
    filt_with_cci = zeros(1,length(closes)); filt_with_cci(filt_cci<-10) = -1; filt_with_cci(filt_cci>10) = 1; 
    filt_against_cci = filt_with_cci*-1; 
    
    filt_with_rsi = zeros(1,length(closes)); filt_with_rsi(filt_rsi<45) = -1; filt_with_rsi(filt_rsi>55) = 1;
    filt_against_rsi = filt_with_rsi*-1; 
    
    filt_with_stoch = zeros(1,length(closes)); filt_with_stoch(filt_stoch<40) = -1; filt_with_stoch(filt_stoch>60) = 1; 
    filt_against_stoch = filt_with_stoch*-1; 
    
    %{
    filt_with_smvg = get_crossovers(closes,filt_smvg); 
    filt_against_smvg = filt_with_smvg*-1; 
    
    filt_with_bbands = get_bbandcrossovers(closes,filt_bbands); 
    filt_against_bbands = filt_with_bbands*-1; 
    
    filt_with_macd = zeros(1,length(filt_macd)); filt_with_macd(filt_macd<-30) = 1; filt_with_macd(filt_macd>30) = -1;
    filt_against_macd = filt_with_macd*-1; 
    %}
    
    all_filt_indis = [filt_with_aroon; filt_against_aroon; filt_with_cci; filt_against_cci; filt_with_rsi; filt_against_rsi; filt_with_stoch; filt_against_stoch]; 
      %  filt_with_smvg; filt_against_smvg; filt_with_bbands; filt_against_bbands; filt_with_macd; filt_against_macd];
    
    all_indis = [with_aroon; against_aroon; with_cci; against_cci; with_rsi; against_rsi; with_stoch; against_stoch; ...
        with_smvg; against_smvg; with_bbands; against_bbands; with_macd; against_macd];
    
    all_indi_names = {'with_aroon', 'against_aroon', 'with_cci', 'against_cci','with_rsi', 'against_rsi', 'with_stoch', 'against_stoch',...
        'with_smvg', 'against_smvg', 'with_bbands', 'against_bbands', 'with_macd', 'against_macd'};
    
    close_ask = dkask{5}; 
    close_bid = dkbid{5}; 
    
    clear allprofit allntrades alltradetimes; 
    for fp1 = 1:size(all_filt_indis,1) 
    for p1=1:size(all_indis,1)
        for p2=1:size(all_indis,1)
            filt_ts = all_filt_indis(fp1,:); 
            entry_ts = all_indis(p1,:);
            exit_ts = all_indis(p2,:); 
            in_buy = false; in_sell = false; entry_price = 0; profit = 0; ntrades = 0; saved_trades = []; entry_index = 0; tradetime = 0; 
            for i=period:length(closes)        
                if ~in_buy && ~in_sell
                    if entry_ts(i) == 1 && filt_ts(i)==1
                        in_buy = true; entry_price = close_ask(i); ntrades = ntrades + 1; entry_index = i; 
                    elseif entry_ts(i) == -1 && filt_ts(i)==-1
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
            allprofit(fp1,p1,p2) = profit; allntrades(fp1,p1,p2) = ntrades; alltradetimes(fp1,p1,p2) = tradetime/ntrades; 
            if ~isempty(saved_trades); allcumsums(fp1,p1,p2,:) = imresize(cumsum(saved_trades),[1,500]); else; allcumsums(fp1,p1,p2,:) = zeros(1,500); end
        end
    end
    end
    
    currprofits(ccy,:,:,:) = allprofit./scales(ccy);
    currntrades(ccy,:,:,:) = allntrades; 
    currcumsums(ccy,:,:,:,:) = allcumsums; 
    currtradetimes(ccy,:,:,:) = alltradetimes; 
end      


for i=1:size(currprofits,1); for j=1:size(currprofits,2); for k=1:size(currprofits,3) ; for el=1:size(currprofits,4) ; corrs(i,j,k,el) = corr2(squeeze(currcumsums(i,j,k,el,:))',1:500) ; end; end; end; end





      