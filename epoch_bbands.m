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
tf=2; 
for ccy=45%:length(currs)
    if(ccy>19)
        cd(['e:\dukas\forex\',timeframes{tf}]) 
    else
        cd(['e:\dukas\dukas_other\',timeframes{tf}]) 
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
    ps = 100; 
    for p=1:length(ps)
    bbands = getbbands(closes,ps(p),2);
    with_bbands = get_bbandcrossovers(closes,bbands); 

    upinds = find(with_bbands==1); upinds(upinds < 201 | upinds > length(closes)-200) = [];
    downinds = find(with_bbands==-1); downinds(downinds < 201 | downinds > length(closes)-200) = []; 
    
    plot(closes) ; hold on ; plot(bbands(1,:),'r') ; plot(bbands(2,:),'g'); 
    plot(upinds,bbands(1,upinds),'ro','LineWidth',3);  plot(downinds,bbands(2,downinds),'go','LineWidth',3); 
    set(gca,'XTickLabel',[]) ; xlabel('time'); ylabel('EURUSD'); title('EURUSD 30 minute chart'); 
    
    %vline(upinds,'r') ; vline(downinds,'g'); 
    
    clear up_epochs
    for i=1:length(upinds)
        up_epochs(i,:) = closes(upinds(i)-200:upinds(i)+200);   
    end
    clear down_epochs
    for i=1:length(downinds)
        down_epochs(i,:) = closes(downinds(i)-200:downinds(i)+200);   
    end
    
    size(down_epochs)
    
    depochs = cumsum(diff(down_epochs(:,201:300),1,2),2)/0.0001;
    uepochs = cumsum(diff(up_epochs(:,201:300),1,2),2)/0.0001;   
    
    [h,pv,ci,stats] = ttest2(uepochs(:,65),depochs(:,65)); 
    
    errorbar(squeeze(mean(depochs(:,:),1)),squeeze(std(depochs(:,:),0,1))/sqrt(788),'g'); hold on ; hline(0,'k'); 
    errorbar(squeeze(mean(uepochs(:,:),1)),squeeze(std(uepochs(:,:),0,1))/sqrt(788),'r');
    ylabel('pips'); xlabel('time (hours)'); set(gca,'XTick',1:25:100,'XTickLabel',round(((1:25:100)*30)/60));
    title('p < 0.001'); 
    
    downs(p,:) = mean(cumsum(diff(up_epochs(:,201:400),1,2),2));
    ups(p,:) = mean(cumsum(diff(down_epochs(:,201:400),1,2),2)); 
    
    
    %plot(mean(cumsum(diff(up_epochs(:,101:200),1,2),2))) ; hold on;
    %plot(mean(cumsum(diff(down_epochs(:,101:200),1,2),2)),'r') ; hline(0,'k'); 
    end
    currdowns(ccy,:,:) = downs./scales(ccy); 
    currups(ccy,:,:) = ups./scales(ccy); 

end      
