clear all ; close all ; 
timeframes = {'minute_0015','minute_0030','minute_0060','minute_0240','minute_1440'};

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

tf=2; 
for ccy=1:length(currs)
   % if(ccy>19)
        cd(['e:\dukas\forex\',timeframes{tf}]) 
  %  else
   %     cd(['e:\dukas\dukas_other\',timeframes{tf}]) 
   % end
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
        
    ps = 60; 
    clear ts; 
    for p=1:length(ps)
    
        
        %chande
        %{
        period=ps(p);
        chande = zeros(1,length(closes)); 
        oc_diffs = highs - lows; 
        for i=period+1:length(closes)
            cdiffs = oc_diffs(i-period:i); 
            pos_sum = sum(cdiffs(cdiffs>0));
            neg_sum = abs(sum(cdiffs(cdiffs<0))); 
            chande(i) = 100*((pos_sum - neg_sum)/(pos_sum+neg_sum));          
        end
        indi = zeros(1,length(closes)); 
        for i=2:length(chande) ; if chande(i) > 30 && chande(i-1)<30 ; indi(i) = -1 ; elseif chande(i) < -30 && chande(i-1) > -30 ; indi(i) = 1; end ; end
        %}
        
        %{
        highlows = zeros(2,length(closes)); 
        for i=ps(p)+1:length(closes)
           highlows(1,i) = min(closes(i-ps(p):i));
           highlows(2,i) = max(closes(i-ps(p):i)); 
        end        
        highlows(1,:) = getmvg(highlows(1,:),round(ps(p)/8));
        highlows(2,:) = getmvg(highlows(2,:),round(ps(p)/8)); 
        indi = zeros(1,length(highlows)); 
        for i=2:length(highlows) ; if closes(i) > highlows(2,i) && closes(i-1) < highlows(2,i-1) ; indi(i) = -1; elseif closes(i) < highlows(1,i) && closes(i-1) > highlows(1,i-1) ; indi(i) = 1 ; end ; end
        %}
        
        %bbands (1.2)
        %bbands = getbbands(closes,ps(p),2.25);
        %indi = zeros(1,length(bbands)); 
        %for i=2:length(bbands) ; if closes(i) > bbands(1,i) && closes(i-1) < bbands(1,i-1) ; indi(i) = -1; elseif closes(i) < bbands(2,i) && closes(i-1) > bbands(2,i-1) ; indi(i) = 1 ; end ; end
        %bbanddiff = bbands(1,:) - bbands(2,:); 
        
        % stds
        %{
        stds = zeros(1,length(closes)); 
        for i=51:length(stds)
           stds(i) = std(closes(i-50:i)); 
        end
        %}
        
        % highlows
        highlows = highs - lows; 
        
        
        %cci (1.2)
        %cci = getcci(lows,highs,closes,ps);
        %indi = zeros(1,length(cci)) ; 
        %for i=2:length(cci) ; if cci(i) > 150 && cci(i-1)<150 ; indi(i) = -1 ; elseif cci(i) < -150 && cci(i-1) > -150 ; indi(i) = 1; end ; end
        
        %rsi (1.1)
        %rsi = getrsi(opens,closes,ps(p)); 
        %indi = zeros(1,length(rsi)); 
        %for i=2:length(rsi) ; if rsi(i) > 60 && rsi(i-1)<60 ; indi(i) = -1 ; elseif rsi(i) < 40 && rsi(i-1) > 40 ; indi(i) = 1; end ; end

        %stoch (1.2)
        %stoch = getstoch(closes,ps(p));
        %indi = zeros(1,length(stoch)); 
        %for i=2:length(stoch) ; if stoch(i) > 90 && stoch(i-1)<90 ; indi(i) = -1 ; elseif stoch(i) < 10 && stoch(i-1) > 10 ; indi(i) = 1; end ; end
        
        %macd (0.3)
        %macd = getmacd(closes,ps(p),ps(p)/2)./scales(ccy); 
        %indi = zeros(1,length(macd)); 
        %for i=2:length(macd) ; if macd(i) > 5 && macd(i-1)<5 ; indi(i) = 1 ; elseif macd(i) < 5 && macd(i-1) > 5 ; indi(i) = -1; end ; end
        
        %aroon (0.9) 
        %aroon = getaroon(closes,ps(p)); 
        %indi = zeros(1,length(aroon)); 
        %for i=2:length(aroon) ; if aroon(i) > 90 && aroon(i-1)<90 ; indi(i) = -1 ; elseif aroon(i) < -90 && aroon(i-1) > -90 ; indi(i) = 1; end ; end
        
        %avg = ((cci/5) +(rsi-50) +(stoch-50) +(chande'))/4;  
        %indi = zeros(1,length(avg)); 
        %for i=2:length(avg) ; if avg(i) > 30 && avg(i-1)<30 ; indi(i) = -1 ; elseif avg(i) < -30 && avg(i-1) > -30 ; indi(i) = 1; end ; end
       
        mvg = getmvg(closes,100); mvgdiff = (closes' - mvg)/scales(ccy); 
        indi = zeros(1,length(mvg));      
        for i=2:length(mvgdiff) ; if mvgdiff(i) > 30 && mvgdiff(i-1) < 30 ; indi(i) = 1 ; elseif mvgdiff(i) < -30 && mvgdiff(i-1) > -30 ; indi(i) = -1 ; end ; end
        
        
        upinds = find(indi==1); upinds(upinds < 201 | upinds > length(closes)-200) = [];
        downinds = find(indi==-1); downinds(downinds < 201 | downinds > length(closes)-200) = []; 

        clear up_epochs up_epochinds
        for i=1:length(upinds)
            up_epochs(i,:) = closes(upinds(i)-200:upinds(i)+200);   
            up_epochinds(i,:) = upinds(i)-200:upinds(i)+200; 
        end
        clear down_epochs down_epochinds
        for i=1:length(downinds)
            down_epochs(i,:) = closes(downinds(i)-200:downinds(i)+200);   
            down_epochinds(i,:) = downinds(i)-200:downinds(i)+200; 
        end

        depochs = cumsum(diff(down_epochs(:,201:300),1,2),2)/scales(ccy);
        uepochs = cumsum(diff(up_epochs(:,201:300),1,2),2)/scales(ccy);   
        
        %{
        [svd,sid] = sort(mean(depochs(:,1:30),2),'ascend'); 
        [svu,siu] = sort(mean(uepochs(:,1:30),2),'descend'); 

        clear up_aroon_epochs up_macd_epochs up_vol_epochs up_spread_epochs up_bband_epochs up_cci_epochs up_rsi_epochs up_std_epochs up_highlow_epochs
        for i=1:size(up_epochinds)
           up_aroon_epochs(i,:) = aroon(up_epochinds(i,:));  
           up_macd_epochs(i,:) = macd(up_epochinds(i,:)); 
           up_vol_epochs(i,:) = vols(up_epochinds(i,:)); 
           up_spread_epochs(i,:) = spreads(up_epochinds(i,:)); 
           up_bband_epochs(i,:) = bbanddiff(up_epochinds(i,:)); 
           up_cci_epochs(i,:)= cci(up_epochinds(i,:)); 
           up_rsi_epochs(i,:)= rsi(up_epochinds(i,:)); 
           up_std_epochs(i,:) = stds(up_epochinds(i,:)); 
           up_highlow_epochs(i,:) = highlows(up_epochinds(i,:)); 
        end
        
        clear down_aroon_epochs down_macd_epochs down_vol_epochs down_spread_epochs down_bband_epochs down_cci_epochs down_rsi_epochs down_std_epochs down_highlow_epochs
        for i=1:size(down_epochinds)
           down_aroon_epochs(i,:) = aroon(down_epochinds(i,:));  
           down_macd_epochs(i,:) = macd(down_epochinds(i,:)); 
           down_vol_epochs(i,:) = vols(down_epochinds(i,:)); 
           down_spread_epochs(i,:) = spreads(down_epochinds(i,:)); 
           down_bband_epochs(i,:) = bbanddiff(down_epochinds(i,:)); 
           down_cci_epochs(i,:) = cci(down_epochinds(i,:)); 
           down_rsi_epochs(i,:) = rsi(down_epochinds(i,:)); 
           down_std_epochs(i,:) = stds(down_epochinds(i,:)); 
           down_highlow_epochs(i,:) = highlows(down_epochinds(i,:)); 
        end
        %}
        %plot(squeeze(mean(down_aroon_epochs(sid(1:200),:),1))); hold on ; plot(squeeze(mean(down_aroon_epochs(sid(end-200:end),:),1)));
        %plot(squeeze(mean(down_macd_epochs(sid(1:200),:),1))); hold on ; plot(squeeze(mean(down_macd_epochs(sid(end-200:end),:),1)));
        %plot(squeeze(mean(down_vol_epochs(sid(1:200),:),1))); hold on ; plot(squeeze(mean(down_vol_epochs(sid(end-200:end),:),1)));
        %plot(squeeze(mean(down_spread_epochs(sid(1:200),:),1))); hold on ; plot(squeeze(mean(down_spread_epochs(sid(end-200:end),:),1)));
        %figure,
        %subplot(1,2,1);
        %plot(squeeze(mean(down_bband_epochs(sid(1:200),:),1))); hold on ; plot(squeeze(mean(down_bband_epochs(sid(end-200:end),:),1)));
        %subplot(1,2,2); 
        %plot(squeeze(mean(up_bband_epochs(siu(1:200),:),1))); hold on ; plot(squeeze(mean(up_bband_epochs(siu(end-200:end),:),1)));

        %{
        nepochs = 100; 
        bband_pre_downs(ccy,:) = (squeeze(mean(down_bband_epochs(sid(1:nepochs),:),1)) - squeeze(mean(down_bband_epochs(sid(end-nepochs:end),:),1))) ./scales(ccy); 
        bband_pre_ups(ccy,:) = (squeeze(mean(up_bband_epochs(siu(1:nepochs),:),1)) - squeeze(mean(up_bband_epochs(siu(end-nepochs:end),:),1))) ./scales(ccy); 
        
        vol_pre_downs(ccy,:) = squeeze(mean(down_vol_epochs(sid(1:nepochs),:),1)) - squeeze(mean(down_vol_epochs(sid(end-nepochs:end),:),1)); 
        vol_pre_ups(ccy,:) = squeeze(mean(up_vol_epochs(siu(1:nepochs),:),1)) - squeeze(mean(up_vol_epochs(siu(end-nepochs:end),:),1)); 
        
        cci_pre_downs(ccy,:) = (squeeze(mean(down_cci_epochs(sid(1:nepochs),:),1)) - squeeze(mean(down_cci_epochs(sid(end-nepochs:end),:),1))) ./scales(ccy); 
        cci_pre_ups(ccy,:) = (squeeze(mean(up_cci_epochs(siu(1:nepochs),:),1)) - squeeze(mean(up_cci_epochs(siu(end-nepochs:end),:),1))) ./scales(ccy); 
        
        rsi_pre_downs(ccy,:) = (squeeze(mean(down_rsi_epochs(sid(1:nepochs),:),1)) - squeeze(mean(down_rsi_epochs(sid(end-nepochs:end),:),1))) ./scales(ccy); 
        rsi_pre_ups(ccy,:) = (squeeze(mean(up_rsi_epochs(siu(1:nepochs),:),1)) - squeeze(mean(up_rsi_epochs(siu(end-nepochs:end),:),1))) ./scales(ccy); 
        
        stds_pre_downs(1,ccy,:) = squeeze(mean(down_std_epochs(sid(1:nepochs),:),1)) ./ scales(ccy); 
        stds_pre_downs(2,ccy,:) = squeeze(mean(down_std_epochs(sid(end-nepochs:end),:),1)) ./scales(ccy); 
        stds_pre_ups(1,ccy,:) = squeeze(mean(up_std_epochs(siu(1:nepochs),:),1)) ./ scales(ccy); 
        stds_pre_ups(2,ccy,:) = squeeze(mean(up_std_epochs(siu(end-nepochs:end),:),1)) ./scales(ccy); 
        
        highlow_pre_downs(ccy,:) = (squeeze(mean(down_highlow_epochs(sid(1:nepochs),:),1)) - squeeze(mean(down_highlow_epochs(sid(end-nepochs:end),:),1))) ./scales(ccy); 
        highlow_pre_ups(ccy,:) = (squeeze(mean(up_highlow_epochs(siu(1:nepochs),:),1)) - squeeze(mean(up_highlow_epochs(siu(end-nepochs:end),:),1))) ./scales(ccy); 

        %}
        %pre_downs(ccy,:) = squeeze(mean(down_aroon_epochs(sid(1:200),:),1)) - squeeze(mean(down_aroon_epochs(sid(end-200:end),:),1));
        %pre_ups(ccy,:) = squeeze(mean(up_aroon_epochs(siu(1:200),:),1)) - squeeze(mean(up_aroon_epochs(siu(end-200:end),:),1));
        
        for i=1:99
           [h,pv,ci,stats] = ttest2(uepochs(:,i),depochs(:,i));  
           ts(p,i) = stats.tstat; 
           
           [h,pv,ci,stats] = ttest(uepochs(:,i));  
           uts(p,i) = stats.tstat; 
           
           [h,pv,ci,stats] = ttest(depochs(:,i));  
           dts(p,i) = stats.tstat; 
                      
        end

        downs(p,:) = mean(cumsum(diff(up_epochs(:,201:400),1,2),2));
        ups(p,:) = mean(cumsum(diff(down_epochs(:,201:400),1,2),2)); 

    end
    currdowns(ccy,:,:) = downs./scales(ccy); 
    currups(ccy,:,:) = ups./scales(ccy); 
    currts(ccy,:,:) = ts; 
    curruts(ccy,:,:) = uts;
    currdts(ccy,:,:) = dts; 
end      

%{
plot(squeeze(mean(cci_pre_downs))) ; hold on; 
plot(squeeze(mean(cci_pre_ups))) ;

plot(squeeze(mean(rsi_pre_downs))) ; hold on; 
plot(squeeze(mean(rsi_pre_ups))) ;

plot(squeeze(mean(bband_pre_downs))) ; hold on; 
plot(squeeze(mean(bband_pre_ups))) ;

plot(squeeze(mean(vol_pre_downs))) ; hold on; 
plot(squeeze(mean(vol_pre_ups))) ;

plot(squeeze(mean(stds_pre_downs))) ; hold on; 
plot(squeeze(mean(stds_pre_ups))) ;

plot(squeeze(mean(highlow_pre_downs))) ; hold on; 
plot(squeeze(mean(highlow_pre_ups))) ;

plot(squeeze(mean(stds_pre_downs(1,:,:),2))) ; hold on ; plot(squeeze(mean(stds_pre_downs(2,:,:),2)));
plot(squeeze(mean(stds_pre_ups(1,:,:),2))) ; hold on ; plot(squeeze(mean(stds_pre_ups(2,:,:),2)));
%}
%{
subplot(1,2,1); 
plot(squeeze(mean(curruts,1))) ; hold on ; hline(0,'k'); 
plot(squeeze(mean(currdts,1)))

subplot(1,2,2) ; plot(squeeze(mean(currups,1)));  hold on ; hline(0,'k'); 
plot(squeeze(mean(currdowns,1))); 
%}
