% predicting brownian motion : you know that it is more likely to reach a
% price level close to where it is currently than one far away...
% simple parameter space search: sell uptick, buy downtick (stop and tp)
% hline strategy...price breaks range, buy or sell. optimize TP/SL with
% parameter space search

% optimize zcurr as well...

close all ; clear all ;
%[retarr,datas] = get_dukascopy('C:\Users\Acer\Documents\fxfiles1m','EURUSD_UTC_5 Mins_Bid_2012.01.01_2014.08.18.csv') ;
%[retarr,datas] = get_dukascopy('C:\Users\Acer\Documents\fxfiles','AUDUSD_UTC_1 Min_Bid_2007.01.01_2014.08.22.csv') ;
[retarr,datas] = get_dukascopy('C:\Users\Acer\Documents\fxfiles','GBPJPY_UTC_1 Min_Bid_2008.08.18_2014.08.22.csv') ;
icount = 1 ;
opens = retarr{1}{3}' ; 
tbuff = 0.1 ; 
buffcount = 1 ;
zcurrthresh = 0 ; 

% larger buffer and larger windowsize means more pips per trade?
for buff = 1:10
    
        %tbuff = tbuff + .0005 ; 
        zcurrthresh = zcurrthresh-.015 ; 
        icount = 1 ; 
    for wsize = 40
        disp(['buff = ',num2str(buff), ' wsize= ',num2str(wsize)]) ; 
    tpthresh = 0.0001 ; 
    for x = 1:20
       % disp(x) ;
        tpthresh = tpthresh + .0004 ;
        stopthresh = .0002 ;
    %    tbuff = tbuff + .0001 ; 
        for y = 1:20
            stopthresh = stopthresh + .0004 ; 
            rval = floor(rand*(size(opens,2)/2)) ; 
            b = opens(rval:floor(rval+(60)*24*5)) ; % 20 days of 1m data
            buy = false ; 
            sell = false ; 
            profits = 0 ; 
            pcount = 1 ; 
            uplarr = [] ;
            ploton = false ;
            squee = [] ;
            tcount = 1 ; 
            range = false ; 
            rangearr = [] ; 
            buythresh = 0 ; 
            sellthresh = 0 ; 
            rangeset = false ; 
            tradecount = 1 ;
            trades = [] ;
            refract = true ; 

                        for i=1000+wsize+1:size(b,2)   
                           squee(size(squee,2)+1) = (max(b(i-wsize:i)) - min(b(i-wsize:i)))/wsize;     
                           zsquee = zscore(squee) ; 
                           zcurr = zsquee(size(zsquee,2)) ; 
                           mx = max(b(i-wsize:i)) ; 
                           mn = min(b(i-wsize:i)) ; 

                           if zcurr < zcurrthresh && refract% && ~buy && ~sell
                               rangearr(size(rangearr,2)+1) = 1 ; 
                                   buythresh = mx + tbuff ;
                                   sellthresh = mn - tbuff ;
                                   refract = false ; 
                           else 
                               rangearr(size(rangearr,2)+1) = 0 ;
                           end             

                           if b(i) > buythresh && ~buy && ~sell && ~refract
                               buy = true ; 
                               entry = b(i) ; 
                           elseif b(i) < sellthresh && ~sell && ~buy && ~refract
                               sell = true ;
                               entry = b(i) ; 
                           end

                           if buy
                               uplarr(size(uplarr,2)+1) = b(i) - entry ; 
                               currupl = uplarr(size(uplarr,2)) ; 
                               if currupl > tpthresh || currupl < -stopthresh
                                   buy = false ; 
                                   uplarr = [] ;
                                   trades(tradecount) = currupl ; 
                                   tradecount = tradecount + 1 ;
                                   refract = true ; 
                               end
                           elseif sell
                               uplarr(size(uplarr,2)+1) = entry - b(i) ; 
                               currupl = uplarr(size(uplarr,2)) ; 
                               if currupl > tpthresh || currupl < -stopthresh
                                   sell = false ; 
                                   uplarr = [] ;
                                   trades(tradecount) = currupl ; 
                                   tradecount = tradecount + 1 ; 
                                   refract = true ; 
                               end
                           end

                           ploton = false ;
                           if ploton && tcount > 1000

                               subplot(2,2,1) ;
                               plot(b(i-100:i)) ; title(['price = ',num2str(b(i))]) ;  
                               if buy
                                   hline(entry,'g') ; 
                               elseif sell
                                   hline(entry,'r') ; 
                               end
                               hline(mx) ; hline(mn) ; 
                               hline(buythresh,'b') ; hline(sellthresh,'b') ; 
                               subplot(2,2,2) ; plot(zsquee) ; title(['zsquee = ',num2str(zsquee(size(zsquee,2))), ' range = ', num2str(mn), ':',num2str(mn), ' = ',num2str(mx-mn)]) ; 
                               subplot(2,2,3) ; imagesc(rangearr(size(rangearr,2)-100:size(rangearr,2)),[0,1]) ; 
                               subplot(2,2,4) ; bar(trades) ; 
                               getframe ; 
                               pause(.01) ;
                           end
                           tcount = tcount + 1 ;
                        end 
                xs(x,y) = mean(trades) ;
        end
    end
    allims(buffcount,icount,:,:) = xs ; 
    icount = icount + 1 ; 
    end
    buffcount = buffcount + 1 ; 
end

for i=1:10 ; 
    subplot(4,3,i) ; 
    imagesc(imfilter(squeeze(mean(allims(i,:,:,:),2)),fspecial('gaussian',3)),[-.1,.1]) ;
end



