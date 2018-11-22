close all ; clear all ;
%[retarr,datas] = get_dukascopy('C:\Users\Acer\Documents\fxfiles1m','EURUSD_UTC_5 Mins_Bid_2012.01.01_2014.08.18.csv') ;
%[retarr,datas] = get_dukascopy('C:\Users\Acer\Documents\fxfiles','AUDUSD_UTC_1 Min_Bid_2007.01.01_2014.08.22.csv') ;
[retarr,datas] = get_dukascopy('C:\Users\Acer\Documents\fxfiles','EURUSD_UTC_5 Mins_Bid_2008.08.18_2014.08.22.csv') ;

opens = retarr{1}{3}' ; 
 
%{
a = rand(1,10000)-.5 ; 
b(1) = a(1) ; 
for i=2:size(a,2)
    b(i) = b(i-1) + a(i) ; 
    
end
%}
rval = rand*size(opens,2)/2 ; 
b = opens(rval:rval+(60/5)*24*120) ; % 30 days of 5m data

wcount = 1 ; 

% predicting brownian motion : you know that it is more likely to reach a
% price level close to where it is currently than one far away...

% simple parameter space search: sell uptick, buy downtick (stop and tp)

buy = false ; 
sell = false ; 
up = false ; 
down = false ; 
upcount = 0 ; 
downcount = 0 ; 

downs = 0 ;
ups = 0 ; 

darrcount = 1 ;
uarrcount = 1 ;

enthresh = 3 ;  
clthresh = 3 ; 

profits = 0 ; 
pcount = 1 ; 

uplarr = [] ;
ploton = false ;

for wsize = 25:25

squee = [] ;

            for i=1000+wsize+1:size(b,2)

               d = b(i) - b(i-1) ; 
               
               %for f=1:wsize
                 squee(size(squee,2)+1) = (max(b(i-wsize:i)) - min(b(i-wsize:i)))/wsize;
                % squee(f,:) = (max(b(i-f:i)) - min(b(i-f:i)))/f;
               %end
               if d > 0 && up
                   upcount = upcount + 1 ; 
               elseif d < 0 && down 
                   downcount = downcount + 1 ; 
               elseif d < 0 && up
                   downcount = 1 ;
               elseif d > 0 && down
                   upcount = 1 ; 
               end 
               if d > 0 
                   up = true ; 
                   downcount = 0 ;
               elseif d < 0 
                   down = true ; 
                   upcount = 0 ;
               end

               if buy
                   uplarr(size(uplarr,2)+1) = b(i) - entry ;
                   if downcount >= clthresh
                       profits(pcount) = uplarr(size(uplarr,2)) ; 
                       pcount = pcount + 1 ;
                       %disp('exiting buy') ; 
                       buy = false ; 
                   end
               elseif sell
                   uplarr(size(uplarr,2)+1) = entry - b(i) ; 
                   if upcount >= clthresh
                       profits(pcount) = uplarr(size(uplarr,2)) ;
                       pcount = pcount + 1 ;
                       %disp('exiting sell')  ;
                       sell = false ; 
                   end    
               end

               if buy==false && sell==false
                   if upcount > enthresh
                       buy = true ;
                       entry = b(i) ;     
                       uplarr = [] ; 
                       %disp('entering buy') ;
                   elseif downcount > enthresh
                       sell = true ; 
                       entry = b(i) ; 
                       uplarr = [] ; 
                       %disp('entering sell') ; 
                   end
               end

               ploton = true ;
               if ploton
                   subplot(2,1,1) ;
                   plot(b(i-100:i)) ; title(['upcount = ',num2str(upcount), ' downcount = ',num2str(downcount)]) ; 
                  % subplot(2,2,2) ; 
                 %  bar(profits) ; 
                   %subplot(2,2,3) ; imagesc(uplarr) ; 
                   subplot(2,1,2) ; plot(zscore(squee)) ; 
                   getframe ; 
                   pause(.1) ;
               end
            end
            
            subplot(6,7,wcount) ; plot(squee(1:2000)) ; wcount = wcount + 1 ;
end