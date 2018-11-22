close all ; clear all ;
[retarr,datas] = get_dukascopy('C:\Users\Acer\Documents\fxfiles1m','EURUSD_UTC_5 Mins_Bid_2012.01.01_2014.08.18.csv') ;
opens = retarr{1}{3}' ; 



for total = 1:100
    total
    %{
a = rand(1,10000)-.5 ; 
b(1) = a(1) ; 
for i=2:size(a,2)
    b(i) = b(i-1) + a(i) ; 
    
end
%}
    rval = rand*size(opens,2)/2 ; 
b = opens(rval:rval+(60/5)*24*30) ; % 30 days of 5m data


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

for et = 1:9
    enthresh = et ; 
    for ct = 1:9
        clthresh = ct ; 
        pcount = 1 ; profits = 0 ; uplarr = [] ; 
           % clthresh = 5 ; 
           % enthresh = 5 ; 
            clear trialprofs
            for i=2:size(b,2)

               d = b(i) - b(i-1) ; 

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

               ploton = false ;
               if ploton
                   subplot(2,2,1) ;
                   plot(b(1:i)) ; title(['upcount = ',num2str(upcount), ' downcount = ',num2str(downcount)]) ; 
                   subplot(2,2,2) ; 
                   bar(profits) ; 
                   subplot(2,2,3) ; imagesc(uplarr) ; 
                   getframe ; 
                   pause(.01) ;
               end
            end
    
        %subplot(3,4,1) ;
        %plot(profits) 
        ps(et,ct) = mean(profits) ;
    end
end
subplot(10,10,total) ;
imagesc(ps) ;
allps(total,:,:) = ps ; 
end
