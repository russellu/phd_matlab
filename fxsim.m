clear all  ; close all ;
%[retarr,datas] = get_dukascopy('c:\users\russ\Documents\dukascopy\eurusd_no_wkds','EURUSD_1 Min_Bid_2009.01.01_2013.08.16.csv') ;

[retarr,datas] = get_dukascopy('C:\Users\Acer\Documents\fxfiles1m','EURUSD_UTC_5 Mins_Bid_2012.01.01_2014.08.18.csv') ;

opens = retarr{1}{3}' ; 
clear totalps allprofs
for scount=1:5
    wcount =1 ;
    for w=2:15:400
        
        clear mps
        for xcount=1:100

            clear allps ;
            clear trades                          
          %  scount = 6; 
           %  w = 10 ;                
          %  ListenChar(2);
            buy = false ; 
            sell = false ; 
            profit = 0 ;
            profits = [0] ; 
            uplarr = [0] ; 
            st = 0 ; 
            entry = 0 ; 
            pcount = 1;  
            i = floor(rand*(size(opens,2)/2))+1000 ;
            istart = i ; 
            repcount = 0 ;
            up = false ;
            down = false ;

            while pcount < 50
                
            %{      
                if GetSecs - st > .00
                
                  
                  subplot(2,2,1) ; 
                   plot(opens(i-150:i)) ; 
                   
                   subplot(2,2,3) ; 
                   bar((profits)) ; title(num2str(sum(profits))) ; %hline(0,'k') ;  
                   subplot(2,2,2) ; 
                   if buy || sell
                       if size(uplarr,2) > 400
                            imagesc(uplarr(size(uplarr,2)-400:size(uplarr,2))) ; 
                       else
                           imagesc(uplarr) ; hline(0) ; 
                       end
                   end
                   getframe ;  
                    

                   st = GetSecs ;    
             
          %}   
                   fsize = w ;
                   if size(uplarr,2) > fsize % if the size of the window is large enough to start monitoring upl
                       if uplarr(size(uplarr,2)) == max(uplarr(size(uplarr,2)-fsize:size(uplarr,2))) % if we are at a local max in the upl count
                            if buy == true % exit the buy
                               buy = false  ;
                               profits(pcount) = profit + uplarr(size(uplarr,2)) ; 
                               uplarr = [0] ; 
                               pcount = pcount + 1 ;            
                            elseif sell == true
                                sell = false ;
                                profits(pcount) = profit + uplarr(size(uplarr,2)) ; 
                                uplarr = [0] ; 
                                pcount = pcount + 1 ; 
                            end

                       end
                   end      
                   
                   % if the trade has not yet been exited, check and update
                   % the upl
                   if buy == true
                      uplarr(size(uplarr,2)+1) =  opens(i)-opens(entry) ;  
                      if uplarr(size(uplarr,2)) < -.9
                          buy = false  ;
                          profits(pcount) = profit + uplarr(size(uplarr,2)) ; 
                          uplarr = [0] ; 
                          pcount = pcount + 1 ; 
                      end
                      if uplarr(size(uplarr,2)) > .9
                          buy = false  ;
                          profits(pcount) = profit + uplarr(size(uplarr,2)) ; 
                          uplarr = [0] ; 
                          pcount = pcount + 1 ; 
                      end
                    elseif sell == true
                      uplarr(size(uplarr,2)+1) = (opens(entry)-opens(i)) ;   
                      if uplarr(size(uplarr,2)) < -.9
                          sell = false ;
                          profits(pcount) = profit + uplarr(size(uplarr,2)) ; 
                          uplarr = [0] ; 
                          pcount = pcount + 1 ; 
                      end
                      if uplarr(size(uplarr,2)) > .9
                          sell = false ;
                          profits(pcount) = profit + uplarr(size(uplarr,2)) ; 
                          uplarr = [0] ; 
                          pcount = pcount + 1 ; 
                      end
                   end
                   i = i+1 ;

          %      end


                if opens(i) - opens(i-1) > 0
                    if up == true ; dircount = dircount + 1 ; else dircount = 0 ; end
                    up = true ; down = false ; 
                elseif opens(i)-opens(i-1) < 0
                    if down == true ; dircount = dircount + 1 ; else dircount = 0 ; end
                    down = true ; up = false ; 
                end


                        if ~buy && ~sell
                            if up && dircount >= scount
                                trades(pcount+1,:) = opens(i-10:i+10) ;
                                sell = true ; entry = i ;
                            elseif down && dircount >= scount
                               trades(pcount+1,:) = opens(i-10:i+10) ;
                               buy = true ; entry = i ;
                            end

                            %{
                            if rand > .5
                                sell = true ; 
                            else 
                                buy = true ; 
                            end
                            %}
                        end    
                        

            end
                close all ; 
             %   ListenChar(1);

            %sum(profits)/size(profits,2) 

            allps = sum(profits) ;
            allsizes = size(profits,2) ; 

            

            x = allps(1) ;
            for i=1:size(allps,2)
                x(i+1) = allps(i) + x(i) ; 


            end
            
           % totalps(scount,wcount) = mean(profits) ; 
            %mean(profits)
            mps(xcount) = mean(profits) ;
            allprofs(xcount,:) = profits ; 
        end 
       % bar(mps) ; title(num2str(mean(mps))) ; 
        
        
            totalps(scount,wcount) = mean(mps) ;
            wcount = wcount + 1; 
            disp(['w = ',num2str(wcount)]) ; 

            
            
            
            
            
        %plot(x) ;
    end
    disp(['scount = ',num2str(scount)]) ; 
end

imagesc(totalps)


% replicate in java and then go to demo
% demo for 2 months + then go live