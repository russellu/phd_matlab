clear all  ; close all ;
%[retarr,datas] = get_dukascopy('c:\users\russ\Documents\dukascopy\eurusd_no_wkds','EURUSD_1 Min_Bid_2009.01.01_2013.08.16.csv') ;

[retarr,datas] = get_dukascopy('C:\Users\Acer\Documents\fxfiles','AUDUSD_UTC_1 Min_Bid_2007.01.01_2014.08.22.csv') ;

opens = retarr{1}{3}' ; 
clear totalps allprofs
hardstop = .9 ; 
sc = 1 ; 
for scount=20:3:300
    wcount =1 ;
    for w=20:3:300
        
        clear mps
        for xcount=1:30

            clear allps ;
            clear trades                          
%            scount = 20+58*3; 
 %           w = 20+73*3 ;                
          %  ListenChar(2);
            buy = false ; 
            sell = false ; 
            profit = 0 ;
            profits = [0] ; 
            uplarr = [0] ;
            entryarr = [0] ; 
            st = 0 ; 
            entry = 0 ; 
            pcount = 1;  
            i = floor(rand*(size(opens,2)/2))+1000 ;
            istart = i ; 
            repcount = 0 ;
            up = false ;
            down = false ;

            while pcount < 15
                entryarr = opens(i-scount:i)-opens(i-(scount+1)) ;
  %{              
                if 1
                
                  
                  subplot(2,2,1) ; 
                   plot(opens(i-150:i)) ; 
                   
                   subplot(2,2,3) ; 
                   bar((profits)) ; title(num2str(sum(profits))) ; %hline(0,'k') ;  
                   subplot(2,2,2) ; 
                   if buy || sell
                       if size(uplarr,2) > 400
                            imagesc(uplarr(size(uplarr,2)-400:size(uplarr,2))) ; 
                       else
                           imagesc(uplarr) ; %hline(0) ; 
                       end
                   end
                   getframe ;  
                   pause(0.001) ;

                   %st = GetSecs ;    
    %}         
             
                   if size(uplarr,2) > w % if the size of the window is large enough to start monitoring upl
                       if uplarr(size(uplarr,2)) == max(uplarr(size(uplarr,2)-w:size(uplarr,2))) % if we are at a local max in the upl count
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
                      if uplarr(size(uplarr,2)) < -hardstop || uplarr(size(uplarr,2)) > hardstop
                          buy = false  ;
                          profits(pcount) = profit + uplarr(size(uplarr,2)) ; 
                          uplarr = [0] ; 
                          pcount = pcount + 1 ; 
                      end
                    elseif sell == true
                      uplarr(size(uplarr,2)+1) = (opens(entry)-opens(i)) ;   
                      if uplarr(size(uplarr,2)) < -hardstop || uplarr(size(uplarr,2)) > hardstop
                          sell = false ;
                          profits(pcount) = profit + uplarr(size(uplarr,2)) ; 
                          uplarr = [0] ; 
                          pcount = pcount + 1 ; 
                      end
                   end
                   i = i+1 ;
      %          end
                    if ~buy && ~sell
                        if entryarr(scount) == max(entryarr)                            
                            sell = true ; entry = i ;
                        elseif entryarr(scount) == min(entryarr)
                            buy = true ; entry = i ;
                        end
                    end                          
            end
            allps = sum(profits) ;
            allsizes = size(profits,2) ;     
            x = allps(1) ;
            for i=1:size(allps,2)
                x(i+1) = allps(i) + x(i) ; 
            end
            mps(xcount) = mean(profits) ;
            allprofs(xcount,:) = profits ; 
        end       
            totalps(sc,wcount) = mean(mps) ;
            wcount = wcount + 1; 
            %disp(['w = ',num2str(wcount)]) ;   
    end
    disp(['scount = ',num2str(scount)]) ; 
    sc = sc + 1 ; 
end
subplot(1,2,1) ; 
imagesc(totalps)
subplot(1,2,2) ; 
imagesc(imfilter(totalps,fspecial('gaussian',5))) ; 


% replicate in java and then go to demo
% demo for 2 months + then go live