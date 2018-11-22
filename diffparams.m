clear all ; close all ; 
cd('C:\shared\mkt') ; 
fid = fopen('ES.txt') ; 
data = textscan(fid,'%s %s %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
fclose(fid) ; 
d5 = data{6} ; d5 = d5(1:5:end) ; 
d5 = d5(length(d5)/2:end) ; 
%d5 = d5(end-50000:end) ; 

% make a moving average:
clear mvg p mean
for p = 1:500  
    for i=p+1:length(d5)
        if i==p+1
            mvg(i,p) = mean(d5(i-p+1:i)) ; 
        else
            mvg(i,p) = mvg(i-1,p) + d5(i)/(p) - d5(i-p)/(p) ;
        end  
    end
end

diffmvg = diff(mvg) ; 

clear allprofits allntrades alltracks
m1s = 1:20:500 ; m2s = 1:20:500 ; 
m1count = 1 ; 
for m1=1:5:300 ; m2count = 1 ; disp(m1) ; 
    for m2=1:5:300 ; diffcount = 1 ; 
        for diffm = 0:0.05:0.6
            trade = false ; entry = 0 ; buy = false ; sell = false ; profit = 0 ; ntrades = 0 ; alltrades = [] ; sellentries = [] ; buyentries = [] ; sellexits = [] ; buyexits = [] ; 
            for i=201:length(d5)-1
                if trade == false
                   % if d5(i) > mvg(i,m1) % buy mode
                        if diffmvg(i-1,m2)<-diffm %d5(i) < mvg(i,m2) && d5(i-1) > mvg(i,m2) % cross down
                            trade = true ; entry = d5(i) ; buy = true ; ntrades = ntrades + 1  ;
                     %   end
                        elseif diffmvg(i-1,m2)>diffm % d5(i) < mvg(i,m1) % sell mode
                   %     if diffmvg(i-1,m1)>diffm %d5(i) > mvg(i,m2) && d5(i-1) < mvg(i,m2) % cross up
                            trade = true ; entry = d5(i) ; sell = true ; ntrades = ntrades + 1 ; 
                    %    end
                        end
                elseif trade == true
                    if buy
                        if diffmvg(i-1,m1)<-diffm 
                            buy = false ; profit = profit + d5(i)-entry ; trade = false ; alltrades(length(alltrades)+1) = d5(i) - entry ; 
                        end
                    elseif sell
                        if  diffmvg(i-1,m1)>diffm 
                            sell = false ; profit = profit + entry-d5(i) ; trade = false ; alltrades(length(alltrades)+1) = entry - d5(i) ; 
                        end                   
                    end
                end    
            end ; 

            allprofits(m1count,m2count,diffcount) = profit ; allntrades(m1count,m2count,diffcount) = ntrades ; alltracks{m1count,m2count,diffcount} = alltrades ;
            allbuyentries{m1count,m2count,diffcount} = buyentries ; allsellentries{m1count,m2count,diffcount} = sellentries ; allbuyexits{m1count,m2count,diffcount} = buyexits ; allsellexits{m1count,m2count,diffcount} = sellexits ; 
            diffcount = diffcount + 1 ; 
        end ;m2count = m2count + 1 ;
    end ; m1count = m1count + 1 ; 
end ;
figure,for i=1:13 ; subplot(4,4,i); imagesc(squeeze(allprofits(:,:,i)./allntrades(:,:,i)),[-2,2]) ; end

figure,plot(cumsum(alltracks{34,1,4}-.25))

