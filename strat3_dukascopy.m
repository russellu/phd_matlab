clear all ; close all ; 
%{
currs = {'AUDJPY','AUDUSD','CHFJPY','EURAUD','EURCAD','EURCHF','EURGBP','EURJPY','EURUSD','GBPCHF','GBPJPY','GBPUSD','NZDUSD','USDCAD','USDCHF','USDJPY'} ; 
curr = 'EURJPY' ; 
cd('C:\Users\Acer\Downloads\market_081916_russell_butler\Market 081916\Forex') ; 
fid = fopen([curr,'.txt']) ; 
data = textscan(fid,'%s %s %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
fclose(fid) ; 
d5 = data{6} ; d5 = d5((end/1.5+5):5:end) ;
%}

cd c:/shared/indices_4 ; ls 
currs = dir('*csv') ; 

clear allprofits alltracks allntrades
for curr=33:length(currs) 
    
fid = fopen(currs(curr).name) ; 
dkdata = textscan(fid,'%s %s %f %f %f %f %f','delimiter',' ') ; 
fclose(fid) ; 
d5 = dkdata{6} ; d5 = d5(1:5:end) ; 

mvgs1 = 1:50:2000 ; 
mvgs2 = 1:50:2000 ; 
% make a moving average:
clear mvg p mean
for p = 1:length(mvgs1)  
    for i=mvgs1(p)+1:length(d5)
        if i==mvgs1(p)+1
            mvg(i,p) = mean(d5(i-mvgs1(p)+1:i)) ; 
        else
            mvg(i,p) = mvg(i-1,p) + d5(i)/(mvgs1(p)) - d5(i-mvgs1(p))/(mvgs1(p)) ;
        end  
    end
end

for m1=1:length(mvgs1) ; disp(['m1=',num2str(m1)]) ;
    for m2=1:length(mvgs2) ; 
        for m3=1:length(mvgs2) ; 
            trade = false ; entry = 0 ; buy = false ; sell = false ; profit = 0 ; ntrades = 0 ; alltrades = [] ; sellentries = [] ; buyentries = [] ; sellexits = [] ; buyexits = [] ; 
            for i=700:length(d5)
                if trade == false
                    if d5(i) > mvg(i,(m1)) % buy mode
                        if d5(i) < mvg(i,(m2)) && d5(i-1) > mvg(i,(m2)) % cross down
                            trade = true ; entry = d5(i) ; buy = true ; ntrades = ntrades + 1  ;
                        end
                    elseif d5(i) < mvg(i,(m1)) % sell mode
                        if d5(i) > mvg(i,(m2)) && d5(i-1) < mvg(i,(m2)) % cross up
                            trade = true ; entry = d5(i) ; sell = true ; ntrades = ntrades + 1 ;  
                        end
                    end
                elseif trade == true
                    if buy
                        if d5(i) > mvg(i,(m3)) && d5(i-1) < mvg(i,(m3))
                            buy = false ; profit = profit + d5(i)-entry ; trade = false ; alltrades(length(alltrades)+1) = d5(i) - entry ; 
                        end
                    elseif sell
                        if d5(i) < mvg(i,(m3)) && d5(i-1) > mvg(i,(m3))
                            sell = false ; profit = profit + entry-d5(i) ; trade = false ; alltrades(length(alltrades)+1) = entry - d5(i) ; 
                        end                   
                    end
                end    
            end  
            allprofits(m1,m2,m3) = profit ; 
            allntrades(m1,m2,m3) = ntrades ; 
            alltracks{m1,m2,m3} = alltrades ;
        end
    end  
end

datas{1} = allprofits ; datas{2} = allntrades ; datas{3} = alltracks ; 
save(['5m_datas_',currs(curr).name,'.mat'],'datas') ; 
end



