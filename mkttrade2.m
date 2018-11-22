clear all ; close all ;
cd C:\shared\indices_2 ; ls 
allfiles = dir('*csv') ; 

fid = fopen('USA500IDXUSD_UTC_30 Mins_Bid_2011.09.18_2016.08.20.csv') ; 
data = textscan(fid,'%s %s %f %f %f %f %f','delimiter',' ') ; 
fclose(fid) ; 
d5 = data{5} ; d5([1:100,end]) = [] ; % 17% gains a yr...over the past 3 years
% two moving averages, up and down (long and short term trends) 

% make a moving average:
clear mvg p mean
for p = 1:300  
    for i=p+1:length(d5)
        if i==p+1
            mvg(i,p) = mean(d5(i-p+1:i)) ; 
        else
            mvg(i,p) = mvg(i-1,p) + d5(i)/(p) - d5(i-p)/(p) ;
        end  
    end
end
%plot(mvg(:,1:10:end),'r');  hold on ; plot(d5,'b') ; 

m1s = 1:2:20 ; m2s = 1:2:50 ; 
% 1st round of testing: find the optimal parameters
mvgcount = 1 ; clear allprofits allntrades alldiffs
for mvgl=1:2:20 ; 
    mvg2count = 1 ; 
    for mvgl2=1:2:50
        trade = false ; entry = 0 ; profit = 0 ; ntrades = 0 ; equity = 150000 ; losses = 0 ; entryinds = [] ; exitinds = [] ; diffs = [] ;
        for i=301:length(d5)   
            if d5(i) > mvg(i,mvgl) && d5(i-1) < mvg(i,mvgl) && trade == false % cross up
                entry = d5(i) ; trade = true ; entryinds(length(entryinds)+1) = i ; 
            end   
            if trade==true
                if d5(i) < mvg(i,mvgl2) && d5(i-1) > mvg(i,mvgl2) % cross down
                    equity = equity + (d5(i)-entry)*(equity/d5(i))-20 ; trade = false ; 
                    if (d5(i)-entry) < 0 ; losses = losses + 1 ; end ; 
                    exitinds(length(exitinds)+1) = i ; diffs(length(diffs)+1) = (d5(i)-entry) ; 
                    ntrades = ntrades + 1 ; 
                elseif d5(i) - entry < -50
                     equity = equity + (d5(i)-entry)*(equity/d5(i))-20 ; trade = false ; 
                     exitinds(length(exitinds)+1) = i ; diffs(length(diffs)+1) = (d5(i)-entry) ; 
                     ntrades = ntrades + 1 ; 
                end               
            end
        end
        allprofits(mvgcount,mvg2count) = equity ; allntrades(mvgcount,mvg2count) = ntrades ; allnlosses(mvgcount,mvg2count) = losses ; allexitinds{mvgcount,mvg2count} = exitinds ; 
        allentryinds{mvgcount,mvg2count} = entryinds ; alldiffs{mvgcount,mvg2count} = diffs ; 
        mvg2count = mvg2count + 1 ; 
    end ; mvgcount = mvgcount + 1 ; 
end
imagesc(allprofits) ; 

plot(cumsum(alldiffs{1,2})) ; xlabel('trades') ; ylabel('cumulative profit') ; 


% 2nd round: test the parameters on different market sections (small time
% periods)  (visualize first) 
%{
x =2 ; y =2; 
plot(d5) ; hold on ; 
plot(mvg(:,x*2),'m') ; plot(mvg(:,y*2),'k') ; 
for i=1:length(allexitinds{x,y}) ; vline(allentryinds{x,y}(i),'b') ; vline(allexitinds{x,y}(i),'r') ; end
%}

% test it on multiple time points (Try other time frames later (hr, 4hr,
for perm=1:1000
trade = false ; entry = 0 ; profit = 0 ; ntrades = 0 ; losses = 0 ; entryinds = [] ; exitinds = [] ; diffs = [] ; equity = 5000000 ; 
x = 1 ; y = 2 ; mvgl = m1s(x) ; mvgl2 = m2s(y) ; 
randind = round(100 + (rand*(length(d5)-15000))) ; randinds = randind:randind+14899 ; 
for i=randinds(1):randinds(end)  
    if d5(i) > mvg(i,mvgl) && d5(i-1) < mvg(i,mvgl) && trade == false % cross up
        entry = d5(i) ; trade = true ; entryinds(length(entryinds)+1) = i ; 
    end   
    if trade==true
        if d5(i) < mvg(i,mvgl2) && d5(i-1) > mvg(i,mvgl2) % cross down
            equity = equity + (d5(i)-entry)*(equity/d5(i))-20 ; trade = false ; 
            if (d5(i)-entry) < 0 ; losses = losses + 1 ; end ; 
            exitinds(length(exitinds)+1) = i ; diffs(length(diffs)+1) = (d5(i)-entry) ; 
            ntrades = ntrades + 1 ; 
        elseif d5(i) - entry < -50
             equity = equity + (d5(i)-entry)*(equity/d5(i))-20 ; trade = false ; 
             exitinds(length(exitinds)+1) = i ; diffs(length(diffs)+1) = (d5(i)-entry) ; 
             ntrades = ntrades + 1 ; 
        end               
    end
end
allperms(perm) = equity ; pdiffs(perm) = d5(randinds(1))/d5(randinds(end)) ; 
end

[sv,si] = sort(allperms,'descend') ; 




