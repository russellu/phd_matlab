clear all ; close all ; 
currs =   {'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDSGD','AUDUSD','CADCHF','CADHKD','CADJPY','CHFJPY',...
    'CHFSGD','EURAUD','EURCAD','EURCHF','EURDKK','EURGBP','EURHKD','EURJPY','EURNZD','EURPLN','EURSEK',...
    'EURSGD','EURTRY','EURUSD','GBPAUD','GBPCAD','GBPCHF','GBPJPY','GBPNZD','GBPUSD','HKDJPY','NZDCAD','NZDCHF',...
    'NZDJPY','NZDUSD','SGDJPY','USDCAD','USDCHF','USDCNH','USDDKK','USDHKD','USDJPY','USDMXN','USDNOK','USDSGD','USDTRY','USDZAR','XAGUSD','XAUUSD','ZARJPY'} ; 
spreads = [ 2.5,     2.5,     0.95,    2.33,    5,       0.9,     1.95,    16.0     1.05,    1.5,     9,       1.9,     1.85,    1.0      4.0      0.85,    18.0,    0.6,     3.6,     20.0,    19.0,    5.25,    9.0,     0.25,    3.35,    3.5,     2.2,     1.6,     5.0,     0.96,    28.0,    3.0,     2.1,     2.1,     1.05,    2.5,     0.9,     1.05,    2.5,     4.9,     3.0,     0.35,    50.0,     22.0,    2.0,     6.0,     100.0,   3.05,    28.0,    0.7   ]; 
cd C:\shared\indices_4 ; ls 
bads = {'CADHKD','EURDKK','EURHKD','EURPLN','EURSEK','EURTRY','HKDJPY','USDCNH','USDDKK','USDHKD','USDMXN','USDNOK','USDTRY','USDZAR','XAUUSD','ZARJPY','GBPAUD'} ;
resnames = dir('results_*') ; badinds = [] ; 
for i=1:length(resnames) ; badfound = false ; 
    for j=1:length(bads)
        if ~isempty(strfind(resnames(i).name,bads{j}))
            badfound = true ; 
        end
    end ; if badfound ; badinds(length(badinds)+1) = i ; end
end
resnames(badinds) = []  ; 

for ccy=1:length(resnames) ;
    ccy_current = resnames(ccy).name ;  
    if (~isempty(strfind(ccy_current,'JPY')) && isempty(strfind(ccy_current,'HKD')) ) || ~isempty(strfind(ccy_current,'XAG')) || ~isempty(strfind(ccy_current,'XAU')) 
        lim1 = -0.4 ; lim2 = 0.4 ; mfactor = 0.01 ; 
    else 
        lim1 = -0.004 ; lim2 = 0.004 ; mfactor = 0.0001 ; 
    end
    results = load(resnames(ccy).name) ; results = results.results ; 
    allprofits = results{1} ; allntrades = results{2} ; alltracks = results{3} ; alltrackprofits = results{4} ; 
    for i=1:size(alltrackprofits,1);for j=1:size(alltrackprofits,2);for k=1:size(alltrackprofits,3);for el=1:size(alltrackprofits,4);corrs(i,j,k,el)=corr2(squeeze(alltrackprofits(i,j,k,el,:))',1:size(alltrackprofits,5));end;end;end;end
    corrs(isnan(corrs)) = 0 ; 
    [sv,si] = sort(corrs(:),'descend') ; 
    [ix,iy,iz,it] = ind2sub(size(corrs),si) ; 
    for i=1:2 ;
        hps(i,:) = imresize(squeeze(alltrackprofits(ix(i),iy(i),iz(i),it(i),:)),[600,1]) ; 
    end
    %figure,plot(mean(hps))
    allhps(ccy,:) = hps(1,:)./mfactor ; 
end

plot(mean(allhps))
