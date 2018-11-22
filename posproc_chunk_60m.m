clear all ; close all ; 
currs =   {'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDSGD','AUDUSD','CADCHF','CADHKD','CADJPY','CHFJPY','CHFSGD','EURAUD','EURCAD','EURCHF','EURDKK','EURGBP','EURHKD','EURJPY','EURNZD','EURPLN','EURSEK','EURSGD','EURTRY','EURUSD','GBPAUD','GBPCAD','GBPCHF','GBPJPY','GBPNZD','GBPUSD','HKDJPY','NZDCAD','NZDCHF','NZDJPY','NZDUSD','SGDJPY','USDCAD','USDCHF','USDCNH','USDDKK','USDHKD','USDJPY','USDMXN','USDNOK','USDSGD','USDTRY','USDZAR','XAGUSD','XAUUSD','ZARJPY'} ; 
spreads = [ 2.5,     2.5,     0.95,    2.33,    5,       0.9,     1.95,    16.0     1.05,    1.5,     9,       1.9,     1.85,    1.0      4.0      0.85,    18.0,    0.6,     3.6,     20.0,    19.0,    5.25,    9.0,     0.25,    3.35,    3.5,     2.2,     1.6,     5.0,     0.96,    28.0,    3.0,     2.1,     2.1,     1.05,    2.5,     0.9,     1.05,    2.5,     4.9,     3.0,     0.35,    50.0,     22.0,    2.0,     6.0,     100.0,   3.05,    28.0,    0.7   ]; 
[spreadsort,spreadsortind] = sort(spreads,'ascend') ; 
currs = currs(spreadsortind) ; spreads = spreads(spreadsortind) ; 
cd c:/shared/back_test; ls 
chunk60s = dir('chunk_60m*') ; 
for c=1:length(chunk60s)
    chunk = load(chunk60s(c).name) ; chunk = chunk.savestuffs ; 
    profits = chunk{1} ; ntrades = chunk{2} ; tracks = chunk{3} ; 
    for i=1:size(tracks,1)
        for j=1:size(tracks,2) 
            for k=1:size(tracks,3)
                for el=1:size(tracks,5)
                    if ~isempty(tracks{i,j,k,1,el})
                        corrs(i,j,k,el) = corr2(squeeze(cumsum(tracks{i,j,k,1,el})),1:length(cumsum(tracks{i,j,k,1,el}))) ; 
                    end
                end
            end
        end
    end
    corrs(isnan(corrs)) = 0 ; 
    mcorrs = squeeze(mean(corrs,4)) ;  
    figure,for i=1:size(mcorrs,3) ; subplot(4,4,i) ; imagesc(squeeze(mcorrs(:,:,i)),[-1,1]) ; end ; title(chunk60s(c).name) ; 
    allcorrs(:,:,:,c) = mcorrs ; 
    
    
end
mallcorrs = squeeze(mean(allcorrs,4)) ; 

