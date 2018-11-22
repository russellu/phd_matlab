clear all ; close all ; 
currs =   {'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDSGD','AUDUSD','CADCHF','CADHKD','CADJPY','CHFJPY','CHFSGD','EURAUD','EURCAD','EURCHF','EURDKK','EURGBP','EURHKD','EURJPY','EURNZD','EURPLN','EURSEK','EURSGD','EURTRY','EURUSD','GBPAUD','GBPCAD','GBPCHF','GBPJPY','GBPNZD','GBPUSD','HKDJPY','NZDCAD','NZDCHF','NZDJPY','NZDUSD','SGDJPY','USDCAD','USDCHF','USDCNH','USDDKK','USDHKD','USDJPY','USDMXN','USDNOK','USDSGD','USDTRY','USDZAR','XAGUSD','XAUUSD','ZARJPY'} ; 
spreads = [ 2.5,     2.5,     0.95,    2.33,    5,       0.9,     1.95,    16.0     1.05,    1.5,     9,       1.9,     1.85,    1.0      4.0      0.85,    18.0,    0.6,     3.6,     20.0,    19.0,    5.25,    9.0,     0.25,    3.35,    3.5,     2.2,     1.6,     5.0,     0.96,    28.0,    3.0,     2.1,     2.1,     1.05,    2.5,     0.9,     1.05,    2.5,     4.9,     3.0,     0.35,    50.0,     22.0,    2.0,     6.0,     100.0,   3.05,    28.0,    0.7   ]; 
[spreadsort,spreadsortind] = sort(spreads,'ascend') ; 
currs = currs(spreadsortind) ; spreads = spreads(spreadsortind) ; 
cd c:/users/acer/documents/back_test; ls 
chunks = dir('chunk_*') ; 
for ch=1:length(chunks)
    chunki = load(chunks(ch).name) ; chunki = chunki.savestuffs ; 
    if (~isempty(strfind(chunks(ch).name,'JPY')) && isempty(strfind(chunks(ch).name,'HKD')) ) || ~isempty(strfind(chunks(ch).name,'XAG')) || ~isempty(strfind(chunks(ch).name,'XAU')) 
        lim1 = -0.4 ; lim2 = 0.4 ; mfactor = 0.01 ; 
    else 
        lim1 = -0.004 ; lim2 = 0.004 ; mfactor = 0.0001 ; 
    end
    disp(['ccy = ',chunks(ch).name,' mfactor = ',num2str(mfactor)]) ; 
    profits = chunki{1} ; ntrades = chunki{2} ; tracks = chunki{3} ; 
    corrs = zeros(size(tracks)) ; 
    for i=1:size(tracks,1)
        for j=1:size(tracks,2)
            for k=1:size(tracks,3)
                for el=1:size(tracks,4)
                    for m=1:size(tracks,5)
                        if ~isempty(tracks{i,j,k,el,m})
                            corrs(i,j,k,el,m) = corr2(squeeze(cumsum(tracks{i,j,k,el,m})),1:length(cumsum(tracks{i,j,k,el,m}))) ;                         
                        end
                    end
                end
            end
        end
    end
    corrs(isnan(corrs)) = 0 ; 
    clear medcorrs ; 
    for i=1:size(corrs,4)
        for j=1:size(corrs,5)
            medcorrs(:,:,:,i,j) = medfilt3(squeeze(corrs(:,:,:,i,j))) ; 
        end
    end
    meanmedcorrs = squeeze(mean(medcorrs(:,:,:,1,:),5)) ;     
    sampinds = {[1,2,3,4,5,6,7,8,9],[1,2,3,4,5,6,7,8,10],[1,2,3,4,5,6,7,9,10],[1,2,3,4,5,6,8,9,10],[1,2,3,4,5,7,8,9,10],[1,2,3,4,6,7,8,9,10],...
        [1,2,3,5,6,7,8,9,10],[1,2,4,5,6,7,8,9,10],[1,3,4,5,6,7,8,9,10],[2,3,4,5,6,7,8,9,10]} ; 
    nonsampinds = [10,9,8,7,6,5,4,3,2,1] ; 
    clear outprofs
    for i=1:length(sampinds)
        outcorrs(:,:,:,:,i) = mean(medcorrs(:,:,:,:,sampinds{i}),5) ; 
        outi = squeeze(outcorrs(:,:,:,:,i)) ; outi(isnan(outi)) = 0 ; 
        [sv,si] = sort(outi(:),'descend') ; [i1,i2,i3,i4] = ind2sub(size(outi),si) ; 
        for s=1:5
            outprofs(i,s) = profits(i1(s),i2(s),i3(s),i4(s),nonsampinds(i)) ; 
        end
    end
    alloutprofs(ch,:,:) = outprofs./mfactor ; figure,bar(mean(outprofs,2)) ; title(['mean=',num2str(mean(mean(outprofs))),', ccy=',chunks(ch).name]) ;     
    finalparams = corrs ; 
    save(['finalparams_',chunks(ch).name],'finalparams') ; 
end


