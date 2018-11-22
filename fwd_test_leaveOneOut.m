currs =   {'AUDCAD','AUDCHF','AUDJPY','AUDNZD','AUDSGD','AUDUSD','CADCHF','CADHKD','CADJPY','CHFJPY','CHFSGD','EURAUD','EURCAD','EURCHF','EURDKK','EURGBP','EURHKD','EURJPY','EURNZD','EURPLN','EURSEK','EURSGD','EURTRY','EURUSD','GBPAUD','GBPCAD','GBPCHF','GBPJPY','GBPNZD','GBPUSD','HKDJPY','NZDCAD','NZDCHF','NZDJPY','NZDUSD','SGDJPY','USDCAD','USDCHF','USDCNH','USDDKK','USDHKD','USDJPY','USDMXN','USDNOK','USDSGD','USDTRY','USDZAR','XAGUSD','XAUUSD','ZARJPY'} ; 
spreads = [ 2.5,     2.5,     0.95,    2.33,    5,       0.9,     1.95,    16.0     1.05,    1.5,     9,       1.9,     1.85,    1.0      4.0      0.85,    18.0,    0.6,     3.6,     20.0,    19.0,    5.25,    9.0,     0.25,    3.35,    3.5,     2.2,     1.6,     5.0,     0.96,    28.0,    3.0,     2.1,     2.1,     1.05,    2.5,     0.9,     1.05,    2.5,     4.9,     3.0,     0.35,    50.0,     22.0,    2.0,     6.0,     100.0,   3.05,    28.0,    0.7   ]; 
[spreadsort,spreadsortind] = sort(spreads,'ascend') ; 
currs = currs(spreadsortind) ; spreads = spreads(spreadsortind) ; 
cd C:\Users\Acer\Documents\back_test ; 
fparams = dir('finalparams*') ; 
for f=1:length(currs) 
    %{
    cd C:\Users\Acer\Documents\back_test ; 
    corrparams = load(fparams(f).name) ; corrparams = corrparams.finalparams ; 
    mcorrs = squeeze(mean(corrparams,5)) ; 
    allmcorrs(:,:,:,:,f) = mcorrs ; 
   %}
    cd C:\Users\Acer\Documents\fwd_test ; 
 %   [sv,si] = sort(mcorrs(:),'descend') ; 
 %   [i1,i2,i3,i4] = ind2sub(size(mcorrs),si) ; 
 %   cstring = strrep(fparams(f).name,'finalparams_chunk_','') ; cstring = strrep(cstring,'.mat','') ; 
    cstring = currs{f} ; 
    ccy=find(strcmpi(cstring,currs)) ; 
    params = [3,19,19,1] ;% allparams(f,:) = params ; 
    ccy_current = dir([currs{ccy},'*.csv']) ; 
    if (~isempty(strfind(currs{f},'JPY')) && isempty(strfind(currs{f},'HKD')) ) || ~isempty(strfind(currs{f},'XAG')) || ~isempty(strfind(currs{f},'XAU')) 
        lim1 = -0.4 ; lim2 = 0.4 ; mfactor = 0.01 ; 
    else 
        lim1 = -0.004 ; lim2 = 0.004 ; mfactor = 0.0001 ; 
    end
    totalcost = (spreads(ccy) + 0.35 + 0.9)*mfactor ; 
    fid = fopen(ccy_current.name) ; 
    dkdata = textscan(fid,'%s %s %f %f %f %f %f','delimiter',' ') ; 
    fclose(fid) ; 
    for idd=1:59 ; 
    d5 = dkdata{6} ; d5 = d5(idd:60:end) ; 
    mvgs1 = (1:2:40).^2 ; 
    pipthreshs = [0,30]*mfactor ; 
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

    m1 = params(1) ; m2 = params(2) ; m3 = params(3) ; pipthresh = params(4) ; 
    trade = false ; entry = 0 ; buy = false ; sell = false ; profit = 0 ; ntrades = 0 ; 
    alltrades = [] ; sellentries = [] ; buyentries = [] ; sellexits = [] ; buyexits = [] ;trackprofits = [] ; 
    for i=1600:length(d5) 
        if mod(i,100)== 0 ; trackprofits(length(trackprofits)+1) = profit ; end
        if trade == false
            if d5(i)-pipthreshs(pipthresh) < mvg(i,(m1)) % buy mode
                if d5(i) < mvg(i,(m2)) && d5(i-1) > mvg(i,(m2)) % cross down
                    trade = true ; entry = d5(i) ; buy = true ; ntrades = ntrades + 1  ;
                end
            elseif d5(i)+ pipthreshs(pipthresh) > mvg(i,(m1)) % sell mode
                if d5(i) > mvg(i,(m2)) && d5(i-1) < mvg(i,(m2)) % cross up
                    trade = true ; entry = d5(i) ; sell = true ; ntrades = ntrades + 1 ;  
                end
            end
        elseif trade == true
            if buy
                if d5(i) > mvg(i,(m3)) && d5(i-1) < mvg(i,(m3))
                    buy = false ; profit = profit + d5(i)-entry-totalcost ; trade = false ; alltrades(length(alltrades)+1) = d5(i) - entry - totalcost ; 
                end
            elseif sell
                if d5(i) < mvg(i,(m3)) && d5(i-1) > mvg(i,(m3))
                    sell = false ; profit = profit + entry-d5(i)-totalcost ; trade = false ; alltrades(length(alltrades)+1) = entry - d5(i) -totalcost ; 
                end
            end
        end
    end
    iddprofits(idd) = profit ; iddntrades(idd) = ntrades ;
    end
    figure,bar(iddprofits) ; title(['ccy=',currs{ccy},', meanprofit=',num2str(mean(iddprofits)),', meanntrades=',num2str(mean(iddntrades)),', meanppt=',num2str(mean(iddprofits)./mean(iddntrades))]) ; 
    allprofs(f,:) = iddprofits./mfactor ; allntrades(f,:)= iddntrades ; 
    
end


avgt = allprofs./allntrades ; 
