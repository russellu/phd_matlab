clear all ; close all ; 
cd c:/shared/indices_4 ; ls 
alldatas = dir('datas_*') ; 

    
fid = fopen('EURPLN_UTC_1 Min_Bid_2011.12.31_2016.09.01.csv') ; 
dkdata = textscan(fid,'%s %s %f %f %f %f %f','delimiter',' ') ; 
fclose(fid) ; 
d5 = dkdata{6} ; d5 = d5(1:30:end) ; 



for ad=1:length(alldatas) 
    dat = load(alldatas(ad).name) ;  
    dat = dat.datas ; 
    profits = dat{1} ; 
    ntrades = dat{2} ; 
    tracks = dat{3} ; 
    allprofits(:,:,:,ad) = profits ; allntrades(:,:,:,ad) = ntrades ; alltracks{ad} = tracks ; 
    divp = profits./ntrades ; divp(isnan(divp)) = 0 ; divp(isinf(divp)) = 0 ; 
    mdiv = squeeze(mean(mean(divp,4),5)) ; 
    lim1 = 0 ; lim2 = 0 ; 
    if ~isempty(strfind(alldatas(ad).name,'JPY'))
        lim1 = -0.15 ; lim2 = 0.15 ; 
    else 
        lim1 = -0.0015 ; lim2 = 0.0015 ; 
    end
    figure,
    for i=1:28 ; subplot(4,7,i) ;imagesc(squeeze(mdiv(:,:,i)),[lim1,lim2]) ; title(i) ; end
    suptitle(alldatas(ad).name) ;    
    
    for i=1:size(tracks,1) ; 
        for j=1:size(tracks,2)
            for k=1:size(tracks,3)          
                if ~isempty(tracks{i,j,k})
                    corrs(i,j,k,ad) = corr2(cumsum(tracks{i,j,k}),1:length(cumsum(tracks{i,j,k}))) ; 
                end
            end
        end
    end

end


for i=1:3 ; figure ; for j=1:28 ; subplot(4,7,j) ; imagesc(squeeze(corrs(:,:,j,i)),[.95,1]) ; title(j) ; end ; suptitle(alldatas(i).name(7:12)) ; end
a = [25,20,21,1] ;
figure,plot(cumsum(alltracks{a(4)}{a(1),a(2),a(3)}))













mvg = 1:25:700 ; 
[mvg(26),mvg(20),mvg(21)]