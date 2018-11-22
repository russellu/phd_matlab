cd c:/users/butr2901/documents/profitcurves30;  ls 
txts=dir('*txt') ; 
for t=1:length(txts)
   txtt = load(txts(t).name) ;  
   figure,plot(txtt') ; title(txts(t).name) ;  
end