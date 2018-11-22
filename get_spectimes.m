clear all ; close all ; 
currs = {'AUDJPY','AUDUSD','CHFJPY','EURAUD','EURCAD','EURCHF','EURGBP','EURJPY','EURUSD','GBPCHF','GBPJPY','GBPUSD','NZDUSD','USDCAD','USDCHF','USDJPY'} ; 
pipsz = [0.01,0.0001,0.01,0.0001,0.0001,0.0001,0.0001,0.01,0.0001,0.0001,0.01,0.0001,0.0001,0.0001,0.0001,0.01] ;
rawtps = 1:5:50 ; 
curr = 'USDCAD' ; 
cd('C:\shared\mkt\Forex') ; 
fid = fopen([curr,'.txt']) ; 
data = textscan(fid,'%s %s %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
fclose(fid) ; 

dtimes = data{1} ; 

for i=1:length(dtimes) ; 
    yrs(i) = str2num(dtimes{i}(end-3:end)) ; 
end





