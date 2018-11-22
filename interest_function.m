clear all ; close all 
houseprice = 200000; 
downpayment = 0.05:0.05:0.5; 
rent = 1000:100:2400; 

for dp=1:length(downpayment)
for rval=1:length(rent)
principal = houseprice - downpayment(dp)*houseprice; 
r = 4.5/1200; 
n = 25*12; 

mp = principal * ((r*((1+r)^n)) / (((1+r)^n)-1)) ;

net(dp,rval) = rent(rval) - mp - 0.50*rent(rval); 
roi(dp,rval) = (net(dp,rval)*12)/(downpayment(dp)*houseprice);

end
end


imagesc(((rent*12)/houseprice)*100,downpayment*100,roi*100); colormap jet; h = colorbar ; title(h,'ROI (%)'); 
xlabel('gross yearly rent as % of purchase price'); 
ylabel('down payment %'); 


