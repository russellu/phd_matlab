clear all ; close all; 

cd(['e:\dukas\forex\minute_0001']) 
namei = dir(['EURUSD','*','Bid*']); disp(namei.name); 
fid = fopen(namei.name) ; 
dkbid = textscan(fid,'%s %f %f %f %f %f %f','delimiter',',','Headerlines',1) ; 
fclose(fid) ; 

d5 = dkbid{5}; 

a = (d5(1320600:1322000)); 
b = cumsum(rand(1,1401)-0.5); 

mvg_a = getmvg(a,150);
mvg_b = getmvg(b,150); 

subplot(1,2,1) ; 
plot(a(151:end)) ; hold on ; plot(mvg_a(151:end),'LineWidth',2); set(gca,'XTickLabel',[],'YTickLabel',[]); xlabel('time'); title('plot #1');
subplot(1,2,2) ;
plot(b(151:end)); hold on ; plot(mvg_b(151:end),'LineWidth',2); set(gca,'XTickLabel',[],'YTickLabel',[]); xlabel('time'); title('plot #2'); 
