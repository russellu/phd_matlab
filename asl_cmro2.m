
alpha = 0.2 ;  
beta = 1.3 ; 
M=0.1 ; 
%%%%% get the cmro2 time course
mcount = 1 ; clear cmro2t
btimes = boldt ; atimes = 3:5 ; 
for M=0.01:.001:0.2
for i=1:6 ;
    masl = squeeze(mean(aslallmtrials(i,:,:,:),3)) ; masl = imresize(masl,[3,29]) ; 
    mbold = squeeze(mean(boldallmtrials(i,:,:,:),3)) ; 
    for j=1:size(mbold,2)
        boldsi = squeeze(mean(mbold(:,j),3)) ; 
        cbfsi = squeeze(mean(masl(:,j),3)) ;    
        cmro2t(mcount,i,:,j) = ((1-boldsi./M).^(1/beta)).*(cbfsi.^(1-alpha/beta)) ; 
    end
end
mcount = mcount + 1 ; 
end

i=10 ; 
for i=1:20  ;
    subplot(4,5,i) ;
    errorbar(squeeze(mean(real(cmro2t(i,:,[1,2,3],:)),2))',squeeze(std(real(cmro2t(i,:,[1,2,3],:)),0,2))'./sqrt(6),'LineWidth',2) ;
    set(gca,'XTick',1:4:29,'XTickLabel',-16:8:40) ;  title(['M = ',num2str(0.01+0.01*i)]) ; hline(1,'k') ; vline(9,'k') ; legend({'unperturbed','5%contrast','60%random'}) ; 
    xlim([0,30]) ; ylabel('fractional change') ; xlabel('time(s)') ; vline(21,'k') ; 
end
suptitle('CMRO2 for various m-values') ; 











