clear all  ; close all ;
pcount = 1 ;
for popsize=2:200
    pcount
pop = rand(popsize,1) ; 
m = mean(pop) ;
clear jacks
for i=1:size(pop,1)
    newpop = pop ;
    newpop(i) = [] ; 
    %mean(newpop)
    %jacks(i) = m*size(pop,1) - (size(pop,1)-1).*mean(newpop) ; 
    jacks(i) = mean(newpop) ; 
end
imeans(pcount) = m ; 
jmeans(pcount) = mean(jacks) ;
pcount = pcount + 1 ; 
end
subplot(1,2,1) ; plot(imeans-.5) ; hline([-.2,.2],'k') ; hline(0) ; title('mean(sample) - .5') ;   hold on ; 
%subplot(1,2,2) ;
plot(jmeans-.5,'r') ; hline([-.2,.2],'k') ; hline(0) ;  title('mean(jacknife(sample)) - .5') ;
suptitle('effect of jacknife of estimate of mean for pseudo random number generation [0,1]') ; 

