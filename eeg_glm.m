% try also combining electrodes which showed high correlation with
% component weight maps 
clear all ; close all
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','jeremie','julie','katrine','lisa','marc','marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','vincent'} ;    

for ss=1:length(subs) ; 
cd c:/shared/allres ; cd(subs{ss}) ; 
beta = load('c:/shared/papsaves/beta_template') ; beta = beta.beta_template ; 
gamma = load('c:/shared/papsaves/gamma_template') ; gamma = gamma.gamma_template ; 
ica = pop_loadset('ica_notch85.set')  ;
ica11 = pop_epoch(ica,{'S 11','S 15','S 13','S 12','S 16'},[-1,2.5]) ; 
beta = imresize(beta,[1,size(ica11.data,2)]) ; gamma = imresize(gamma,[1,size(ica11.data,2)]) ; 
clear filtact filtactb ; 
for i=1:64 ; 
    filtact(i,:,:) = eegfiltfft(squeeze(ica11.icaact(i,:,:))',ica.srate,35,110) ; 
    filtactb(i,:,:) = eegfiltfft(squeeze(ica11.icaact(i,:,:))',ica.srate,15,25) ; 
end
times = ica11.times ; 
filtact = filtact.^2 ;%- (filtactb.^2) ; 
filtactb = filtactb.^2 ; 
sumact = squeeze(sum(filtact,3)) ; zsums = zscore(sumact,0,2) ;
for i=1:64 ; goods{i} = find(zsums(i,:)<5) ; end
ideal = zeros(1,length(ica11.times)) ; 
ideal(257:768) = 1 ; 
for i=1:64;
    for j=1:size(filtact,2) ;
        %datcorrmat(i,j) = corr2(squeeze(filtdat(i,j,:)),ideal') ; 
        actcorrmat(i,j) = corr2(squeeze(filtact(i,j,:)),gamma') + (corr2(squeeze(filtactb(i,j,:)),beta'))./4 ; 
    end
end
zelecs = abs(zscore(ica.icawinv,1,1)) ; 
maxz = max(zelecs,[],1) ; bads = find(maxz>5) ;  
for i=1:64 ; 
   [ersp(i,:,:),itc,powbase,times,freqs,~,~] = newtimef(squeeze(ica11.icaact(i,:,goods{i})),ica11.pnts,[ica11.xmin,ica11.xmax],ica11.srate,...
                                               0,'plotitc','off','plotersp','off','freqs',[1,120],'nfreqs',60,'winsize',round(ica11.srate/4)) ;    
end
mcorrs = squeeze(mean(actcorrmat,2)) ; mcorrs(bads) = -1 ; 
[sv,si] = sort(mcorrs,'descend') ; 
figure;
for i=1:64 ; 
    subplot(8,8,i) ; 
    imagesc(squeeze(ersp(si(i),:,:)),[-3,3]) ; 
    %title(num2str(mean(actcorrmat(si(i),:),2))) ;
    title(num2str(mcorrs(si(i)))) ; 
    a(i) = std(mean(ersp(si(i),freqs>30 & freqs<120,times>0),2),0,3) ; 
    allersp(ss,i,:,:) = ersp(si(i),:,:) ; 
end
allsis(ss,:) = si ; 
suptitle(subs{ss}) ; 
end



%{
templ = squeeze(mean(mean(allersp(:,1:2,:,:)))) ; 
t = find(times<2.5) ; 
plot(mean(templ(find(freqs>15 & freqs<25),t))) ; hold on ; 
plot(mean(templ(find(freqs>40 & freqs<110),t))*3,'r') ; hline(0,'k') ; 

beta_template = mean(templ(find(freqs>15 & freqs<25),t)) ; 
gamma_template = mean(templ(find(freqs>40 & freqs<110),t)) ; 
%}