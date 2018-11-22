function newgrad = final_bcg(grad) 
% remove bcg

grad.data(32,:) = rand(1,size(grad.data,2))*0.001 ; 
% epoch the BCG
bdat = grad.data ; bdat = eegfiltfft(bdat,250,0.5,125) ;
nbdat = eegfiltfft(grad.data,250,3,125) ; 
subdat = zeros(size(bdat)) ; 
[weights,sphere] = runica(bdat,'maxsteps',128) ; 
acts = weights*sphere*bdat ; 

clear wskews ; 
wsize = 250*20 ; 
for i=1:5 ; jcount = 1 ; 
    for j=1:wsize:size(acts,2)-wsize
        wskews(i,jcount) = skewness(acts(i,j:j+wsize)) ; 
        jcount = jcount + 1 ; 
    end
end
skews = median((wskews),2) ; 
maxind = find(abs(skews(1:5))==max(abs(skews(1:5)))) ; 
if skews(maxind)<0 ; polarity = -1 ; else polarity = 1 ; end
template = acts(maxind,:)*polarity ; figure,plot(template) ; 



[pks,locs] = findpeaks(template,'MINPEAKDISTANCE',150) ; 

for i=2:length(locs)-2 ;
    allepochs(:,i,:) =  nbdat(:,locs(i)-100:locs(i)+150) ; 
    subepochs(:,i,:) =  bdat(:,locs(i)-100:locs(i)+150) ; 
    inds(i,:) = locs(i)-100:locs(i)+150 ; 
end

for i=1:size(allepochs,2) ; disp(i) ; 
    repi = repmat(allepochs(:,i,:),[1,size(allepochs,2),1]) ; 
    sdiffs = sum(sum(sqrt((allepochs - repi).^2),1),3) ; 
    [sv,si(i,:)] = sort(sdiffs,'ascend') ; 
end

ts = bdat ; 
for i=2:size(si,2)
    x = squeeze(mean(subepochs(:,si(i,2:28),:),2)) ; 
    ts(:,inds(i,:)) = ts(:,inds(i,:)) - x ; 
end
newgrad = grad ; 
ts = ts + eegfiltfft(grad.data,grad.srate,0,0.5) ; 
newgrad.data = ts ; 

end