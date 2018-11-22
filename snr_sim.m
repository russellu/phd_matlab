clear all; close all; 
% keep keep noise constant, vary amplitude and distance
ampls = 0.1:0.1:10;
dists = 1:.2:20; 
allas = zeros(length(ampls),length(dists),9991); 
for i=1:length(ampls)
    for j=1:length(dists)
        boxcar = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]; 
        a = (sin(1:0.1:1000)*ampls(i))./dists(j).^2; 
        rands = rand(1,length(a)); 
        boxcar = imresize(boxcar,[1,length(rands)],'nearest'); 
        a = a.*boxcar;
        a = a+rands; 
        allas(i,j,:) = a; 
        tasks = a(boxcar>0); 
        rests = a(boxcar==0); 
        nelems = min([length(tasks),length(rests)]); 
        ftask = log(abs(fft(tasks(1:nelems))));
        frest = log(abs(fft(rests(1:nelems)))); 
        inds = 60:90; 
        fdiff = ftask-frest; 
        testvals(i,j) = mean(fdiff(inds)); 
    end
end

subplot(2,2,1); 
imagesc(dists,ampls,testvals) ; xlabel('distance'); ylabel('amplitude'); title('SNR as a function of amplitude and distance');
subplot(2,2,2); 
plot(dists,testvals(10,:)) ; hold on ; plot(dists,testvals(end,:)); hline(0,'k'); 
ylabel('snr'); xlabel('distance'); xlim([1,max(dists)]); vline(5,'r'); legend({'gamma','alpha'});title('SNR as a function of distance for two amplitudes');
subplot(2,2,3);
plot(squeeze(allas(10,5,1:1500))) ; hold on ; plot(squeeze(allas(10,15,1:1500))) ; title('gamma at 2mm and 10mm');
legend({'2mm','10mm'}) ; ylim([-2,2]);vline([650,1340],'k'); 
subplot(2,2,4);
plot(squeeze(allas(50,5,1:1500))); hold on ; plot(squeeze(allas(50,15,1:1500))) ; title('alpha at 2mm and 10mm'); 
legend({'2mm','10mm'});ylim([-2,2]); vline([650,1340],'k');



