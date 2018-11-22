[x,y] = meshgrid(-10:.1:10,-10:.1:10) ; 
% grating color vs non color for cortical activation
icount = 1 ; 
colormap gray
for i=.1:.1:1
    subplot(10,1,icount) ; 
    [x,y] = meshgrid(-10:i:300,-10:i:300) ; 
    grate = sin(x) ; 
    grate = grate(1:50,1:50) ; 
    imagesc(grate) ; set(gca,'XTickLabel','') ; set(gca,'YTickLabel','') ; 
    icount = icount + 1 ; 
end