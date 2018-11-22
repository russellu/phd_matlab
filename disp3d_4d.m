function disp3d_4d(input4d)
% display a 4d image as multiple subplots in time
% right now assumes it is a 3d+t or 3d+weight image, and scrolls in the 
% z-direction
for i=1:size(input4d,4) ; 
    mins(i) = min(min(min(squeeze(input4d(:,:,:,i))))) ; 
    maxes(i) = max(max(max(squeeze(input4d(:,:,:,i))))) ; 
end

squaresz = ceil(sqrt(size(input4d,4))) ; 
for i=1:100
    for j=1:size(input4d,3)
        for kcount=1:size(input4d,4) 
            subplot(squaresz,squaresz,kcount) ; 
            imagesc(squeeze(input4d(:,:,j,kcount))) ; 
            set(gca,'XTickLabel',[],'YTickLabel',[]) ; 
            title(num2str(kcount)) ; 
        end
        getframe ;
    end
end

end