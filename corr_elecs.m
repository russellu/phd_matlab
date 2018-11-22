cd E:\coupling_all_data\resutes\russell


ls

mnicolors = load_untouch_nii('mni_color_in_t1.nii.gz'); 
t1colors = load_untouch_nii('dilsegcoords.nii.gz'); 

for i=1:65
   [cx(i),cy(i),cz(i)] = centmass3(mnicolors.img==i);  
   [tx(i),ty(i),tz(i)] = centmass3(t1colors.img==i);  

end




plot(cx,tx,'r.') ; hold on; 
plot(cy,ty,'g.') ; hold on; 
plot(cz,tz,'b.') ; hold on; 
xlabel('mni-transformed coords (mm)');
ylabel('UTE scan coords (mm)'); 
legend({'x','y','z'});

[c,p] = corr(tx',cx'); 
title(['rho=',num2str(c)]); 
