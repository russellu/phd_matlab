clear all ;close all ;
cd c:/shared/MRE ; ls
t1 = load_untouch_nii('RF_Test_MRE_2017_05_25_2_WIP_3D_T1_0.8mm_SENSE_2_1.nii');
vc = load_untouch_nii('rf_vc4_in_t1.nii.gz'); 
vc.img(1,1,:) = -11; vc.img(1,2,:) = 38; 
mre = load_untouch_nii('rf_procpval_ImagBulk_in_t1.nii.gz'); 
mkdir imgs
cd imgs; 
for i=70:100
    %subplottight(5,6,i-69); 
    f = figure,
    plotoverlayIntensity2D(flipud(mre.img(:,:,i)),flipud(mat2gray(vc.img(:,:,i))),flipud(vc.img(:,:,i)),270);
    saveas(f,['img_',num2str(i),'.png']);
end

clear imgs
pngs = dir('*png'); 
for i=1:length(pngs)
   imgs(:,:,:,i) = imread(pngs(i).name);  
    
end
filename = 'giff.gif'; 
for i=1:31
    %subplot(1,2,1); 
    imshow(imgs(:,:,:,i));
    %subplot(1,2,2); 
    %imshow(flipud(imrotate(vc.img(:,:,69+i),90))); set(gca,'XTick',[],'YTick',[]); colormap jet; 
    frame = getframe ;
     im = frame2im(frame); 
           [imind,cm] = rgb2ind(im,256); 

       % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
     
     
    pause(0.01); 
end
