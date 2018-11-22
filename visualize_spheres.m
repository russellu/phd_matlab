clear all ; close all; 
cd E:\allfmris\sub_alex ; 

t1 = load_untouch_nii('fs_t1.nii.gz'); 


coord = [100,128,128]; 
radius = 10; 

[xg,yg,zg] = ndgrid(-coord(1):size(t1.img,1)-coord(1)-1,-coord(2):size(t1.img,2)-coord(2)-1,-coord(3):size(t1.img,3)-coord(3)-1); 

sphere = sqrt(xg.^2 + yg.^2 + zg.^2) < radius; 

%disp3d(double(t1.img).*double(sphere));

t1.img = double(t1.img) + double(sphere*255); 
save_untouch_nii(t1,'lyes_sphere.nii.gz');
