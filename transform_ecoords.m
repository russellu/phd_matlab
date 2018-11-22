%clear all ; close all ; 
cd c:/shared/freesurfer_segs/sub1_valerie
cd mri
% get the freesurfer file:
t1 = load_untouch_nii('C:\shared\freesurfer_segs\sub1_alex\mri\T1.mgz.nii.gz') ; 
t1.img(:,25:end,:) = 0 ; 
[nasx,nasy,nasz] = centmass3(t1.img==50) ; 
[lpx,lpy,lpz] = centmass3(t1.img==150) ; 
[rpx,rpy,rpz] = centmass3(t1.img==250) ; 

fids = [nasx,nasy,nasz;lpx,lpy,lpz;rpx,rpy,rpz] ; 

coords = load_untouch_nii('C:\shared\freesurfer_segs\sub1_alex\mri\coords.nii.gz') ; clear elocs
for i=1:65 ; 
    [elocs(i,1),elocs(i,2),elocs(i,3)] = centmass3(coords.img==i) ; 
end
elocs(:,1) = 256-elocs(:,1) ; 
temp = elocs(:,2) ; elocs(:,2) = elocs(:,3) ; elocs(:,3) = temp ; 
elocs(:,3) = 256 - elocs(:,3) ; 
sub = load('C:\brainstorm_db\Protocol01\anat\alex\subjectimage_T1.mat') ; 
elocs(66:68,:) = fids ; 
conv = cs_convert(sub,'voxel','scs',elocs(:,1:3)') ; 
%theta = -pi/4 ; 
%conv = conv*[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
%conv = elocs ; 
elecorder = load('C:\shared\freesurfer_segs\elecorder.mat') ; elecorder = elecorder.elecorder ;
%elecorder{66} = 'NAS' ; elecorder{67} = 'LPA' ; elecorder{68} = 'RPA' ; 
%conv = elocs ; conv(66:68,:) = fids ; 
cd ..
fid = fopen('locs.txt','wt');
for i=1:65 ; fprintf(fid, '%s %f%f%f\n', elecorder{i}, conv(i,1),conv(i,2),conv(i,3)) ; end
fclose(fid);









