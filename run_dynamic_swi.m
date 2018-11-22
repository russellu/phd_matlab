clear all ; close all; 
cd e:/switest/sorted/tr10s; 
basepath = 'e:/switest/sorted/tr10s/';
Path1 = [basepath,'SWIp_10sec_Series0501_WIP_SWIp_10sec_20180719153814_501_e1.nii'];
Path2 = [basepath,'SWIp_10sec_Series0501_WIP_SWIp_10sec_20180719153814_501_e2.nii'];
Path3 = [basepath,'SWIp_10sec_Series0501_WIP_SWIp_10sec_20180719153814_501_e3.nii'];


allnii = dynamicSWIp(Path1,Path2,Path3); 

swiref = load_untouch_nii('SWIp_10sec_Series0501_WIP_SWIp_10sec_20180719153814_501_e1.nii'); 
swiref.img = allnii; 
save_untouch_nii(swiref,'dyn_swi.nii.gz'); 
