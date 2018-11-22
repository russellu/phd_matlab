clear all ;close all ; 
cd C:\shared\lastute
subs=dir('*') ; subs(1:2) = []  ; 
for sb=1:length(subs)
    cd(['c:/shared/lastute/',subs(sb).name]) ; 
    dilsegs = load_untouch_nii('dilsegcoords.nii.gz') ; 
    disp(numel(unique(dilsegs.img)));
end