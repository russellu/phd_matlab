revdir = 'E:\imeka_data\imeka_florence\revEEG-fMRI_3mm_MB6_Series1201\revEEG-fMRI_3mm_MB6_Series1201' ; % directory with rev image
revname = 'Mag (0001)' ; % name of rev DICOM
forwarddir = 'E:\imeka_data\imeka_florence\EEG-fMRI_3mm_MB6_Series1101\EEG-fMRI_3mm_MB6_Series1101' ; % directory with fwd image
forwardname = 'Mag (0001)' ; % name of fwd DICOM
name = 'vec.txt' ; % txt file name (--datain arg. for topup)
savedir = 'E:\imeka_data\imeka_florence' ; 

% for the forward encoding
cd(forwarddir) ; 
disp(['loading: ',forwardname]) ; 
dinfo = dicominfo(forwardname) ; 
trainl = double(dinfo.EchoTrainLength) ; 
waterfatshift = dinfo.Private_2001_1022 ; 
wfs = 434.215 ; 
echospacing = ((1000*waterfatshift) / (wfs*(trainl+1))) ;
dwelltime = echospacing * trainl ; 
topup = echospacing * (trainl-1) ; 

% for the reverse encoding
cd(revdir) ; 
disp(['loading: ',revname]) ; 
dinfo2 = dicominfo(revname) ; 
trainl2 = double(dinfo2.EchoTrainLength) ; 
waterfatshift2 = dinfo2.Private_2001_1022 ; 
wfs2 = 434.215 ; 
echospacing2 = ((1000*waterfatshift2)/(wfs2*(trainl2+1))) ;
dwelltime2 = echospacing2 * trainl2 ; 
topup2 = echospacing2 * (trainl2-1) ; 

vec = [0,1,0,topup/1000 ; 0,-1,0,topup2/1000] ; 
cd(savedir) ; 
dlmwrite(name,vec,' ')

