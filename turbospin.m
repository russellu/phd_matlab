% BRIEF EXPLANATION: (in my own words)
% turbo spin echo or fast spin echo works by manipulating the k-space
% through applying multiple phase encoding gradients within a single TR, 
% allowing us to fill different y lines of K-space within a single TR. we
% can do this by reversing the phase encoding gradient before each
% consecutive echo, and applying a new one. thus, if we use 8 echos within
% a single TR, we can fill 8 different y lines of k space within a single
% TR.

% EFFECTIVE TE: effective TE is the time in ms between the largest echo and
% the start of the pulse sequence. in turbo spin echo, effective TE is
% controlled by how we arrange our phase encoding gradients, we place our
% lowest strength phase encoding gradient directly before the echo we get
% at our effective TE, because that produces the largest echo (due to less 
% dephasing). 

% advantage compared to regular SE sequence: faster filling of k space
% disadvantage: TE is not as precise, so the contrast may not be as good or
% precise as if you used one echo per TR

% k space diagram: fill top to bottom, effective te around 1/2 of way in

% initialize pulse/gradient/DET vectors
pulses = zeros(1,200) ; % RF pulse
xgrad = zeros(1,200) ; % frequency encoding
ygrad = zeros(1,200) ; % phase encoding
zgrad = zeros(1,200) ; % slice selection
DET = zeros(1,200) ; % detector

% 8 refocusing pulse plus first 90 degree pulse
pulses(11) = 1 ; pulses(21) = 2 ; pulses(41) = 2 ; pulses(61) = 2 ; pulses(81) = 2 ; pulses(101) = 2 ; pulses(121) = 2 ; pulses(141) = 2 ; pulses(161) = 2 ; 
% z gradients (slice selection)
zgrad(10:13) = 1 ; zgrad(13:14)=-1 ; 
zgrad(20:23) = 1 ; zgrad(40:43) = 1 ; zgrad(60:63) = 1 ;  zgrad(80:83) = 1 ; zgrad(100:103) = 1 ; zgrad(120:123) = 1 ; zgrad(140:143) = 1 ; zgrad(160:163) = 1 ; 
% y gradients (phase encoding)
ygrad(24:27) = 6 ; ygrad(44:47) = 5 ; ygrad(64:67) = 2 ; ygrad(84:87) = 1 ; ygrad(104:107) = 3 ; ygrad(124:127) = 4 ; ygrad(144:147) = 6 ; ygrad(164:167) = 8 ; 
ygrad(35:38) = -6 ; ygrad(55:58) = -5 ; ygrad(75:78) = -2 ; ygrad(95:98) = -1 ; ygrad(115:118) = -3 ; ygrad(135:138) = -4 ; ygrad(155:158) = -6 ; ygrad(175:178) = -8 ; 
% x gradients (frequency encoding)
xgrad(28:31) = 1 ; xgrad(48:51) = 1 ; xgrad(68:71) = 1 ; xgrad(88:91) = 1 ; xgrad(108:111) = 1 ; xgrad(128:131) = 1 ; xgrad(148:151) = 1 ; xgrad(168:171) = 1 ; 
% detector
DET(30) = 1 ; DET(50) = 1 ; DET(70) = 1 ; DET(90) = 1 ; DET(110) = 1 ; DET(130) = 1 ; DET(150) = 1 ; DET(170) = 1 ; 

zgrad = zgrad - 4 ; 
ygrad = ((ygrad./6)) - 8;
xgrad = xgrad - 12 ; 
DET = DET - 16 ;
figure,
hold on
plot(pulses,'LineWidth',2) ; 
plot(zgrad,'LineWidth',2) ; 
plot(ygrad,'LineWidth',2) ; 
plot(xgrad,'LineWidth',2) ; 
plot(DET,'LineWidth',2) ; 
vline(find(pulses==2),'r') ; 
vline(find(DET==1-16),'g') ; 
suptitle('green lines = detector, red lines = 180deg RF pulse, row1=RF,row2=slice select,row3=phase encode,row4=frequency encode,row5=detector') ;

% k-space diagram
[x,y] = meshgrid(-127:127,-127:127) ; 


