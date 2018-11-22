subs = {'alex','dina','genevieve','jeremie','karl','russell','sukhman','valerie'} ; 
for sub=1:length(subs) 

cd(['c:/shared/badger_mri/',subs{sub},'/nii/allmelodic']) ; ls 
mix = load('melodic_mix') ; 
names = {
'reg_topup_mc_retino_allstims_01.txt' 
'reg_topup_mc_retino_allstims_02.txt'
'reg_topup_mc_retino_gamma_01.txt'
'reg_topup_mc_retino_gamma_02.txt'
'reg_topup_mc_retino_movie.txt'
'reg_topup_mc_retino_rest.txt'
} ; 

a = 1:735:735*5 ;
ts1 = a(1):a(1)+734 ; 
ts2 = a(2):a(2)+734 ; 
ts3 = a(3):a(3)+734 ; 
ts4 = a(4):a(4)+734 ; 
ts6 = size(mix,1)-449:size(mix,1) ; 

newmix1 = mix(ts1,:) ; 
newmix2 = mix(ts2,:) ; 
newmix3 = mix(ts3,:) ; 
newmix4 = mix(ts4,:) ; 
newmix6 = mix(ts6,:) ; 

dlmwrite(names{1},newmix1,' ') ; 
dlmwrite(names{2},newmix2,' ') ; 
dlmwrite(names{3},newmix3,' ') ; 
dlmwrite(names{4},newmix4,' ') ; 
dlmwrite(names{6},newmix6,' ') ; 

if sub ~= 5
    ts5 = a(5):a(5)+734 ; 
    newmix5 = mix(ts5,:) ; 
    dlmwrite(names{5},newmix5,' ') ; 
end

end
