cd C:\Users\butr2901\Downloads ;

[v,fs] = mp3read('Celtic Music - For the King.mp3') ; 

v2 = eegfiltfft(v',fs,2000,20000) ; 
p = audioplayer(v2,fs) ; 

f = abs(fft(v(:,1))) ; 

