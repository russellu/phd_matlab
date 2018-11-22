clear all
 
cd('c:/shared/allres/russell')
EEG = pop_loadset('ica_notch85.set');
EEG = eeg_checkset(EEG);
 
%this removes only the bad components
EEG = pop_subcomp( EEG, [1   2   3   4   5   6   7   8   9  10  11  12  13  15  16  17  19  20  21  22  23  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64], 0);
EEG = eeg_checkset(EEG);
clear t1 t2 t3 m1
t1=round(EEG.srate);
t2=round(EEG.srate+EEG.srate/2)  ; 
t3=round(2.5.*EEG.srate)   ;
m1=0.15
 
 
EEG1 = pop_epoch( EEG, {'S 12'}, [-1 2.5]);
EEG1 = eeg_checkset(EEG1);
clear mod
for i=1:EEG1.trials
clear r t
r=double(EEG1.data(:,1:t1+1,i));
t=double(EEG1.data(:,t2:t3,i));
r=abs(eegfiltfft(r,EEG.srate,60,80));
t=abs(eegfiltfft(t,EEG.srate,60,80));
rest(i,:)=mean(r,2);
task(i,:)=mean(t,2);
end
subplot(1,3,1)
topoplot(mean(task)-mean(rest),EEG.chanlocs,'maplimits',[0 m1])
title('contrast 5%')
mod1=mean(task)-mean(rest);
EEG1 = pop_epoch( EEG, {'S 11'}, [-1 2.5]);
EEG1 = eeg_checkset(EEG1);
clear mod
for i=1:EEG1.trials
clear r t
r=double(EEG1.data(:,1:t1+1,i));
t=double(EEG1.data(:,t2:t3,i));
r=abs(eegfiltfft(r,EEG.srate,60,80));
t=abs(eegfiltfft(t,EEG.srate,60,80));
rest(i,:)=mean(r,2);
task(i,:)=mean(t,2);
end
subplot(1,3,2)
topoplot(mean(task)-mean(rest),EEG.chanlocs,'maplimits',[0 m1])
title('Unp')
mod2=mean(task)-mean(rest);
 
EEG1 = pop_epoch( EEG, {'S 16'}, [-1 2.5]);
EEG1 = eeg_checkset(EEG1);
clear mod
for i=1:EEG1.trials
clear r t
r=double(EEG1.data(:,1:t1+1,i));
t=double(EEG1.data(:,t2:t3,i));
r=abs(eegfiltfft(r,EEG.srate,60,80));
t=abs(eegfiltfft(t,EEG.srate,60,80));
rest(i,:)=mean(r,2);
task(i,:)=mean(t,2);
end
subplot(1,3,3)
topoplot(mean(task)-mean(rest),EEG.chanlocs,'maplimits',[0 m1])
title('rand 60%')
mod3=mean(task)-mean(rest);
 
clear ch
ch=find(mod2>prctile(mod2,75))
[H P]=ttest(mod2(ch),mod1(ch),0.05,'right')
[H P]=ttest(mod2(ch),mod3(ch),0.05,'right')
 