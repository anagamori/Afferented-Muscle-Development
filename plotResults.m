close all
clear all
clc

codeFolder = '/Users/akira/Documents/Github/Afferented-Muscle-Development';
dataFolder = '/Volumes/DATA2/TwoMuscleSystemData/CombinedModel/Model Development';

Fs = 1000;
t = 0:1/Fs:8;
amp_temp = 0.1:0.1:1;
trialN_temp = [1 5 8];
for i = 1:3 %:length(amp_temp)
    trialN = trialN_temp(i);
    cd (dataFolder)
    load(['output_Song_' num2str(trialN)])
    load(['output_Tsianos_' num2str(trialN)])
    load(['output_MU_' num2str(trialN)])
    cd (codeFolder)
    
    maxForce1(i) =  mean(output_1.Force_total(6*Fs:7*Fs));
    maxForce2(i) = mean(output_2.Force_total(6*Fs:7*Fs));
    maxForce3(i) = mean(output_3.Force_total(6*Fs:7*Fs));
    CoV(i) = std(output_3.Force_total(4*Fs:5*Fs))/mean(output_3.Force_total(4*Fs:5*Fs));
    Force_Song = output_1.Force_total(5*Fs:7*Fs)-mean(output_1.Force_total(5*Fs:7*Fs));
    [pxx_Song,f] = pwelch(Force_Song,[],[],0:0.1:25,Fs,'power');
    Force_MU = output_3.Force_total(2*Fs:7*Fs)-mean(output_3.Force_total(2*Fs:7*Fs));
    [pxx_MU,f] = pwelch(Force_MU,[],[],0:0.1:25,Fs,'power');
    
    figure(2)
    plot(f,pxx_MU./sum(pxx_MU))
    hold on 
    xlabel('Frequency (Hz)','FontSize',14)
    ylabel('Proportion of Total Power','FontSize',14)

end

legend('10 %MVC','50 %MVC','80 %MVC')
% figure(1)
% plot(amp_temp,maxForce1./maxForce1(end))
% hold on
% plot(amp_temp,maxForce2./maxForce2(end))
% hold on
% plot(amp_temp,maxForce3./maxForce3(end))
% legend('Song','Tsianos','Motor Unit')
% xlabel('Activation','FontSize',14)
% ylabel('%MVC','FontSize',14)

% figure(2)
% plot(f,pxx_Song./sum(pxx_Song))
% hold on 
% plot(f,pxx_MU./sum(pxx_MU))
% legend('Song','Motor Unit')
% xlabel('Frequency (Hz)','FontSize',14)
% ylabel('Power','FontSize',14)
