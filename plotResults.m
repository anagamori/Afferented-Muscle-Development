close all
clear all
clc

codeFolder = '/Users/akira/Documents/Github/Afferented-Muscle-Development';
dataFolder = '/Volumes/DATA2/TwoMuscleSystemData/CombinedModel/Model Development';

Fs = 1000;
t = 0:1/Fs:8;
amp_temp = 0.1:0.1:1;

for i = 1:length(amp_temp)
    trialN = i;
    cd (dataFolder)
    load(['output_Song_' num2str(trialN)])
    load(['output_Tsianos_' num2str(trialN)])
    load(['output_MU_' num2str(trialN)])
    cd (codeFolder)
    
    maxForce1(i) =  mean(output_1.Force_total(6*Fs:7*Fs));
    maxForce2(i) = mean(output_2.Force_total(6*Fs:7*Fs));
    maxForce3(i) = mean(output_3.Force_total(6*Fs:7*Fs));
    CoV(i) = std(output_3.Force_total(4*Fs:5*Fs))/mean(output_3.Force_total(4*Fs:5*Fs));

    
end

figure(1)
plot(amp_temp,maxForce1./maxForce1(end))
hold on
plot(amp_temp,maxForce2./maxForce2(end))
hold on
plot(amp_temp,maxForce3./maxForce3(end))
legend('Song','Tsianos','Motor Unit')
xlabel('Activation','FontSize',14)
ylabel('%MVC','FontSize',14)