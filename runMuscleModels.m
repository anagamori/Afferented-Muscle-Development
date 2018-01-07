close all
clear all
clc

codeFolder = '/Users/akira/Documents/Github/Afferented-Muscle-Development';
dataFolder = '/Volumes/DATA2/TwoMuscleSystemData/CombinedModel/Model Development';
modelParameter.L0 = 6.8;
modelParameter.mass = 0.15;
load('cor_factor')
modelParameter.cor_factor = cor_factor;
simulationParameter.Lce = 1;

Fs = 1000;
t = 0:1/Fs:8;
amp_temp = 0.1:0.1:1;

for i = length(amp_temp)
amp = amp_temp(i);
trialN = i+11;
input = [zeros(1,1*Fs) amp*[0:1/Fs:1] amp*ones(1,5*Fs),zeros(1,1*Fs)];
%input = 0.3*sin(2*pi*1.*t-pi/2)+0.4;
    
output_1 = muscleModel_Song(t,Fs,input,modelParameter,simulationParameter);
output_2 = muscleModel_Tsianos(t,Fs,input,modelParameter,simulationParameter);
output_3 = muscleModel_Combined(t,Fs,input,modelParameter,simulationParameter,1);

maxForce1(i) =  mean(output_1.Force_total(6*Fs:7*Fs));
maxForce2(i) = mean(output_2.Force_total(6*Fs:7*Fs));
maxForce3(i) = mean(output_3.Force_total(6*Fs:7*Fs));
CoV(i) = std(output_3.Force_total(4*Fs:5*Fs))/mean(output_3.Force_total(4*Fs:5*Fs));

%%
figure(1)
plot(t,output_1.Force_total)
hold on 
plot(t,output_2.Force_total)
hold on 
plot(t,output_3.Force_total)
hold on 
%plot(t,output_3.Force_total)

% cd (dataFolder)
% save(['output_Song_' num2str(trialN)],'output_1')
% save(['output_Tsianos_' num2str(trialN)],'output_2')
% save(['output_MU_' num2str(trialN)],'output_3')
% cd (codeFolder)
end

figure(2)
plot(amp_temp,maxForce1./maxForce1(end))
hold on
plot(amp_temp,maxForce2./maxForce2(end))
hold on
plot(amp_temp,maxForce3./maxForce3(end))
legend('Song','Tsianos','Motor Unit')

% figure(4)
% plot(amp_temp,CoV)

