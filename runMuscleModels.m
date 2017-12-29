close all
clear all
clc

modelParameter.L0 = 6.8;
modelParameter.mass = 0.15;
load('cor_factor')
modelParameter.cor_factor = cor_factor;
simulationParameter.Lce = 1;

Fs = 1000;
t = 0:1/Fs:5;
amp_temp = 0.1:0.1:1;

for i = 1 %:length(amp_temp)
amp = 1; %amp_temp(i);
input = [zeros(1,1*Fs) amp*ones(1,3*Fs),zeros(1,1*Fs+1)];
    
output_1 = muscleModel_Loeb(t,Fs,input,modelParameter,simulationParameter);
output_2 = muscleModel_Combined(t,Fs,input,modelParameter,simulationParameter,1);
%output_3 = muscleModel_Combined(t,Fs,input,modelParameter,simulationParameter,2);

maxForce(i) = max(output_1.Force_total);
%%
figure(1)
plot(t,output_1.Force_total)
hold on 
plot(t,output_2.Force_total)
hold on 
%plot(t,output_3.Force_total)
end

figure(2)
plot(amp_temp,maxForce./maxForce(end))


figure(3)
plot(t,output_2.spike_train(450,:))
hold on
plot(t,output_2.Af(450,:))
hold on
plot(t,output_2.FF(450,:))
