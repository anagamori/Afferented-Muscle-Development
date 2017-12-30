close all
clear all
clc

modelParameter.L0 = 6.8;
modelParameter.mass = 0.15;
load('cor_factor')
modelParameter.cor_factor = cor_factor;
simulationParameter.Lce = 1;

Fs = 1000;
t = 0:1/Fs:8;
amp_temp = 0.1:0.1:1;

for i = 1:length(amp_temp)
amp = amp_temp(i);
input = [zeros(1,1*Fs) amp*[0:1/Fs:1] amp*ones(1,5*Fs),zeros(1,1*Fs)];
    
output_1 = muscleModel_Loeb(t,Fs,input,modelParameter,simulationParameter);
output_2 = muscleModel_Combined(t,Fs,input,modelParameter,simulationParameter,1);
%output_3 = muscleModel_Combined(t,Fs,input,modelParameter,simulationParameter,2);

maxForce1(i) =  mean(output_1.Force_total(6*Fs:7*Fs));
% halfForce1 = maxForce1(i)/2;
% [~,loc1] = min(abs(halfForce1-output_1.Force_total(1*Fs:2*Fs)));
% t_0_50_1(i) = t(loc1+1*Fs)-1;

maxForce2(i) = mean(output_2.Force_total(6*Fs:7*Fs));
% halfForce2 = maxForce2(i)/2;
% [~,loc2] = min(abs(halfForce2-output_2.Force_total(1*Fs:2*Fs)));
% t_0_50_2(i) = t(loc2+1*Fs)-1;
CoV(i) = std(output_2.Force_total(4*Fs:5*Fs))/mean(output_2.Force_total(4*Fs:5*Fs));

%%
figure(1)
plot(t,output_1.Force_total)
hold on 
plot(t,output_2.Force_total)
hold on 
%plot(t,output_3.Force_total)
end

figure(2)

plot(amp_temp,maxForce1./maxForce1(end))
hold on
plot(amp_temp,maxForce2./maxForce2(end))

% figure(3)
% plot(amp_temp,t_0_50_1)
% hold on 
% plot(amp_temp,t_0_50_2)

figure(4)
plot(amp_temp,CoV)


figure(3)
plot(t,output_2.spike_train(450,:))
hold on
plot(t,output_2.Af(450,:))
hold on
plot(t,output_2.FF(450,:))

