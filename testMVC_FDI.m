%--------------------------------------------------------------------------
% testMVC_FDI.m
% Author: Akira Nagamori
% Last update: 2/21/18
% Code descriptions
% Ojbective: Test maximal force level achieved by the twitch-based model
%       * tendon not included due to high sampling rate requirment
%   
%--------------------------------------------------------------------------

close all
clear all
clc

codeFolder = '/Users/akira/Documents/Github/Afferented-Muscle-Development';
dataFolder = '/Volumes/DATA2/TwoMuscleSystemData/CombinedModel/Model Development';
modelParameter.optimalLength = 3;
modelParameter.mass = 0.0001;

Fs = 1000;
t = 0:1/Fs:8;
amp_temp = 0.1:0.1:1;


amp = 1;

input = [zeros(1,1*Fs) amp*[0:1/Fs:1] amp*ones(1,5*Fs),zeros(1,1*Fs)];
    
output = muscleModel_Combined_FDI(t,Fs,input,modelParameter,1);

maxForce = mean(output.Force_total(6*Fs:7*Fs))
CoV = std(output.Force_total(4*Fs:5*Fs))/mean(output.Force_total(4*Fs:5*Fs))


figure(1)
plot(t,output.Force_total)


