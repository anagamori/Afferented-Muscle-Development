%--------------------------------------------------------------------------
% motorUnitModel_Cisi.m (10/13/17)
% The model of motor units described by Cisi & Kohn (J. Computational Neuroscience, 2008)
%   The motor unit model is based on the two-compartment model
%
%--------------------------------------------------------------------------
clear all
close all
clc

Fs = 10000;
t = 0:1/Fs:3;

%--------------------------------------------------------------------------
% Model Parameters
% Table2
modelParameters.r_s = 82.5/2*1e-4; %77.5/2*1e-4 Somatic compartment radius [cm]
modelParameters.l_s = 82.5*1e-4; %77.5*1e-4 Somatic compartment length [cm]
modelParameters.R_m_s = 1.05; %1.15 Somatic memberance specific resistance [kOhm*cm^2]

modelParameters.r_d = 62.5/2*1e-4; %41.5/2*1e-4 Dendritic compartment radius [cm]
modelParameters.l_d = 6789*1e-4; %5500*1e-4 Dendritic compartment length [cm]
modelParameters.R_m_d = 10.65; %14.4 Dendritic memberance specific resistance [kOhm*cm^2]

Ca_Threshold = 2.5; %% [mV]
axonThreshold = 18; % [mA] (page 526)

modelParameters.rheobase = 6.5; %6.5 [nA]

modelParameters.g_Na_bar = 30; %30 % Maximal conductances of sodium current [mS/cm^2]
modelParameters.g_Kf_bar = 4; %4 % Maximal conductances of potassium current  [mS/cm^2]
modelParameters.g_Ks_bar = 34; %16 % [mS/cm^2]

% State variables
time2peak = 0.6; % time to peak for state variables [ms]
modelParameters.alpha_m_bar = 22; %22 peak value of alpha [1/ms]
modelParameters.beta_m_bar = 13; %13 peak value of beta [1/ms]
modelParameters.alpha_h_bar = 0.5; %0.5 peak value of alpha [1/ms]
modelParameters.beta_h_bar = 4; %4 peak value of beta [1/ms]
modelParameters.alpha_n_bar = 1.5; %1.5 peak value of alpha [1/ms]
modelParameters.beta_n_bar = 0.1; % peak value of beta [1/ms]
modelParameters.alpha_q_bar = 1.5; % peak value of alpha [1/ms]
modelParameters.beta_q_bar = 0.05; % 0.025 % peak value of beta [1/ms]

%
RefractoryPeriod_bar = 5; % absolute refractory period [ms]

%--------------------------------------------------------------------------
amp = 0.0011;
duration = 10;

Input = zeros(1,length(t));
Input(1*Fs:2*Fs) =  amp;

output = motorUnitModel_Cisi_function(t,Input,Fs,modelParameters);

[pks,locs] = findpeaks(output.V_s,Fs,'MinPeakHeight',100);
length(pks)

figure(1)
subplot(3,1,1)
plot(t,Input)
%xlim([0 0.01])
subplot(3,1,2)
plot(t,output.V_d)
xlabel('Time (s)')
ylabel('V_d')
%xlim([0 0.01])
subplot(3,1,3)
plot(t,output.V_s)
xlabel('Time (s)')
ylabel('V_s')
%xlim([0 0.01])

figure(2)
subplot(4,1,1)
plot(t,output.m)
ylabel('m')
subplot(4,1,2)
plot(t,output.h)
ylabel('h')
subplot(4,1,3)
plot(t,output.n)
ylabel('n')
subplot(4,1,4)
plot(t,output.q)
ylabel('q')

