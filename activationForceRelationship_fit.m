%--------------------------------------------------------------------------
% activationForceRelationship_fit.m
% Author: Akira Nagamori
% Last update: 12/25/17
% Code descriptions
% Ojbective: fit the relationship between activation (in Hz) and force to
% the Af relationship (Af function) formulated by Loeb's model
%   - Define the minimum and peak firing rate of a given motor unit
%   - Compute f_0.5 (firing rate at which half of maximum force is achieved)
%   - Sweep firing rate from minimum (1 Hz) up to 3*f_0.5 and obtain mean
%   force level from summation of twitches 
%   - Find parameters (a,b) of a scaling function that best-fits the Af
%   relationship formulated in the Loeb's model (Song et al. 2008)
%--------------------------------------------------------------------------

close all
clear all
clc
%--------------------------------------------------------------------------
% motor unit parameters
N_MU = 450; % number of motor units
i_MU = 1:N_MU; % index for motor units

Ur = 0.8; % recruitment threshold for the lastly recruited motor unit
Ur_1 = 0.01; % reruitment threshold for the first unit
f_RT = fit([1 N_MU]',[Ur_1 Ur]','exp1');
coeffs_f_RT = coeffvalues(f_RT);
U_th = coeffs_f_RT(1)*exp(coeffs_f_RT(2)*i_MU); % the resulting recruitment threshold for individual units

MFR_MU = 8; %muscle_parameter.MFR; %minimum firing rate constant for all motoneurons
PFR1_MU = 20; %the peak firing rate of the first recruited motoneuron in unit of impulse/sec
PFRD_MU = 40; %the desired difference in peak firing rates between the first and last units in unit of impulse/sec
RTEn_MU = U_th(end); %recruitment threshold of the last motor unit
PFR_MU = PFR1_MU + PFRD_MU * (U_th./RTEn_MU); %peak firing rate
FR_half = PFR_MU./2; % firing rate at which half of maximum tension is achieved

CT_n = 20;
FR_half_n = FR_half(end);
CT = 1.5*(CT_n*FR_half_n)./FR_half;
CT = CT - (CT(end)-CT_n);
CT = CT/1000;
RT = CT;

Fs = 1000;
t = 0:1/Fs:5;
t_temp = 0:1/Fs:3;
%--------------------------------------------------------------------------
% simulation parameters
Lce = 1;
Y = 1;
testingUnit =  450;
FR = 1:1:3*FR_half(testingUnit); %PFR_MU(testingUnit);
%FR = 0:10:100; %PFR_MU(testingUnit);

%[twitch,T1,T2] = twitch_function(1,Af,Lce,CT(testingUnit),RT(testingUnit),Fs);
%% Obtain non-corrected activation-force relationship 
% twitch amplitude = 1
for i = 1:length(FR)
    f_env = FR(i)/FR_half(testingUnit);
    spikeTrain_temp = spikeTrainGenerator(t_temp,Fs,FR(i));
    spikeTrain = [zeros(1,1*Fs) spikeTrain_temp zeros(1,1*Fs)];
    if testingUnit <= 383
        Af = Af_slow_function(f_env,Lce,Y);
    else
        Af = Af_fast_function(f_env,Lce,Y);
    end
    Af_vec(i) = Af;
    [twitch,T1(i),T2(i)] = twitch_function(f_env,Af,Lce,CT(testingUnit),RT(testingUnit),Fs);
    force_temp = conv(spikeTrain,twitch);
    force = force_temp(1:length(t));    
    figure(1)
    plot(twitch)
    hold on       
    figure(2)
    plot(t,force)
    meanForce(i) = mean(force(2.5*Fs:3*Fs));  
    hold on 
end

figure(3)
subplot(2,1,1)
plot(FR/FR_half(testingUnit),T1)
hold on
subplot(2,1,2)
plot(FR/FR_half(testingUnit),T2)
hold on
    
figure(11)
plot(FR/FR_half(testingUnit),meanForce)
xlabel('Frequency (f_{0.5})','FontSize',14)
ylabel('Force (AU)','FontSize',14)

%% Fit a scaling function to Af for muscle length = 1*L0
f_env_test = FR/FR_half(testingUnit);

L = Lce;
n_f0 = 2.1;
if testingUnit <= 383
    n_f1 = 5;
else
    n_f1 = 3.3;
end
n_f = n_f0 + n_f1*(1/L-1);
a_f = 0.56;
Af_temp = 1 - exp(-(Y*f_env_test/(a_f*n_f)).^n_f);

a = 0.1:0.01:10;
b = 0.1:0.01:10;
for j = 1:length(a)
    for k = 1:length(b)
    Af_test = 1 - exp(-(Y*f_env_test/(a(j)*b(k))).^b(k));
    Af_test = Af_test./f_env_test;
    force_new = meanForce.*Af_test;
    force_new = force_new./force_new(end);
    
    error(k,j) = sum(abs(Af_temp - force_new));
    end
end

figure(5)
surf(a,b,error)

[minimum_error,loc] = min(error(:));
[row,col] = ind2sub(size(error),loc);
a(col)
b(row)

Af_best = (1 - exp(-(Y*f_env_test/(a(col)*b(row))).^b(row)))./f_env_test;
force_best = meanForce.*Af_best;
force_best = force_best./force_best(end);

Af_new = Af_temp./f_env_test;
force_af = meanForce.*Af_new;
force_af = force_af./force_af(end);

figure(4)
plot(f_env_test,Af_vec)
xlabel('Frequency (f_{0.5})','FontSize',14)
ylabel('Activation','FontSize',14)
hold on
plot(f_env_test,force_best)
legend('Original Af','Fitted Af')
% hold on 
% plot(f_env_test,Af_new)

figure(6)
plot(f_env_test,Af_new)
xlabel('Frequency (f_{0.5})','FontSize',14)
ylabel('Activation','FontSize',14)
% best fit for  
% a = 0.46
% b = 2.17

%% Fit a scaling function to Af at different muscle lengths
clear error
clear force_new
clear Af_test
clear loc
Lce = 1.1;
for i = 1:length(FR)
    f_env = FR(i)/FR_half(testingUnit);
    spikeTrain_temp = spikeTrainGenerator(t_temp,Fs,FR(i));
    spikeTrain = [zeros(1,1*Fs) spikeTrain_temp zeros(1,1*Fs)];
    if testingUnit <= 383
        Af = Af_slow_function(f_env,Lce,Y);
    else
        Af = Af_fast_function(f_env,Lce,Y);
    end
    Af_vec(i) = Af;
    [twitch,T1(i),T2(i)] = twitch_function(f_env,Af,Lce,CT(testingUnit),RT(testingUnit),Fs);
    force_temp = conv(spikeTrain,twitch);
    force = force_temp(1:length(t));    
    figure(1)
    plot(twitch)
    hold on       
    figure(2)
    plot(t,force)
    meanForce(i) = mean(force(2.5*Fs:3*Fs));  
    hold on 
end

L = Lce;
n_f0 = 2.1;
if testingUnit <= 383
    n_f1 = 5;
else
    n_f1 = 3.3;
end
n_f = n_f0 + n_f1*(1/L-1);
a_f = 0.56;
Af_temp = 1 - exp(-(Y*f_env_test/(a_f*n_f)).^n_f);

a = 0.46;
b = 2.17;
c_temp = 0.1:0.01:10;

for j = 1:length(c_temp)
    c = b + c_temp(j)*(1/L-1);
    Af_test = 1 - exp(-(Y*f_env_test/(a*c)).^c);
    Af_test = Af_test./f_env_test;
    force_new = meanForce.*Af_test;
    force_new = force_new./force_new(end);
    
    error(j) = sum(abs(Af_temp - force_new));
end

figure()
plot(c_temp,error)
[~,loc] = min(error);
c_temp(loc)

c = b + c_temp(loc)*(1/L-1);
Af_best = (1 - exp(-(Y*f_env_test/(a*c)).^c))./f_env_test;
force_best = meanForce.*Af_best;
force_best = force_best./force_best(end);

figure()
plot(f_env_test,Af_vec)
hold on
plot(f_env_test,force_best)

% c = 6.1000 for slow
% c = 4
%% function used
function [twitch,T1,T2_temp] = twitch_function(f_env,Af,Lce,CT,RT,Fs)
T1 = CT*Lce^2+(CT*1/4)*f_env;
T2_temp = (RT + (RT*1/4)*Af)/Lce;
T2 = T2_temp/1.68;
t_twitch = 0:1/Fs:1;
f_1 = t_twitch./T1.*exp(1-t_twitch./T1);
f_2 = t_twitch./T2.*exp(1-t_twitch./T2);

twitch = [f_1(1:round(T1*Fs+1)) f_2(round(T2*Fs+1):end)];
twitch = twitch(1:length(t_twitch));

end

function Af = Af_slow_function(f_eff,L,Y)
a_f = 0.56;
n_f0 = 2.1;
n_f1 = 5;
n_f = n_f0 + n_f1*(1/L-1);
Af = 1 - exp(-(Y*f_eff/(a_f*n_f))^n_f);
end

function Af = Af_fast_function(f_eff,L,S)
a_f = 0.56;
n_f0 = 2.1;
n_f1 = 3.3;
n_f = n_f0 + n_f1*(1/L-1);
Af = 1 - exp(-(S*f_eff/(a_f*n_f))^n_f);
end

function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end