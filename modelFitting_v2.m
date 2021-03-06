%--------------------------------------------------------------------------
% Author: Akira Nagamori
% Last update: 12/18/2017
% Model descriptions:
%   Generate model parameters given the distributino of peak firing rates
%   across units, muscle architecture
%--------------------------------------------------------------------------


close all
clear all
clc

%--------------------------------------------------------------------------
% muscle architecture parameters

density = 1.06; % muscle density [g/cm^3]
L0 = 6.8; % optimal muscle length [cm]
mass = 0.15; % muscle mass [kg]
PCSA = (mass*1000)/(density*L0); % PCSA of muscle
sigma = 22.5; % specific tension
F0 = PCSA*sigma;
N_fibers = 272850; % number of muscle fibers
N_MU = 450; % number of motor units
i_MU = 1:N_MU; % index for motor units

%--------------------------------------------------------------------------
% Find peak tetanus amplitude
RP_MU = 100; %range of peak tension across motor untis in unit of fold
b_MU = log(RP_MU)/N_MU; %coefficient to establish a range of twich force values
P_MU = exp(b_MU*i_MU); %force generated by a motor unit as a function of its recruitment threshold

fPCSA = PCSA/N_fibers; % fractional PCSA occupied by a single fiber
ni = P_MU/sum(P_MU)*N_fibers; % innervation number of each unit assuming that peak tension generated by each unit is proportional to innervation number
Pi = sigma*fPCSA*ni; % peak tetanus amplitude of each unit

% load('cor_factor')
% Pi = Pi.*cor_factor';
%-------------------------

F_pcsa_slow = 0.5; % fractional PSCA of slow-twitch motor units (0-1)
[error, index_slow] = min(abs(cumsum(Pi) - F0*F_pcsa_slow));

% Find recruitment threshold for individual units using exponential fit
Ur = 0.8; % recruitment threshold for the lastly recruited motor unit
Ur_1 = 0.01; % reruitment threshold for the first unit
f_RT = fit([1 N_MU]',[Ur_1 Ur]','exp1');
coeffs_f_RT = coeffvalues(f_RT);
U_th = coeffs_f_RT(1)*exp(coeffs_f_RT(2)*i_MU); % the resulting recruitment threshold for individual units

%--------------------------------------------------------------------------
% Find peak firing rates of individual units
%  In this formulation, the peak firing rate of a unit increases as
%  recruitment threshold increases (Loeb's formulation)
%  The minimum firing rate for individual units is constant as opposed to
%  Loeb's formulation

MFR_MU = 8; %muscle_parameter.MFR; %minimum firing rate constant for all motoneurons
g_e_MU = 1; %missing parameter from the paper

PFR1_MU = 20; %the peak firing rate of the first recruited motoneuron in unit of impulse/sec
PFRD_MU = 40; %the desired difference in peak firing rates between the first and last units in unit of impulse/sec
RTEn_MU = U_th(end); %recruitment threshold of the last motor unit
PFR_MU = PFR1_MU + PFRD_MU * (U_th./RTEn_MU); %peak firing rate
FR_half = PFR_MU./2; % firing rate at which half of maximum tension is achieved

%--------------------------------------------------------------------------
% Find contraction time of individual units
%   Using Loeb's formulation that contraction time of individual units is
%   proportional to 1/FR_half (Brown & Loeb 2000 IV; Botterman et al., 1996)
CT_n = 11;
FR_half_n = FR_half(end);
CT = 3*(CT_n*FR_half_n)./FR_half;
CT = CT - (CT(end)-CT_n);
CT = CT/1000;
RT = CT.*1.3719;
cv_MU = 0; %ISI variability as per coefficient of variation (=mean/SD)

%--------------------------------------------------------------------------
Lce = 1;
Y = 1;
S = 0.96;
Af = 1;
Fs = 1000;
t_twitch = 0:1/Fs:1;
t_temp = 0:1/Fs:8;
t = 0:1/Fs:10;
count = 0;
testedUnits = [1 50 100 150 200 250 300 350 400 450];
U = [zeros(1,Fs) ones(1,8*Fs) zeros(1,1*Fs+1)];

for k = 1:length(testedUnits)
    
    unitN = testedUnits(k);
    FR = 8:0.5:PFR_MU(unitN);
    f_env_temp = FR./FR_half(unitN);
    
    %cor_factor*
    %twitch = cor_factor(unitN).*Pi(unitN).*t_twitch./CT(unitN).*exp(1-t_twitch./CT(unitN));
    meanForce = zeros(1,length(FR));
    for j = 1:length(FR)
        f_env = f_env_temp(j).*U;
        f_int = 0;
        f_int_dot = 0;
        f_eff = 0;
        f_eff_dot = 0;
        f_int_vec = zeros(1,length(t));
        f_eff_vec = zeros(1,length(t));
        Tf = zeros(1,length(t));
        spikeTrain_temp = spikeTrainGenerator(t_temp,Fs,FR(j));
        spikeTrain = [zeros(1,1*Fs) spikeTrain_temp zeros(1,1*Fs)];
        force = zeros(1,length(t));
        for i = 1:length(t)
            spike = zeros(1,length(t));
            [f_int,f_int_dot,Tf(i)] = f_function(f_int,f_env(i),f_env(i),f_eff_dot,CT(unitN),RT(unitN),Af,Lce,Fs);
            [f_eff,f_eff_dot,~] = f_function(f_eff,f_int,f_env(i),f_eff_dot,CT(unitN),RT(unitN),Af,Lce,Fs);
            if unitN < index_slow
                Af = Af_slow_function(f_eff,Lce,Y);
            else
                Af_fast_function(f_eff,Lce,S);
            end
            f_int_vec(i) = f_int;
            f_eff_vec(i) = f_eff;
            if spikeTrain(i) == 1
                spike(i) = 1;
                twitch = Pi(unitN).*t_twitch./Tf(i).*exp(1-t_twitch./Tf(i));
                force_temp = conv(spike,twitch);
                force = force + force_temp(1:length(t));
            end
        end
        
        
        %         force = conv(spikeTrain,twitch);
        %
        meanForce(j) = mean(force(5*Fs:7*Fs));
        %
        figure(1)
        plot(t,force)
        hold on
    end
    
    figure(2)
    plot(FR,meanForce)
    hold on
    
    maxForce = meanForce(end);
    ratio(k) = Pi(unitN)/maxForce;
    
end

figure(3)
plot(ratio)

x = testedUnits';
y = ratio';
f = fit(x,y,'smoothingspline');
figure(4)
plot(f,x,y)

cor_factor = feval(f,i_MU);
save('cor_factor','cor_factor')

function [f_out,f_out_dot,Tf] = f_function(f_out,f_in,f_env,f_eff_dot,CT,RT,Af,Lce,Fs)
T_f1 = CT*3/4;
T_f2 = CT*1/4;
T_f3 = RT*3/4;
T_f4 = RT*1/4;

if f_eff_dot >= 0
    Tf = T_f1*Lce^2+T_f2*f_env;
else
    Tf = (T_f3 + T_f4*Af)/Lce;
end
f_out_dot = (f_in - f_out)/Tf;
f_out = f_out_dot*1/Fs + f_out;

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