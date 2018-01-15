%--------------------------------------------------------------------------
% muscleModel_Combined.m
% Author: Akira Nagamori
% Last update: 12/26/207
%--------------------------------------------------------------------------

function output = muscleModel_Combined_v3(t,Fs,input,modelParameter,simulationParameter,recruitmentType)
%--------------------------------------------------------------------------
% Define model parameters
%--------------------------------------------------------------------------
% F0
density = 1.06;
L0 = modelParameter.L0; % optimal muscle length [cm]
mass = modelParameter.mass; % muscle mass [kg]
PCSA = (mass*1000)/(density*L0); % PCSA of muscle
sigma = 31.8; % specific tension
F0 = PCSA * sigma; % maximal force

%--------------------------------------------------------------------------
% Motor unit architecture
N_fibers = 272850; % number of muscle fibers
N_MU = 450; % number of motor units
i_MU = 1:N_MU; % index for motor units

%--------------------------------------------------------------------------
% Peak tetanic force
RP_MU = 30; %range of peak tension across motor untis in unit of fold
b_MU = log(RP_MU)/N_MU; %coefficient to establish a range of twich force values
P_MU = exp(b_MU*i_MU); %force generated by a motor unit as a function of its recruitment threshold
fPCSA = PCSA/N_fibers; % fractional PCSA occupied by a single fiber
ni = P_MU/sum(P_MU)*N_fibers; % innervation number of each unit assuming that peak tension generated by each unit is proportional to innervation number
Pi = sigma*fPCSA*ni; % peak tetanus amplitude of each unit
Pi_half = Pi./2; % half peak tetanus amplitude
cor_factor = modelParameter.cor_factor;
Pi = cor_factor; % adjusted peak twitch amplitude


%--------------------------------------------------------------------------
% Recruitment
F_pcsa_slow = 0.5; % fractional PSCA of slow-twitch motor units (0-1)
[~, index_slow] = min(abs(cumsum(Pi) - F0*F_pcsa_slow));

% Find recruitment threshold for individual units using exponential fit
Ur = 0.8; % recruitment threshold for the lastly recruited motor unit
Ur_1 = 0.01; % reruitment threshold for the first unit
f_RT = fit([1 N_MU]',[Ur_1 Ur]','exp1');
coeffs_f_RT = coeffvalues(f_RT);
U_th = coeffs_f_RT(1)*exp(coeffs_f_RT(2)*i_MU); % the resulting recruitment threshold for individual units

%--------------------------------------------------------------------------
% Minimum and peak firing rate
MFR_MU = 8; %muscle_parameter.MFR; %minimum firing rate constant for all motoneurons
PFR1_MU = 20; %the peak firing rate of the first recruited motoneuron in unit of impulse/sec
PFRD_MU = 40; %the desired difference in peak firing rates between the first and last units in unit of impulse/sec
RTEn_MU = U_th(end); %recruitment threshold of the last motor unit
PFR_MU = PFR1_MU + PFRD_MU * (U_th./RTEn_MU); %peak firing rate
FR_half = PFR_MU./2; % firing rate at which half of maximum tension is achieved

g_e = (PFR_MU(end)-MFR_MU)/(1-Ur);

%--------------------------------------------------------------------------
% Contraction time and half relaxation time
CT_n = 20;
FR_half_n = FR_half(end);
CT = 1.5*(CT_n*FR_half_n)./FR_half;
CT = CT - (CT(end)-CT_n);
CT = CT/1000;
RT = CT;
t_twitch = 0:1/Fs:2;
%--------------------------------------------------------------------------
cv_MU = 0; %ISI variability as per coefficient of variation (=mean/SD)
U = input;
U_eff = 0;
T_U = 0.03;

%--------------------------------------------------------------------------
c_y = 0.35;
V_y = 0.1;
T_y = 0.2;
T_s = 0.043;
a_f = 0.56;
n_f0 = 2.1;
n_f1_slow = 5;
n_f1_fast = 3.3;
a_f_FF = 0.52;
n_f0_FF = 1.97;
n_f1_FF_slow = 5.1;
n_f1_FF_fast = 3.28;
beta_slow = 2.3;
omega_slow = 1.12;
rho_slow = 1.62;
beta_fast = 1.55;
omega_fast = 0.75;
rho_fast = 2.12;
Vmax_slow = -7.88;
cv0_slow= 5.88;
cv1_slow = 0;
Vmax_fast = -9.15;
cv0_fast = -5.7;
cv1_fast = 9.18;
av0_slow = -4.7;
av1_slow = 8.41;
av2_slow = -5.34;
bv_slow = 0.35;
av0_fast = -1.53;
av1_fast = 0;
av2_fast = 0;
bv_fast = 0.69;


%--------------------------------------------------------------------------
% initialize parameters
Lce = simulationParameter.Lce;
Vce = 0;

spike_time = zeros(N_MU,1);
spike_train = cell(1,N_MU);
force = cell(1,N_MU);
for j = 1:N_MU
    spike_train{j} = zeros(1,length(t));
    force{j} = zeros(1,length(t));
end
Y = zeros(1,length(t));
S = zeros(1,length(t));
%--------------------------------------------------------------------------
% Simulation
for i = 1:length(t)
    U_eff_dot = (U(i) - U_eff)/T_U;
    U_eff = U_eff_dot*1/Fs + U_eff;
    if recruitmentType == 1 % linear increase in firing rate up to Ur
        FR_MU = (PFR_MU-MFR_MU)./(1-U_th).*(U_eff-U_th) + MFR_MU;
    elseif recruitmentType == 2 % equal gain across units and saturation
        FR_MU = g_e*(U(i)-U_th)+ MFR_MU;
    end
    FR_MU(FR_MU<8) = 0;
    f_env = FR_MU./FR_half;
    
    parfor n = 1:N_MU
        force_temp = force{n};
        force_half = force_temp(i)/Pi_half(n);
        if n <= index_slow
            Y_dot = (1-c_y*(1-exp(-abs(Vce)/V_y))-Y(n))/T_y;
            Y(n) = Y_dot*1/Fs + Y(n);
            n_f_slow = n_f0 + n_f1_slow*(1/Lce-1);
            Af = 1 - exp(-(Y(n)*force_half/(a_f*n_f_slow))^n_f_slow);
            n_f_FF_slow = n_f0_FF + n_f1_FF_slow*(1/Lce-1);
            FF = 1 - exp(-(Y(n)*f_env(n)/(a_f_FF*n_f_FF_slow))^n_f_FF_slow);
            FF = FF/f_env(n);
            FL = exp(-abs((Lce^beta_slow - 1)/omega_slow)^rho_slow);
            if Vce < 0
                FV = (Vmax_slow - Vce)/(Vmax_slow + (cv0_slow + cv1_slow*Lce)*Vce);
            elseif Vce >= 0
                FV = (bv_slow - (av0_slow + av1_slow*Lce + av2_slow*Lce^2)*Vce)/(bv_slow+Vce);
            end
            
        else
            if force_half < 0.1
                a_s = 1.76;
            else
                a_s = 0.96;
            end
            S_dot = (a_s-S(n))/T_s;
            S(n) = S_dot*1/Fs + S(n);
            n_f_fast = n_f0 + n_f1_fast*(1/Lce-1);
            Af = 1 - exp(-(S(n)*force_half/(a_f*n_f_fast))^n_f_fast);
            n_f_FF_fast = n_f0_FF + n_f1_FF_fast*(1/Lce-1);
            FF = 1 - exp(-(S(n)*f_env(n)/(a_f_FF*n_f_FF_fast))^n_f_FF_fast);
            FF = FF/f_env(n);
            FL = exp(-abs((Lce^beta_fast - 1)/omega_fast)^rho_fast);
            if Vce < 0
                FV = (Vmax_fast - Vce)/(Vmax_fast + (cv0_fast + cv1_fast*Lce)*Vce);
            elseif Vce >= 0
                FV = (bv_fast - (av0_fast + av1_fast*Lce + av2_fast*Lce^2)*Vce)/(bv_fast+Vce);
            end
        end
        if FR_MU(n) > MFR_MU
            
            spike_train_temp = zeros(1,length(t));
            if ~any(spike_train{n}) % when the motor unit fires at the first time
                spike_train_temp(i) = 1;
                spike_train{n} = spike_train_temp+spike_train{n}; % add a spike to the vector
                mu = 1/FR_MU(n);
                Z = randn(1);
                Z(Z>3.9) = 3.9;
                Z(Z<-3.9) = -3.9;
                spike_time_temp = (mu + mu*cv_MU*Z)*Fs;
                spike_time(n) = round(spike_time_temp) + i;
                T1 = CT(n)*Lce^2+(CT(n)*1/2)*Af;
                T2_temp = (RT(n) + (RT(n)*1/2)*Af)/Lce;
                T2 = T2_temp/1.68;
                f_1 = t_twitch./T1.*exp(1-t_twitch./T1);
                f_2 = t_twitch./T2.*exp(1-t_twitch./T2);
                twitch_temp = [f_1(1:round(T1*Fs+1)) f_2(round(T2*Fs+1):end)];
                twitch_temp = twitch_temp(1:1*Fs);
                twitch =  Pi(n).*twitch_temp*FF*FL*FV;
                force_temp = conv(spike_train_temp,twitch);
                force{n} = force{n} + force_temp(1:length(t));
            else
                if spike_time(n) == i % when the motor unit fires
                    spike_train_temp(i) = 1;
                    spike_train{n} = spike_train_temp+spike_train{n};
                    % update mean firing rate of the motor unit given the
                    % current value of input
                    mu = 1/FR_MU(n); % interspike interval
                    Z = randn(1);
                    Z(Z>3.9) = 3.9;
                    Z(Z<-3.9) = -3.9;
                    spike_time_temp = (mu + mu*cv_MU*Z)*Fs; % interspike interval
                    spike_time(n) = round(spike_time_temp) + i; % spike time of next spike
                    
                    T1 = CT(n)*Lce^2+(CT(n)*1/2)*Af;
                    T2_temp = (RT(n) + (RT(n)*1/2)*Af)/Lce;
                    T2 = T2_temp/1.68;
                    f_1 = t_twitch./T1.*exp(1-t_twitch./T1);
                    f_2 = t_twitch./T2.*exp(1-t_twitch./T2);
                    twitch_temp = [f_1(1:round(T1*Fs+1)) f_2(round(T2*Fs+1):end)];
                    twitch_temp = twitch_temp(1:1*Fs);
                    twitch =  Pi(n).*twitch_temp*FF*FL*FV;
                    force_temp = conv(spike_train_temp,twitch);
                    force{n} = force{n} + force_temp(1:length(t));
                elseif i > spike_time(n) + round(1/FR_MU(n)*Fs)
                    spike_train_temp(i) = 1;
                    spike_train{n} = spike_train_temp+spike_train{n};
                    spike_time(n) = i;
                    mu = 1/FR_MU(n); % interspike interval
                    Z = randn(1);
                    Z(Z>3.9) = 3.9;
                    Z(Z<-3.9) = -3.9;
                    spike_time_temp = (mu + mu*cv_MU*Z)*Fs; % interspike interval
                    spike_time(n) = round(spike_time_temp) + i; % spike time of next spike
                    
                    T1 = CT(n)*Lce^2+(CT(n)*1/2)*Af;
                    T2_temp = (RT(n) + (RT(n)*1/2)*Af)/Lce;
                    T2 = T2_temp/1.68;
                    f_1 = t_twitch./T1.*exp(1-t_twitch./T1);
                    f_2 = t_twitch./T2.*exp(1-t_twitch./T2);
                    twitch_temp = [f_1(1:round(T1*Fs+1)) f_2(round(T2*Fs+1):end)];
                    twitch_temp = twitch_temp(1:1*Fs);
                    twitch =  Pi(n).*twitch_temp*FF*FL*FV;
                    force_temp = conv(spike_train_temp,twitch);
                    force{n} = force{n} + force_temp(1:length(t));
                end
            end
        end
    end
end




%Force = sum(force);
%output.Force_total = Force;
output.force = force;
output.spike_train = spike_train;
% output.f_env = Outputfenv;
% output.force_half = Outputforce_half;
% output.Af = OutputAf;
% output.FF = OutputFF;
% output.S = S_af;
% output.Y = Y_af;

end

