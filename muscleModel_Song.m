%--------------------------------------------------------------------------
% muscleModel_Loeb.m
% Author: Akira Nagamori
% Last update: 12/20/207
%--------------------------------------------------------------------------

function output = muscleModel_Song(t,Fs,input,modelParameter,simulationParameter)

% model parameters
density = 1.06;
L0 = modelParameter.L0; % optimal muscle length [cm]
mass = modelParameter.mass; % muscle mass [kg]
PCSA = (mass*1000)/(density*L0); % PCSA of muscle
sigma = 31.8; % specific tension
F0 = PCSA * sigma; % maximal force

Ur = 0.8; % fractional activation level at which all motor units for a given muscle are recruited (0-1)
F_pcsa_slow = 0.5; % fractional PSCA of slow-twitch motor units (0-1)
F_pcsa_fast = 1-F_pcsa_slow; % fractional PSCA of fast-twitch motor units (0-1)
U1_th = 0.01; % threshold for slow-twitch fiber
U2_th = Ur*F_pcsa_slow;

% activation-frequency relationship (Brown and Loeb 2000)
f_half = 8.5; % frequency at which the motor unit produces half of its maximal isometric force
fmin = 0.5*f_half; % minimum firing frequency of slow-twitch fiber
fmax = 2*f_half; % maximum firing frequency of slow-twitch fiber

f_half_fast = 34;% frequency at which the motor unit produces half of its maximal isometric force
fmin_fast = 0.5*f_half_fast; % minimum firing frequency of fast-twitch fiber
fmax_fast = 2*f_half_fast; % maximum firing frequency of fast-twitch fiber

% simulation parameters
% 'isometric' = Vce = 0
Lce = simulationParameter.Lce;
Vce = 0;
U = input;

%--------------------------------------------------------------------------
% parameter initilization
U_eff = 0;
f_int_slow = 0;
f_eff_slow = 0;
f_eff_slow_dot = 0;
f_int_fast = 0;
f_eff_fast = 0;
f_eff_fast_dot = 0;

Af_slow = 0;
Af_fast = 0;

Y = 0;
S = 0;

%--------------------------------------------------------------------------
% storing variables
f_env_slow_vec = zeros(1,length(t));
f_int_slow_vec = zeros(1,length(t));
f_eff_slow_vec = zeros(1,length(t));
f_eff_slow_dot_vec = zeros(1,length(t));
f_env_fast_vec = zeros(1,length(t));
f_int_fast_vec = zeros(1,length(t));
f_eff_fast_vec = zeros(1,length(t));
f_eff_fast_dot_vec = zeros(1,length(t));

Af_slow_vec = zeros(1,length(t));
Af_fast_vec = zeros(1,length(t));

Y_vec = zeros(1,length(t));
S_vec = zeros(1,length(t));

Force = zeros(1,length(t));
Force_slow = zeros(1,length(t));
Force_fast = zeros(1,length(t));

%--------------------------------------------------------------------------
% simulation
for i = 1:length(t)
        if U(i) >= U_eff
            T_U = 0.03;
        elseif U(i) < U_eff
            T_U = 0.15;
        end
      
        U_eff_dot = (U(i) - U_eff)/T_U;
        U_eff = U_eff_dot*1/Fs + U_eff;
        
        if U_eff < U1_th
            W1 = 0;
        elseif U_eff < U2_th
            W1 = (U_eff - U1_th)/(U_eff - U1_th);
        else
            W1 = (U_eff - U1_th)/((U_eff - U1_th) + (U_eff - U2_th));
        end
        if U_eff < U2_th 
            W2 = 0;
        else
            W2 = (U_eff - U2_th)/((U_eff - U1_th) + (U_eff - U2_th));
        end
        
        % firing frequency input to second-order excitation dynamics of
        % slow-twitch fiber
        if U(i) >=  U1_th
            f_env_slow = (fmax-fmin)/(1-U1_th).*(U(i)-U1_th)+fmin;
            f_env_slow = f_env_slow/f_half;
        else
            f_env_slow = 0;
        end
        
        [f_int_slow,~] = f_slow_function(f_int_slow,f_env_slow,f_env_slow,f_eff_slow_dot,Af_slow,Lce,Fs);
        [f_eff_slow,f_eff_slow_dot] = f_slow_function(f_eff_slow,f_int_slow,f_env_slow,f_eff_slow_dot,Af_slow,Lce,Fs);
        
        if U(i) >= U2_th
            f_env_fast = (fmax_fast-fmin_fast)/(1-U2_th).*(U(i)-U2_th)+fmin_fast;
            f_env_fast = f_env_fast/f_half_fast;
        else
            f_env_fast = 0;
        end
        
        [f_int_fast,~] = f_fast_function(f_int_fast,f_env_fast,f_env_fast,f_eff_fast_dot,Af_fast,Lce,Fs);
        [f_eff_fast,f_eff_fast_dot] = f_fast_function(f_eff_fast,f_int_fast,f_env_fast,f_eff_fast_dot,Af_fast,Lce,Fs);
        
        Y = yield_function(Y,Vce,Fs);
        S = sag_function(S,f_eff_fast,Fs);
        Af_slow = Af_slow_function(f_eff_slow,Lce,Y);
        Af_fast = Af_fast_function(f_eff_fast,Lce,S);
        
        % activation dependent force of contractile elements
        Fce = U_eff*(W1*Af_slow + W2*Af_fast);
        if Fce < 0
            Fce = 0;
        elseif Fce > 1
            Fce = 1;
        end
        
        
        f_env_slow_vec(i) = f_env_slow;
        f_int_slow_vec(i) = f_int_slow;
        f_eff_slow_vec(i) = f_eff_slow;
        f_eff_slow_dot_vec(i) = f_eff_slow_dot;
        f_env_fast_vec(i) = f_env_fast;
        f_int_fast_vec(i) = f_int_fast;
        f_eff_fast_vec(i) = f_eff_fast;
        f_eff_fast_dot_vec(i) = f_eff_fast_dot;
        
        Af_slow_vec(i) = Af_slow;
        Af_fast_vec(i) = Af_fast;
        
        Y_vec(i) = Y;
        S_vec(i) = S;
        
        Force(i) = Fce*F0;
        Force_slow(i) = W1*Af_slow*F0;
        Force_fast(i) = W2*Af_fast*F0;
end
    
output.Force_total = Force;
output.Force_slow = Force_slow;
output.Force_fast = Force_fast;
output.f_env_slow = f_env_slow_vec;
output.f_env_fast = f_env_fast_vec;
output.f_int_slow = f_int_slow_vec;
output.f_int_fast = f_int_fast_vec;
output.f_eff_slow = f_eff_slow_vec;
output.f_eff_fast = f_eff_fast_vec;
output.S = S_vec;
output.Y = Y_vec;
output.Af_slow = Af_slow_vec;
output.Af_fast = Af_fast_vec;

%--------------------------------------------------------------------------
% function used in simulation
function Y = yield_function(Y,V,Fs)
c_y = 0.35;
V_y = 0.1;
T_y = 0.2;

Y_dot = (1-c_y*(1-exp(-abs(V)/V_y))-Y)/T_y;
Y = Y_dot*1/Fs + Y;
end

function S = sag_function(S,f_eff,Fs)
if f_eff < 0.1
    a_s = 1.76;
else
    a_s = 0.96;
end
T_s = 0.043;
S_dot = (a_s-S)/T_s;
S = S_dot*1/Fs + S;
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
Af = 1 - exp(-((S*f_eff)/(a_f*n_f))^n_f);
end

function [f_out,f_out_dot,Tf] = f_slow_function(f_out,f_in,f_env,f_eff_dot,Af,Lce,Fs)
% T_f1 = 0.0484;
% T_f2 = 0.032;
% T_f3 = 0.0664;
% T_f4 = 0.0356;

T_f1 = 0.0343;
T_f2 = 0.0227;
T_f3 = 0.047;
T_f4 = 0.0252;

if f_eff_dot >= 0
    Tf = T_f1*Lce^2+T_f2*f_env;
else
    Tf = (T_f3 + T_f4*Af)/Lce;
end
f_out_dot = (f_in - f_out)/Tf;
f_out = f_out_dot*1/Fs + f_out;

end

function [f_out,f_out_dot,Tf] = f_fast_function(f_out,f_in,f_env,f_eff_dot,Af,Lce,Fs)
% T_f1 = 0.0121;
% T_f2 = 0.008;
% T_f3 = 0.0166;
% T_f4 = 0.0089;

T_f1 = 0.0206;
T_f2 = 0.0136;
T_f3 = 0.0282;
T_f4 = 0.0151;

if f_eff_dot >= 0
    Tf = T_f1*Lce^2+T_f2*f_env;
else
    Tf = (T_f3 + T_f4*Af)/Lce;
end
f_out_dot = (f_in - f_out)/Tf;
f_out = f_out_dot*1/Fs + f_out;

end

end