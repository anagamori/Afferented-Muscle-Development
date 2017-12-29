function output = motorUnitModel_Cisi_function(t,input,Fs,modelParameters)
%--------------------------------------------------------------------------
% user-defined parameters that vary with motor unit size
r_s = modelParameters.r_s; % Somatic compartment radius [cm]
l_s = modelParameters.l_s; % Somatic compartment length [cm]
R_m_s = modelParameters.R_m_s; % Somatic memberance specific resistance [kOhm*cm^2]

r_d = modelParameters.r_d; % Dendritic compartment radius [cm]
l_d = modelParameters.l_d; % Dendritic compartment length [cm]
R_m_d =  modelParameters.R_m_d; % Dendritic memberance specific resistance [kOhm*cm^2]

g_Na_bar_temp = modelParameters.g_Na_bar; % Maximal conductances of sodium current [mS/cm^2]
g_Kf_bar_temp = modelParameters.g_Kf_bar; % Maximal conductances of fast potassium current  [mS/cm^2]
g_Ks_bar_temp = modelParameters.g_Ks_bar; % Maximal conductances of slow potassium current [mS/cm^2]

alpha_m_bar = modelParameters.alpha_m_bar; % peak value of alpha [1/ms]
beta_m_bar = modelParameters.beta_m_bar; % peak value of beta [1/ms]
alpha_h_bar = modelParameters.alpha_h_bar; % peak value of alpha [1/ms]
beta_h_bar = modelParameters.beta_h_bar; % peak value of beta [1/ms]
alpha_n_bar = modelParameters.alpha_n_bar; % peak value of alpha [1/ms]
beta_n_bar = modelParameters.beta_n_bar; % peak value of beta [1/ms]
alpha_q_bar = modelParameters.alpha_q_bar; % peak value of alpha [1/ms]
beta_q_bar = modelParameters.beta_q_bar; % peak value of beta [1/ms]

% constant parameters across motor unit size
E_Na = 120; % Sodium equilibrium potential [mV] (page 6)
E_Kf = -10; % Potassium equilibrium potential [mV] (page 6)
E_Ks = -10; % Potassium equilibrium potential [mV] (page 6)
E_Ca = 140; % Calcium equilibrium potential [mV]
E_l = 0; % Leak equilibrium potential [mV] (page 6)

C_m = 1.0; % Membrane specific capacitance [microF/cm^2] (page 526)
R_i = 70; % Cytoplasm resistivity [Ohm*cm] (page 526)

pulseDuration = 0.6;

% parameters defined by user-defined parameters amd constant parameters
C_s = 2*pi*r_s*l_s*C_m; % Membrane capacitance [microF] (Eq. 8)
g_ls = (2*pi*r_s*l_s)/R_m_s; % Leakage conductance [mS] (Eq. 6)

C_d = 2*pi*r_d*l_d*C_m; % Membrane capacitance [microF] (Eq. 7)
g_ld = (2*pi*r_d*l_d)/R_m_d; % Leakage conductance [mS] (Eq. 5)

% Coupling term
g_c = 2/((R_i*l_d)/(pi*r_d^2)+(R_i*l_s)/(pi*r_s^2)); % Coupling conductance [mS] (Eq. 4)
g_c = g_c/1000;

g_Na_bar = g_Na_bar_temp*2*pi*r_s*l_s; %30e-3*2*pi*r_s*l_s; % Maximal conductances of sodium current [mS/cm^2]
g_Kf_bar = g_Kf_bar_temp*2*pi*r_s*l_s; %4e-3*2*pi*r_s*l_s; % Maximal conductances of potassium current  [mS/cm^2]
g_Ks_bar = g_Ks_bar_temp*2*pi*r_s*l_s; %16e-3*2*pi*r_s*l_s; % [mS/cm^2]

% Voltage threshold 
R_s = 1/g_ls;
R_d = 1/g_ld;
R_c = 1/g_c;
Rm = ((R_s+R_c)*R_d)/(R_c+R_s+R_d); % Theoretical Neuroscience Chapter 6: Model Neurons II: The Rall Model
rheobase = modelParameters.rheobase*1e-9;

V_threshold = Rm*1e3*rheobase; 
V_threshold = V_threshold*1000; %[mV]
%--------------------------------------------------------------------------
% simulation parameter initilization 
V_d = 0; %zeros(1,length(t)); % (page 527)
V_s = 0; %zeros(1,length(t)); % (page 527)
V_d_vec = zeros(1,length(t));
V_s_vec = zeros(1,length(t));
I_ion_vec = zeros(1,length(t));

alpha_m = zeros(1,length(t));
beta_m = ones(1,length(t))*beta_m_bar;
alpha_h = ones(1,length(t))*alpha_h_bar;
beta_h = zeros(1,length(t));
alpha_n = zeros(1,length(t));
beta_n = ones(1,length(t))*beta_n_bar;
alpha_q = zeros(1,length(t));
beta_q = ones(1,length(t))*beta_q_bar;


m = zeros(1,length(t));
h = ones(1,length(t));
n = zeros(1,length(t));
q = zeros(1,length(t));

I_syn_d = zeros(1,length(t)); 
I_syn_s = zeros(1,length(t)); 
I_inj_d = zeros(1,length(t)); 
I_inj_s = input;

for i = 1:length(t)
    if i > 2
        if V_s_vec(i-1) >= V_threshold && V_s_vec(i-2) < V_threshold
            alpha_m(i:i+pulseDuration*Fs/1000) = alpha_m_bar;
            beta_m(i:i+pulseDuration*Fs/1000) = 0;
            alpha_h(i:i+pulseDuration*Fs/1000) = 0;
            beta_h(i:i+pulseDuration*Fs/1000) = beta_h_bar;
            alpha_n(i:i+pulseDuration*Fs/1000) = alpha_n_bar;
            beta_n(i:i+pulseDuration*Fs/1000) = 0;
            alpha_q(i:i+pulseDuration*Fs/1000) = alpha_q_bar;
            beta_q(i:i+pulseDuration*Fs/1000) = 0;
            m_temp = stateVariableDynamics_1(m(i-1),alpha_m_bar,beta_m_bar,Fs,pulseDuration);
            h_temp = stateVariableDynamics_2(h(i-1),alpha_h_bar,beta_h_bar,Fs,pulseDuration);
            n_temp = stateVariableDynamics_1(n(i-1),alpha_n_bar,beta_n_bar,Fs,pulseDuration);
            q_temp = stateVariableDynamics_1(q(i-1),alpha_q_bar,beta_q_bar,Fs,pulseDuration);
            m(i:i+length(m_temp)-1) = m_temp;
            h(i:i+length(h_temp)-1) = h_temp;
            n(i:i+length(n_temp)-1) = n_temp;  
            q(i:i+length(q_temp)-1) = q_temp;  
        end
    end
    
    I_ion = g_Na_bar*m(i)^3*h(i)*(V_s - E_Na) + g_Kf_bar*n(i)^4*(V_s - E_Kf) + g_Ks_bar*q(i)^2*(V_s-E_Ks);
    % Membrane current due to the voltage-dependent ionic conductances
    
    dV_d = (-I_syn_d(i) - g_ld*(V_d - E_l) - g_c*(V_d - V_s) + I_inj_d(i))/C_d; % dV in volts
    dV_d = dV_d*1000;
    % dendritic membrane potential
    V_d = dV_d*1/Fs + V_d;
    dV_s = (-I_syn_s(i) - g_ls*(V_s - E_l) - g_c*(V_s - V_d) - I_ion + I_inj_s(i))/C_s; % dV in volts
    dV_s = dV_s*1000;
    % somatic membrane potential
    V_s = dV_s*1/Fs + V_s;
    
    
    V_d_vec(i) = V_d;
    V_s_vec(i) = V_s;
    I_ion_vec(i) = I_ion;
    
end

output.V_d = V_d_vec;
output.V_s = V_s_vec;
output.I_ion = I_ion_vec;
output.m = m;
output.h = h;
output.n = n;
output.q = q;


function y = stateVariableDynamics_1(x0,alpha,beta,Fs,pulseDuration)
pulseDuration_step = pulseDuration*Fs/1000;
t_temp = 0:1/Fs:0.5;

y1 = 1 + (x0 - 1)*(exp(-alpha*1000.*t_temp));
y2 = y1(pulseDuration_step)*exp(-beta*1000.*t_temp);
y = [y1(1:pulseDuration_step) y2(2:end)];

end

function y = stateVariableDynamics_2(x0,alpha,beta,Fs,pulseDuration)
pulseDuration_step = pulseDuration*Fs/1000;
t_temp = 0:1/Fs:0.5;

y1 = x0*exp(-beta*1000.*t_temp);
y2 = 1 + (y1(pulseDuration_step)-1)*(exp(-alpha*1000.*t_temp));
y = [y1(1:pulseDuration_step) y2(2:end)];
end

end