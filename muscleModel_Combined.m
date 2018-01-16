%--------------------------------------------------------------------------
% muscleModel_Combined.m
% Author: Akira Nagamori
% Last update: 12/26/207
%--------------------------------------------------------------------------

function output = muscleModel_Combined(t,Fs,input,modelParameter,simulationParameter,recruitmentType)
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
%--------------------------------------------------------------------------
cv_MU = 0; %ISI variability as per coefficient of variation (=mean/SD)
U = input;
U_eff = 0;
T_U = 0.03;
%--------------------------------------------------------------------------
% initialize parameters
Y_af_temp = 0;
S_af_temp = zeros(N_MU,1);
Lce = simulationParameter.Lce;
Vce = 0;

spike_time = zeros(N_MU,1);
spike_train = zeros(N_MU,length(t));
force = zeros(N_MU,length(t));
FR_MU = zeros(N_MU,length(t));
S_af = zeros(N_MU,length(t));
Y_af = zeros(N_MU,length(t));
Af = zeros(N_MU,1);
FF = zeros(N_MU,1);

Outputfenv = zeros(N_MU,length(t));
Outputforce_half= zeros(N_MU,length(t));
OutputAf = zeros(N_MU,length(t));
OutputFF = zeros(N_MU,length(t));

Fs_MU = 1000;
count = 1;
index_MU = 1:Fs/Fs_MU:length(t);

%--------------------------------------------------------------------------
% Simulation
for i = 1:length(t)
    U_eff_dot = (U(i) - U_eff)/T_U;
    U_eff = U_eff_dot*1/Fs + U_eff;
    %----------------------------------------------------------------------
    % Determine firing rate of individual motor units given an input (U)
    if recruitmentType == 1 % linear increase in firing rate up to Ur
        FR_MU(:,i) = (PFR_MU-MFR_MU)./(1-U_th).*(U_eff-U_th) + MFR_MU;
    elseif recruitmentType == 2 % equal gain across units and saturation 
        FR_MU(:,i) = g_e*(U(i)-U_th)+ MFR_MU;
    end
    % Make FR of units that is below 8 Hz be zero
    FR_MU(FR_MU<8) = 0;
    % Determine f_env
    f_env = FR_MU(:,i)./FR_half';
    % Determine force_half
    force_half = force(:,i)./Pi_half';
    % Store
    Outputfenv(:,i) = f_env;
    Outputforce_half(:,i) = force_half;
    % Yielding
    Y_af_temp = yield_function(Y_af_temp,Vce,Fs);
    Y_af(1:index_slow,i) = Y_af_temp;
    for f_temp = 1:N_MU
        if FR_MU(f_temp,:) > PFR_MU(f_temp)
            FR_MU(f_temp,:) = PFR_MU(f_temp);
        end
        if f_temp <= index_slow
            Af(f_temp) = Af_slow_function(force_half(f_temp),Lce,Y_af(f_temp,i));
            OutputAf(f_temp,i)= Af(f_temp);
            FF(f_temp) = frequency2Force_slow_function(f_env(f_temp),Lce,Y_af(f_temp,i));
            OutputFF(f_temp,i)= FF(f_temp);
        else
            % Sag
            S_af_temp(f_temp) = sag_function(S_af_temp(f_temp),force_half(f_temp),Fs);
            S_af(f_temp,i) = S_af_temp(f_temp);
            
            Af(f_temp) = Af_fast_function(force_half(f_temp),Lce,S_af(f_temp,i));
            OutputAf(f_temp,i)= Af(f_temp);
            FF(f_temp) = frequency2Force_fast_function(f_env(f_temp),Lce,S_af(f_temp,i));
            OutputFF(f_temp,i)= FF(f_temp);
        end
    end
    % Af
    if i == count       
        for n = 1:length(find(FR_MU(:,i)>=MFR_MU)) % loop through motor units whose firing rate is greater than minimum firing rate defined by the user                     
            spike_train_temp = zeros(1,length(t));
            if ~any(spike_train(n,:)) % when the motor unit fires at the first time
                spike_train(n,i) = 1; % add a spike to the vector
                spike_train_temp(i) = 1;
                mu = 1/FR_MU(n,i);
                Z = randn(1);
                Z(Z>3.9) = 3.9;
                Z(Z<-3.9) = -3.9;
                spike_time_temp = (mu + mu*cv_MU*Z)*Fs;                               
                spike_time_temp2 = round(spike_time_temp) + i;
                spike_time(n) = interp1(index_MU,index_MU,spike_time_temp2,'nearest');
                
                [twitch_temp,~,~] = twitch_function(Af(n),Lce,CT(n),RT(n),Fs);
                twitch =  Pi(n).*twitch_temp*FF(n);
                force_temp = conv(spike_train_temp,twitch);
                force(n,:) = force(n,:) + force_temp(1:length(t));
            else
                if spike_time(n) == i % when the motor unit fires
                    spike_train(n,i) = 1;
                    spike_train_temp(i) = 1;
                    % update mean firing rate of the motor unit given the
                    % current value of input
                    mu = 1/FR_MU(n,i); % interspike interval
                    Z = randn(1);
                    Z(Z>3.9) = 3.9;
                    Z(Z<-3.9) = -3.9;
                    spike_time_temp = (mu + mu*cv_MU*Z)*Fs; % interspike interval
                    spike_time_temp2 = round(spike_time_temp) + i;
                    spike_time(n) = interp1(index_MU,index_MU,spike_time_temp2,'nearest');
                    
                    [twitch_temp,~,~] = twitch_function(Af(n),Lce,CT(n),RT(n),Fs);
                    twitch =  Pi(n).*twitch_temp*FF(n);
                    force_temp = conv(spike_train_temp,twitch);
                    force(n,:) = force(n,:) + force_temp(1:length(t));
                elseif i > spike_time(n) + round(1/FR_MU(n,i)*Fs)
                    spike_train(n,i) = 1;
                    spike_train_temp(i) = 1;
                    spike_time(n) = i;
                    mu = 1/FR_MU(n,i); % interspike interval
                    Z = randn(1);
                    Z(Z>3.9) = 3.9;
                    Z(Z<-3.9) = -3.9;
                    spike_time_temp = (mu + mu*cv_MU*Z)*Fs; % interspike interval
                    spike_time_temp2 = round(spike_time_temp) + i;
                    spike_time(n) = interp1(index_MU,index_MU,spike_time_temp2,'nearest');
                    
                    % force_temp = conv(spike_train_temp,g.*twitch(n,:)).*FL_temp.*FV_temp.*Af;
                    [twitch_temp,~,~] = twitch_function(Af(n),Lce,CT(n),RT(n),Fs);
                    twitch =  Pi(n).*twitch_temp*FF(n);
                    force_temp = conv(spike_train_temp,twitch);
                    force(n,:) = (force(n,:) + force_temp(1:length(t)));
                end
            end
            if n <= index_slow
                FL_temp = FL_slow_function(Lce);
                if Vce > 0
                    FV_temp = FVecc_function(Lce,Vce);
                else
                    FV_temp = FVcon_function(Lce,Vce);
                end
            else
                FL_temp = FL_fast_function(Lce);
                if Vce > 0
                    FV_temp = FVecc_fast_function(Lce,Vce);
                else
                    FV_temp = FVcon_fast_function(Lce,Vce);
                end
            end
            force(n,i) = force(n,i)*FL_temp*FV_temp;
        end
        
        count = count + Fs/Fs_MU;
    end
end

Force = sum(force);
output.Force_total = Force;
output.force = force;
output.FR = FR_MU;
output.spike_train = spike_train;
output.f_env = Outputfenv;
output.force_half = Outputforce_half;
output.Af = OutputAf;
output.FF = OutputFF;
output.S = S_af;
output.Y = Y_af;


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
        Af = 1 - exp(-(S*f_eff/(a_f*n_f))^n_f);
    end
    
    function FF = frequency2Force_slow_function(f_env,L,Y)
        a_f = 0.52;
        n_f0 = 1.97;
        n_f1 = 5.1;
        n_f = n_f0 + n_f1*(1/L-1);
        FF = 1 - exp(-(Y*f_env/(a_f*n_f))^n_f);
        if f_env == 0
            FF = 0;
        else
            FF = FF/f_env;
        end
    end

    function FF = frequency2Force_fast_function(f_env,L,S)
        a_f = 0.52;
        n_f0 = 1.97;
        n_f1 = 3.28;
        n_f = n_f0 + n_f1*(1/L-1);
        FF = 1 - exp(-(S*f_env/(a_f*n_f))^n_f);
        if f_env == 0
            FF = 0;
        else
            FF = FF/f_env;
        end
    end

    function [twitch,T1,T2_temp] = twitch_function(Af,Lce,CT,RT,Fs)
        T1 = CT*Lce^2+(CT*1/2)*Af;
        T2_temp = (RT + (RT*1/2)*Af)/Lce;
        T2 = T2_temp/1.68;
        t_twitch = 0:1/Fs:2;
        f_1 = t_twitch./T1.*exp(1-t_twitch./T1);
        f_2 = t_twitch./T2.*exp(1-t_twitch./T2);
        
        twitch = [f_1(1:round(T1*Fs+1)) f_2(round(T2*Fs+1):end)];
        twitch = twitch(1:1*Fs);
        
    end



    function FL = FL_slow_function(L)
        %---------------------------
        % force length (F-L) relationship for slow-tiwtch fiber
        % input: normalized muscle length and velocity
        % output: F-L factor (0-1)
        %---------------------------
        beta = 2.3;
        omega = 1.12;
        rho = 1.62;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
    end

    function FL = FL_fast_function(L)
        %---------------------------
        % force length (F-L) relationship for fast-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-L factor (0-1)
        %---------------------------
        beta = 1.55;
        omega = 0.75;
        rho = 2.12;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
    end

    function FVcon = FVcon_function(L,V)
        %---------------------------
        % concentric force velocity (F-V) relationship for slow-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        Vmax = -7.88;
        cv0 = 5.88;
        cv1 = 0;
        
        FVcon = (Vmax - V)/(Vmax + (cv0 + cv1*L)*V);
    end

    function FVcon = FVcon_fast_function(L,V)
        %---------------------------
        % concentric force velocity (F-V) relationship for fast-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        Vmax = -9.15;
        cv0 = -5.7;
        cv1 = 9.18;
        
        FVcon = (Vmax - V)/(Vmax + (cv0 + cv1*L)*V);
    end

    function FVecc = FVecc_function(L,V)
        %---------------------------
        % eccentric force velocity (F-V) relationship for slow-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        av0 = -4.7;
        av1 = 8.41;
        av2 = -5.34;
        bv = 0.35;
        FVecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end

    function FVecc = FVecc_fast_function(L,V)
        %---------------------------
        % eccentric force velocity (F-V) relationship for fast-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        av0 = -1.53;
        av1 = 0;
        av2 = 0;
        bv = 0.69;
        FVecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end

end
