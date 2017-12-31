close all
%clear all
clc

U_eff = 0:0.01:1;

Fpcsa_slow = 0.4;
Fpcsa_fast = 1-Fpcsa_slow;
U_r = 0.8;
U1_th = 0.01;
U2_th = U_r*Fpcsa_slow;


Y = 1;
S = 0.96;
Lce = 1.1;

for i = 1:length(U_eff)
    if U_eff(i) >=  U1_th
            f_env_slow = 1.5/(1-U1_th).*(U_eff(i)-U1_th)+0.5;
           
        else
            f_env_slow = 0;
    end
    if U_eff(i) >= U2_th
        f_env_fast = 1.5/(1-U2_th).*(U_eff(i)-U2_th)+0.5;
       
    else
        f_env_fast = 0;
    end
    Af_slow(i) = Af_slow_function(f_env_slow,Lce,Y);
    Af_fast(i) = Af_fast_function(f_env_fast,Lce,S);
    
    if U_eff(i) < U1_th
        rFpcsa_slow(i) = 0;
        W1(i) = 0;
    elseif U_eff(i) >= U1_th &&  U_eff(i) < U2_th
        rFpcsa_slow(i) = Fpcsa_slow/(1-U1_th)*(U_eff(i)-U1_th);
        W1(i) = 1.56*U_eff(i)^2-1.2*U_eff(i)+0.884;
    elseif U_eff(i) >= U2_th
        rFpcsa_slow(i) = Fpcsa_slow/(1-U1_th)*(U_eff(i)-U1_th);
        W1_slow_temp_1 = 1.56*U2_th^2-1.2*U2_th+0.884;
        A = -1.99*U2_th^4+1.899*U2_th^3-1.748*U2_th^2+1.186*U2_th+0.1298;
        tau_s = 0.26;
        W1_slow_temp_2 = A*(1-exp(-(U_eff(i)-U2_th)/tau_s));
        W1(i) = W1_slow_temp_1 + W1_slow_temp_2;
    end
    if U_eff(i) < U2_th
        rFpcsa_fast(i) = 0;
        W2(i) = 0;
    elseif U_eff(i) >= U2_th && U_eff(i) < U_r       
        rFpcsa_fast(i) = Fpcsa_fast/(1-U2_th)*(U_eff(i)-U2_th);
        B = 0.59*U2_th + 0.39;
        C = -98.82*U2_th^5+155.7*U2_th^4-85.96*U2_th^3+19.32*U2_th^2-2.62*U2_th+1.02;
        W2(i) = ((U_eff(i)-B)/C)^2+0.65;
    elseif U_eff(i) >= U_r
        rFpcsa_fast(i) = Fpcsa_fast/(1-U2_th)*(U_eff(i)-U2_th);
        B = 0.59*U2_th + 0.39;
        C = -98.82*U2_th^5+155.7*U2_th^4-85.96*U2_th^3+19.32*U2_th^2-2.62*U2_th+1.02;
        W2_temp_1 = ((U_r-B)/C)^2+0.65;
        D = 0.35;
        tau_f = -0.27*U2_th+0.25;
        W2_temp_2 = D*(1-exp(-(U_eff(i)-U_r)/tau_f));
        W2(i) = W2_temp_1 + W2_temp_2;
    end
    
    Act(i) = rFpcsa_slow(i)*W1(i)*Af_slow(i)+rFpcsa_fast(i)*W2(i)*Af_fast(i);
end

figure(1)
subplot(2,1,1)
plot(U_eff,rFpcsa_slow)
subplot(2,1,2)
plot(U_eff,rFpcsa_fast)

figure(2)
subplot(2,1,1)
plot(U_eff,W1)
ylim([0.5 1])
subplot(2,1,2)
plot(U_eff,W2)
ylim([0.5 1])

figure(3)
subplot(2,1,1)
plot(U_eff,Af_slow)
subplot(2,1,2)
plot(U_eff,Af_fast)

figure(4)
plot(U_eff,Act./Act(end))


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