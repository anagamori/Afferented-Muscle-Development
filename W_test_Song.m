close all
clear all
clc


U_eff = 0:0.01:1;
U1_th = 0.001;
U2_th = 0.32;

for i = 1:length(U_eff)
    if U_eff(i) < U1_th
        W1(i) = 0;
    elseif U_eff(i) < U2_th
        W1(i) = (U_eff(i) - U1_th)/(U_eff(i) - U1_th);
    else
        W1(i) = (U_eff(i) - U1_th)/((U_eff(i) - U1_th) + (U_eff(i) - U2_th));
    end
    if U_eff(i) < U2_th
        W2(i) = 0;
    else
        W2(i) = (U_eff(i) - U2_th)/((U_eff(i) - U1_th) + (U_eff(i) - U2_th));
    end
    
    Act(i) = U_eff(i)*(W1(i)+W2(i));
end

figure(1)
subplot(3,1,1)
plot(U_eff,W1)
subplot(3,1,2)
plot(U_eff,W2)
subplot(3,1,3)
plot(U_eff,Act)