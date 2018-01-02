close all
clear all
clc

Fpcsa_slow_temp = [0.2 0.4 0.5 0.6 0.8];
for i = 2 %:length(Fpcsa_slow_temp)
    Fpcsa_slow = Fpcsa_slow_temp(i);
    Fpcsa_fast = 1-Fpcsa_slow;
    U_r = 0.8;
    U1_th = 0.01;
    U2_th = U_r*Fpcsa_slow;
    
    W1_slow_temp_1 = 1.56*U2_th^2-1.2*U2_th+0.884;
    
    x = [U2_th 0.4 0.6 1]'
    y = [W1_slow_temp_1 0.75 0.9 1]'
    
    eq = ['a*(1-exp(-(x-' num2str(U2_th) ')/b))'];
    
    f = fit(x,y,eq)
    
    
    coeffvals = coeffvalues(f);
    xx = U2_th:0.01:1;
    yy = W1_slow_temp_1+coeffvals(1)*(1-exp(-(xx-U2_th)/coeffvals(2)));
    
    figure(1)
    plot(xx,yy)
    hold on
    
    A(i) = coeffvals(1);
    tau(i) = coeffvals(2);
    U2_th_vec(i) = U2_th;
   
end

% f_A = fit(U2_th_vec',A','poly4');
% 
% f_tau = fit(U2_th_vec',tau','poly1');
