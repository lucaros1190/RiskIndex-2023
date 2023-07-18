
% Plot of the Risk Index for D. suzukii over temperature
% Developed by Luca Rossini on 09/03/2023
% e-mail: luca.rossini@unitus.it

% Parameters environment

    % Briere: Egg-Pupa stages

a_EP = 1.5909 * 10^(-4);
T_L_EP = 2.0919;
T_M_EP = 32.088;
m_EP = 4.0;

    % Briere: Pupa-Adult stage

a_PA = 2.3699 * 10^(-4);
T_L_PA = 4.0;
T_M_PA = 33.164;
m_PA = 4.0;

    % Briere: Adult survival

a_Sur = 6.8417 * 10^(-5);
T_L_Sur = -3.0;
T_M_Sur = 30.034;
m_Sur = 2.50;

    % Fertility

alpha = 659.06;
gamma = 88.53;
lambda = 52.32;
delta = 6.06;
tau = 22.87;
T_low = -3;
T_max = 39;

    % Mortality

k_mort = 1.0;
T_MAX_mort = 23.4265085; 
rho_mort = -5.5455493;

    % Sex ratio

SR = 0.5;

% Calculate the values of the RiskIndex to plot

    % Define the array to store results
RI = [];
threshold = [];
temp = [];

    % The for loop to evaluate the RiskIndex

for j = 1 : 4000
    
    % Development rates - The if/else statements are to avoid complex
    % numbers in the RiskIndex!

    i=j*0.01;
    
    GeL = a_EP * i * (i - T_L_EP) * (T_M_EP - i)^(1 / m_EP);

    if (i > T_L_EP && i < T_M_EP)
        GeL = a_EP * i * (i - T_L_EP) * (T_M_EP - i)^(1 / m_EP);
    else
        GeL = 0.001;
    end


    GP = a_PA * i * (i - T_L_PA) * (T_M_PA - i)^(1 / m_PA);

    if (i > T_L_PA && i < T_M_PA)
        GP = a_PA * i * (i - T_L_PA) * (T_M_PA - i)^(1 / m_PA);
    else
        GP = 0.001;
    end

    GA = a_Sur * i * (i - T_L_Sur) * (T_M_Sur - i)^(1 / m_Sur);

    if (i > T_L_Sur && i < T_M_Sur)
        GA = a_Sur * i * (i - T_L_Sur) * (T_M_Sur - i)^(1 / m_Sur);
    else
        GA = 0.001;
    end
    

    % Fertility rate

    beta = alpha * (((gamma + 1) / (pi * (lambda^(2 * gamma + 2)))) * ((lambda^2) ...
        - ( ((i - tau)^2) + (delta^2)) )^gamma);
    
    % Mortality rate

    Survival = k_mort * exp(1 + ((T_MAX_mort - i) / rho_mort) - ...
        exp((T_MAX_mort - i) / rho_mort));

    M = 1 - Survival;
    

    RiskIndex = -(beta * (GA + M - 1) * SR * GP * GeL^4) / ((GeL + M)^4 * ...
        (GP + M) * (GA + M));

    RI = [RI, RiskIndex];
    threshold = [threshold, 1];
    temp = [temp, i];


end

figure

hold on 
plot(temp, RI, LineWidth=2)
plot(temp, threshold, 'Color', 'red', 'LineStyle', '--', LineWidth=2)
legend('Risk Index', 'Threshold for stability')
title('Risk Index for Drosophila suzukii')
xlabel('Temperature (°C)')
ylabel('Risk Index')
