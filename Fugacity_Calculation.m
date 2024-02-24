%% Defining Constants

close all; clear all; clc;
warning off;

% making variables needed for EoS global so they can be accessed 
% from within the function at the bottom of the code
global P1 T1 structure amix bmix R a1 a2 a12 b1 b2 y1_glob y2_glob


% general parameters/constants; number 1 assigned to CO, 2 assigned to H2
R = 8.314; % J/mol*K
Tc1 = 132.9; Tc2 = 33.2; % in terms of Kelvin
Pc1 = 34.96*10^5; Pc2 = 12.97*10^5; % in terms of Pascal
n_total = 10; % assume basis for total moles for later calculations

% Van der Waals parameters properties
a1 = (27/64)*(((R*Tc1)^2)/(Pc1)); % a_CO: corrective term for CO molecules 
% attractive forces
b1 = (R*Tc1)/(8*Pc1); % b_CO: corrective constant for CO molecule volumes

a2 = (27/64)*(((R*Tc2)^2)/(Pc2)); % a_H2: corrective term for H2 molecules 
% attractive forces
b2 = (R*Tc2)/(8*Pc2); % b_H2: corrective constant for H2 molecule volumes

a12 = ((a1*a2)^1/2); % intermolecular constant accounting for the 
% attractive forces of CO on H2 and vice-versa


% creating range of composition, pressure, and temperature values to test
y1_glob = linspace(0,1,50);
y2_glob = linspace(1,0,50);
P1 = linspace(12000,10000000,50);
T1 = linspace(10,650,50);

% setting variables for the length of my values I'm looking to test;
% not necessary but makes further calculations easy to update if I decide
% to change the range of values I'm looking to solve for
P1_length = length(P1);
T1_length = length(T1);
y1_length = length(y1_glob);

% creating relative pressure values for all pressure values being tested
for k = 1:length(P1)
    Pr_CO(k) = (P1(k)/Pc1);
    Pr_H2(k) = (P1(k)/Pc2);
end

% creating relative temperature values for all pressure values being tested
for s = 1:length(T1)
    Tr_CO(s) = T1(s)/Tc1;
    Tr_H2(s) = T1(s)/Tc2;
end

%% General VdW Structure for Organization


% creating structure organized into composition, pressure, temperature,
% amix, bmix, and molar volume to compare later 
% This can later be used to tie  together relations between compositions, 
% pressures, and temperatures--do this by taking one of the three as 
% variable at multiple combinations of constant values for the other 2

for k = 1:length(y1_glob)
    % % composition loop
    structure.y1(k) = (y1_glob(k));
    structure.y2(k) = (1-y1_glob(k));
    y1 = y1_glob(k);
    y2 = y2_glob(k);
    
    % mixing properties at all composition values
    amix = (y1^2)*a1 + (y2^2)*a2 + 2*y1*y2*a12;
    bmix = y1*b1 + y2*b2;

    % storing amix and bmix data into a structure per each composition
    structure.amix(k) = amix;
    structure.bmix(k) = bmix;
end


%% Solving for Volumes at These Defined Ranges


% Using van der waals to solve for intensive volumes of CO and H2 
% at all combinations of T & P (does not change with molar composition)

% Since 50 T and 50 P values, should give a 50x50 matrix with 2500
% combinations of T, P, and v values
vol_CO = zeros([50 50]);
vol_H2 = zeros([50 50]);

% initial condition for the fsolve function used to solve for molar volume
x0 = rand(1,1);

% loop solving for the molar volume at each T and P combination (2500
% values) using the van der waals equation of state and fsolve

for i = 1:length(P1)
    % % loop running through each pressure value
    P = P1(i);

    % storing pressure values in structure
    structure.P(i) = P1(i);
    
    for j = 1:length(T1)
        % % loop running through each temperature value with respect to
        % % each pressure
        T = T1(j);

        % storing temperature values in structure
        structure.T(j) = T1(j);

        % setting expression of van der waals equal to 0 with molar volume
        % defined as the constant (one expression for each component)
        f_CO = @(v1) R*T/(v1-b1) - a1/(v1^2) - P;
        f_H2 = @(v2) R*T/(v2-b2) - a2/(v2^2) - P;


        % using fsolve as an equation solver for the above expressions
        % should output molar volumes for each (T, P) state
        [v_CO, r_CO, exitf] = fsolve(f_CO,x0);
        [v_H2, r_H2, exitf] = fsolve(f_H2,x0);

        % storing molar volumes in structure
        structure.vol_CO(i,j) = v_CO;
        structure.vol_H2(i,j) = v_H2;

        % molar volume only changes with temperature and pressure (seen
        % by the fact that its only variables are T, P, a1, and b1);
        % therefore it does not need to be sorted by composition
    end
    
end

%% Using Function to Find Fugacity-Related Values with Previously-Found Volumes


% function that outputs structure containing ln(phi), log(phi), phi, 
% and fugacity of CO and H2 respectively;
% details of the function are listed within the function itself 
% at the bottom of the code 
[values] = EoS (P1_length, T1_length, y1_length);

%% Organizing Lee-Kesler Appendix Values


% % Using Lee-Kessler correlations to solve for 2 different points.
% % The first point is close to ideal, for CO: 
% % P1_ideal = P1(6) = 10.312 * 10^5 Pa; T1_ideal = T1(22) = 284.286 K
% % for H2: P1_ideal = P1(2) = 2.15 * 10^5 Pa; T1_ideal = T1(10) = 127.5
% % The next point does not behave ideally: P2 = P1(26) = 5.1079 * 10^6 Pa
% % T2 = T1(6) = 75.31 K;

P1CO_ideal = P1(6); P1H2_ideal = P1(2); T1CO_ideal = T1(22); T1H2_ideal = T1(10);
P2_nonideal = P1(26); T2_nonideal = T1(6);

% Relative P and T for CO and H2 at ideal point
Pr_LK_CO1 = P1CO_ideal/Pc1; Tr_LK_CO1 = T1CO_ideal/Tc1;
Pr_LK_H21 = P1H2_ideal/Pc2; Tr_LK_H21 = T1H2_ideal/Tc2;

% Relative P and T for CO and H2 at non;-ideal point
Pr_LK_CO2 = P2_nonideal/Pc1; Tr_LK_CO2 = T2_nonideal/Tc1;
Pr_LK_H22 = P2_nonideal/Pc2; Tr_LK_H22 = T2_nonideal/Tc2;

% Creating Lee Kessler structure to organize values of appendix:

global LK

% general structure: Lee_Kess.CO.(Pr, Tr, log(phi)0, log(phi)1)
% general structure: Lee_Kess.H2.(Pr, Tr, log(phi)0, log(phi)1)
% each value in the parentheses will be a separate array of different 
% values within the specific molecules category

% Pr organized into: 
% first row is the ideal point, and the 2nd row is the second point
% columns 1-3 are from the C.7 Appendix: P_low, Pr, P_high
% since C.7 and C.8 have the same intervals and both , I will use the same index
% values for interpolation for both

% relevant relative pressure and temperature for interpolation and 
% solving for the log(phi) of each component using the Lee-Kesler 
% Equation of State

LK.CO.Pr = [0.25, Pr_LK_CO1, 0.5; % ideal Pr of CO: P_low, Pr, P_high
            1.4, Pr_LK_CO2, 1.5]; % second Pr of CO: P_low, Pr, P_high
LK.CO.Tr = [2, Tr_LK_CO1, 2.25; % ideal Tr of CO: T_low, Tr, T_high
            0.55, Tr_LK_CO2, 0.6]; % second Tr of CO: T_low, Tr, T_high
LK.H2.Pr = [0.1, Pr_LK_H21, 0.25; % ideal Pr of H2: P_low, Pr, P_high
            3, Pr_LK_H22, 4]; % second Pr of H2: P_low, Pr, P_high
LK.H2.Tr = [3.5, Tr_LK_H21, 4; % ideal Tr of H2: T_low, Tr, T_high
            2.25, Tr_LK_H22, 2.5]; % second Tr of H2: T_low, Tr, T_high

% matching log(phi)^0, and log(phi)^1 values to each respective combination
% of Pr and Tr
% In order to interpolate for the log(phi) at each value, I only need the
% respective values from C.7 and C.8 of: 
% log(phi) of [(P_low, T_low), (P_low, T_high), and (P_high, T_high)]


% first row is ideal, second row is the second point
% columns 1-4 are log(phi) of:
% (P_low, T_low), (P_low, T_high), (P_high, T_low), and (P_high, T_high) from C.7
% columns 5-8 are log(phi) of:
% (P_low, T_low), (P_low, T_high), (P_high, T_low), and (P_high, T_high) from C.8
LK.CO.AppC = [-0.003, -0.002, -0.006, -0.003,      0.008, 0.008, 0.017, 0.016;
              -1.947, -1.608, -1.969, -1.630,      -2.214, -1.692, -2.218, -1.695];

LK.H2.AppC = [0, 0, 0.001, 0.001,             0.002, 0.002, 0.006, 0.005;
              -0.012, -0.002, -0.013, 0,      0.9, 0.083, 0.116, 0.108];

%% Interpolations to Find log(phi)^0 and log(phi)^1 at Specified Points


% interpolating all data for the Lee_Kesler function

% y1 will be the first interpolation of the set, following this pattern:
% Interpolate between (P_low, T_low) and (P_low, T_high) to find 
% log(phi) of (P_low, Tr); 
% use interpolation function defined at bottom by setting:
% yl to log_phi (P_low, T_low), yh to log_phi (P_low, T_high), 
% xl = T_low, x = Tr at desired point, and xh = T_high

% for C.7
y1(1,1) = interp(LK.CO.AppC(1,1),LK.CO.AppC(1,2),LK.CO.Tr(1,1), ...
          LK.CO.Tr(1,2),LK.CO.Tr(1,3)); % CO_ideal
y1(1,2) = interp(LK.CO.AppC(2,1),LK.CO.AppC(2,2),LK.CO.Tr(2,1), ...
          LK.CO.Tr(2,2),LK.CO.Tr(2,3)); % CO_second
y1(1,3) = interp(LK.H2.AppC(1,1),LK.H2.AppC(1,2),LK.H2.Tr(1,1), ...
          LK.H2.Tr(1,2),LK.H2.Tr(1,3)); % H2_ideal
y1(1,4) = interp(LK.H2.AppC(2,1),LK.H2.AppC(2,2),LK.H2.Tr(2,1), ...
          LK.H2.Tr(2,2),LK.H2.Tr(2,3)); % H2_second

% for C.8
y1(2,1) = interp(LK.CO.AppC(1,4),LK.CO.AppC(1,5),LK.CO.Tr(1,1), ...
          LK.CO.Tr(1,2),LK.CO.Tr(1,3)); % CO_ideal
y1(2,2) = interp(LK.CO.AppC(2,4),LK.CO.AppC(2,5),LK.CO.Tr(2,1), ...
          LK.CO.Tr(2,2),LK.CO.Tr(2,3)); % CO_second
y1(2,3) = interp(LK.H2.AppC(1,4),LK.H2.AppC(1,5),LK.H2.Tr(1,1), ...
          LK.H2.Tr(1,2),LK.H2.Tr(1,3)); % H2_ideal
y1(2,4) = interp(LK.H2.AppC(2,4),LK.H2.AppC(2,5),LK.H2.Tr(2,1), ...
          LK.H2.Tr(2,2),LK.H2.Tr(2,3)); % H2_second


% y2 will be the second interpolation of the set, following this pattern:
% Interpolate between (P_high, T_low) and (P_high, T_high) to find 
% log(phi) of (P_high, Tr); yl = log_phi (P_high, T_low) 
% yh = log_phi (P_high, T_high), xl = T_low, x = Tr at desired point, 
% and xh = T_high

% for C.7
y2(1,1) = interp(LK.CO.AppC(1,3),LK.CO.AppC(1,4),LK.CO.Tr(1,1), ...
          LK.CO.Tr(1,2),LK.CO.Tr(1,3)); % CO_ideal
y2(1,2) = interp(LK.CO.AppC(2,3),LK.CO.AppC(2,4),LK.CO.Tr(2,1), ...
          LK.CO.Tr(2,2),LK.CO.Tr(2,3)); % CO_second
y2(1,3) = interp(LK.H2.AppC(1,3),LK.H2.AppC(1,4),LK.H2.Tr(1,1), ...
          LK.H2.Tr(1,2),LK.H2.Tr(1,3)); % H2_ideal
y2(1,4) = interp(LK.H2.AppC(2,3),LK.H2.AppC(2,4),LK.H2.Tr(2,1), ...
          LK.H2.Tr(2,2),LK.H2.Tr(2,3)); % H2_second

% for C.8
y2(2,1) = interp(LK.CO.AppC(1,7),LK.CO.AppC(1,8),LK.CO.Tr(1,1), ...
          LK.CO.Tr(1,2),LK.CO.Tr(1,3)); % CO_ideal
y2(2,2) = interp(LK.CO.AppC(2,7),LK.CO.AppC(2,7),LK.CO.Tr(2,1), ...
          LK.CO.Tr(2,2),LK.CO.Tr(2,3)); % CO_second
y2(2,3) = interp(LK.H2.AppC(1,7),LK.H2.AppC(1,7),LK.H2.Tr(1,1), ...
          LK.H2.Tr(1,2),LK.H2.Tr(1,3)); % H2_ideal
y2(2,4) = interp(LK.H2.AppC(2,7),LK.H2.AppC(2,8),LK.H2.Tr(2,1), ...
          LK.H2.Tr(2,2),LK.H2.Tr(2,3)); % H2_second


%% Lee-Kesler Function to Find log(phi) of CO and H2 at 2 Different Points


% using Lee-Kesler function at bottom of code to solve for the log(phi) of
% CO and H2 at an ideal point and a secondary point

[LK.CO.logphi, LK.H2.logphi] = Lee_Kesler(y1,y2);

% Since fugacity is the product of the fugacity coefficient and pressure,
% raise 10 to a factor of each log(phi) value (isolating each
% phi/coefficient), then multiply the result by the pressure used for that
% specific Lee-Kesler correlation.

LK.CO.fugacity(1) = (10^(LK.CO.logphi(1)))*P1(6);
% Lee Kesler fugacity at ideal (P = 10.312 * 10^5 Pa; T = 284.286 K)
LK.CO.fugacity(2) = (10^(LK.CO.logphi(2)))*P1(26);
% Lee Kesler fugacity at ideal (P = 5.1079 * 10^6 Pa; T = 75.31 K)
LK.H2.fugacity(1) = (10^(LK.H2.logphi(1)))*P1(2);
% Lee Kesler fugacity at ideal (P = 2.15 * 10^5 Pa; T = 127.5 K)
LK.H2.fugacity(2) = (10^(LK.H2.logphi(2)))*P1(26);
% Lee Kesler fugacity at ideal (P = 5.1079 * 10^6 Pa; T = 75.31 K)


%% Legend Labels


% defining names/labels to be used next for the legends of the plots

for i = 1:P1_length
legend_text_PrCO(i) = {sprintf("Pr = %6.2f",Pr_CO(i))};
legend_text_PrH2(i) = {sprintf("Pr = %6.2f",Pr_H2(i))};
legend_text_TrCO(i) = {sprintf("Tr = %6.2f",Tr_CO(i))};
legend_text_TrH2(i) = {sprintf("Tr = %6.2f",Tr_H2(i))};
legend_text_y1(i) = {sprintf("y1 = %6.2f",y1_glob(i))};
legend_text_y2(i) ={sprintf("y2 = %6.2f",y2_glob(i))};
end

% labels for the legends of Lee-Kessler values
legend_text_PrCO(P1_length+1) = {'Lee Kessler: Pr = 0.2950; Tr = 2.1391'};
legend_text_PrCO(P1_length+2) = {'Lee Kessler: Pr = 1.4611; Tr = 0.5666'};

legend_text_PrH2(P1_length+1) = {'Lee Kessler: Pr = 0.1664; Tr = 3.84'};
legend_text_PrH2(P1_length+2) = {'Lee Kessler: Pr = 3.9383; Tr = 2.2683'};

legend_text_TrCO(P1_length+1) = {'Lee Kessler: Pr = 0.2950; Tr = 2.1391'};
legend_text_TrCO(P1_length+2) = {'Lee Kessler: Pr = 1.4611; Tr = 0.5666'};

legend_text_TrH2(P1_length+1) = {'Lee Kessler: Pr = 0.1664; Tr = 3.84'};
legend_text_TrH2(P1_length+2) = {'Lee Kessler: Pr = 3.9383; Tr = 2.2683'};


%% Plotting log(phi) vs Reduced Pressure


% Plotting log(phi) vs changing pressures (constant T and composition)
% using 13 isotherms at y1(16) = 0.3061 and y1(34) = 0.6939

figure(1)
hold on;
sgtitle('Changing Pressures (Constant Temperature and Composition)')

for n = 2:4:P1_length
    
    % first tile of changing pressure; plots the change in the fugacity
    % coefficient of CO as a function of pressure at y1_CO = 0.3061 across
    % 10 different Temperatures. 
    % Plotting log_phi(CO) against: P1(T1, y1_CO), P2(T2, y1_CO), 
    % P3(T3, y1_CO), ... , P10(T10, y1_CO)
    hold on;
    subplot(2,2,1)
    plot(Pr_CO,values.y1(16).log_phi.CO(n,:),'DisplayName',legend_text_TrCO{n});
   
    if n == 50
    % Lee-Kesler value at Pr = 0.2950; Tr = 2.1391;
    plot(Pr_CO(6), LK.CO.logphi(1),'.r', 'MarkerSize',12,...
        'DisplayName',legend_text_TrCO{P1_length+1})
    % Lee-Kesler value at Pr = 1.4611; Tr = 0.5666;
    plot(Pr_CO(26), LK.CO.logphi(2),'.b', 'MarkerSize',12,...
        'DisplayName',legend_text_TrCO{P1_length+2})
    end
    
    xlim([0 1.6]); xlabel('Pr'); ylabel('log(phi)'); ylim([-3 1]);
    title('Carbon Monoxide Isotherms (y_C_O = 0.3061)');
    legend('show');
    
    % second tile of changing pressure; plots the change in the fugacity
    % coefficient of H2 as a function of pressure at y1_H2 = 0.6939 across
    % 10 different Temperatures. 
    % Plotting log_phi(H2) against: P1(T1, y1_H2), P2(T2, y1_H2), 
    % P3(T3, y1_H2), ... , P10(T10, y1_H2)
    hold on;
    subplot(2,2,2) 
    plot(Pr_H2,values.y1(16).log_phi.H2(n,:),'DisplayName',legend_text_TrH2{n})
   
    if n == 50
    % Lee-Kesler value at Pr = 0.1664; Tr = 3.84;
    plot(Pr_H2(2), LK.H2.logphi(1), '.r', 'MarkerSize', 12,...
        'DisplayName',legend_text_TrH2{P1_length+1})
    % Lee-Kesler value at Pr = 3.9383; Tr = 2.2683;    
    plot(Pr_H2(26), LK.H2.logphi(2), '.b', 'MarkerSize', 12,...
        'DisplayName',legend_text_TrH2{P1_length+2})
    end
   
    xlim([0 4]); xlabel('Pr'); ylabel('log(phi)'); ylim([-2.5 0.5]);
    title('Hydrogen Isotherms (y_H_2 = 0.6939)')
    legend('show');
   
    % third tile of changing pressure; plots the change in the fugacity
    % coefficient of CO as a function of pressure at y2_CO = 0.6939 across
    % 10 different Temperatures. 
    % Plotting log_phi(CO) against: P1(T1, y2_CO), P2(T2, y2_CO), 
    % P3(T3, y2_CO), ... , P10(T10, y2_CO)
    hold on;
    subplot(2,2,3)
    plot(Pr_CO,values.y1(35).log_phi.CO(n,:),'DisplayName',legend_text_TrCO{n})
    
    if n == 50
    % Lee-Kesler value at Pr = 0.2950; Tr = 2.1391;
    plot(Pr_CO(6), LK.CO.logphi(1), '.r', 'MarkerSize', 12,...
        'DisplayName',legend_text_TrCO{P1_length+1})
    % Lee-Kesler value at Pr = 1.4611; Tr = 0.5666;
    plot(Pr_CO(26), LK.CO.logphi(2), '.b', 'MarkerSize', 12,...
        'DisplayName',legend_text_TrCO{P1_length+2})
    % legend('Lee Kessler: Pr = 0.2950; Tr = 2.1391',...
    % 'Lee Kessler: Pr = 1.4611; Tr = 0.5666')
    end
    
    xlim([0 1.5]); xlabel('Pr'); ylabel('log(phi)'); ylim([-3 0.55]);
    title('Carbon Monoxide Isotherms (y_C_O = 0.6939)')
    legend('show')
    
    % fourth tile of changing pressure; plots the change in the fugacity
    % coefficient of H2 as a function of pressure at y2_H2 = 0.3061 across
    % 10 different Temperatures. 
    % Plotting log_phi(H2) against: P1(T1, y2_H2), P2(T2, y2_H2), 
    % P3(T3, y2_H2), ... , P10(T10, y2_H2)
    hold on;
    subplot(2,2,4)
    plot(Pr_H2,values.y1(35).log_phi.H2(n,:),'DisplayName',legend_text_TrH2{n})

    if n == 50
    % Lee-Kesler value at Pr = 0.1664; Tr = 3.84;
    plot(Pr_H2(2), LK.H2.logphi(1), '.r', 'MarkerSize', 12,...
        'DisplayName',legend_text_TrH2{P1_length+1})
    % Lee-Kesler value at Pr = 3.9383; Tr = 2.2683;    
    plot(Pr_H2(26), LK.H2.logphi(2), '.b', 'MarkerSize', 12,...
        'DisplayName',legend_text_TrH2{P1_length+2})
    % legend('Lee Kessler: Pr = 0.1664; Tr = 3.84',...
    % 'Lee Kessler: Pr = 3.9383; Tr = 2.2683')
    end
   
    xlim([0 4]); xlabel('Pr'); ylabel('log(phi)'); ylim([-3 1]);
    title('Hydrogen Isotherms (y_H_2 = 0.3061)')
    legend('show');
    hold on;
end
hold off;


%% Plotting log(phi) vs Reduced Temperature


% Plotting log(phi) vs changing temperatures (constant P and composition)
% using 13 isobars at y1(16) = 0.3061 and y1(34) = 0.6939

figure(2)
hold on;
sgtitle('Changing Temperatures (Constant Pressure and Composition)')
for n = 2:4:T1_length 

if n~=10 % gets rid of an ugly line LOL

    % first tile of changing temperature; plots the change in the fugacity
    % coefficient of CO as a function of temp at y1_CO = 0.3061 across
    % 10 different Pressures. 
    % Plotting log_phi(CO) against: T1(P1, y1_CO), T2(P2, y1_CO), 
    % T3(P3, y1_CO), ... , T10(P10, y1_CO)
    hold on;
    subplot(2,2,1)
    plot(Tr_CO,values.y1(16).log_phi.CO(:,n),'DisplayName',legend_text_PrCO{n})
    
    if n == 50
    % Lee-Kesler value at Pr = 0.2950; Tr = 2.1391;
    plot(Tr_CO(22), LK.CO.logphi(1), '.r', 'MarkerSize', 12,...
        'DisplayName',legend_text_PrCO{P1_length+1})
    % Lee-Kesler value at Pr = 1.4611; Tr = 0.5666;
    plot(Tr_CO(6), LK.CO.logphi(2), '.b', 'MarkerSize', 12,...
        'DisplayName',legend_text_PrCO{P1_length+2})
    % legend('Lee Kessler: Pr = 0.2950; Tr = 2.1391',...
    % 'Lee Kessler: Pr = 1.4611; Tr = 0.5666')
    end
    
    xlim([0.5 4]); xlabel('Tr'); ylabel('log(phi)'); ylim([-2.6 0.7]);
    title('Carbon Monoxide Isobars (y_C_O = 0.3061)')
    legend('show');

    % second tile of changing temperature; plots the change in the fugacity
    % coefficient of H2 as a function of temp at y1_H2 = 0.6939 across
    % 10 different Pressures. 
    % Plotting log_phi(H2) against: T1(P1, y1_H2), T2(P2, y1_H2), 
    % T3(P3, y1_H2), ... , T10(P10, y1_H2)
    hold on;
    subplot(2,2,2)
    plot(Tr_H2,values.y1(16).log_phi.H2(:,n),'DisplayName',legend_text_PrH2{n})
    
    if n == 50
    % Lee-Kesler value at Pr = 0.1664; Tr = 3.84;
    plot(Tr_H2(10), LK.H2.logphi(1), '.r', 'MarkerSize', 12,...
        'DisplayName',legend_text_PrH2{P1_length+1})
    % Lee-Kesler value at Pr = 3.9383; Tr = 2.2683;    
    plot(Tr_H2(6), LK.H2.logphi(2), '.b', 'MarkerSize', 12,...
        'DisplayName',legend_text_PrH2{P1_length+2})
    % legend('Lee Kessler: Pr = 0.1664; Tr = 3.84',...
    % 'Lee Kessler: Pr = 3.9383; Tr = 2.2683')
    end
    
    xlim([0.5 7]); xlabel('Tr'); ylabel('log(phi)'); ylim([-0.08 0.02]);
    title('Hydrogen Isobars (y_H_2 = 0.6939)')
    legend('show');
    
    % third tile of changing temperature; plots the change in the fugacity
    % coefficient of CO as a function of temp at y2_CO = 0.6939 across
    % 10 different Pressures. 
    % Plotting log_phi(CO) against: T1(P1, y2_CO), T2(P2, y2_CO), 
    % T3(P3, y2_CO), ... , T10(P10, y2_CO)
    hold on;
    subplot(2,2,3)
    plot(Tr_CO,values.y1(35).log_phi.CO(:,n),'DisplayName',legend_text_PrCO{n})
    
    if n == 50
    % Lee-Kesler value at Pr = 0.2950; Tr = 2.1391;
    plot(Tr_CO(22), LK.CO.logphi(1), '.r', 'MarkerSize', 12,...
        'DisplayName',legend_text_PrCO{P1_length+1})
    % Lee-Kesler value at Pr = 1.4611; Tr = 0.5666;
    plot(Tr_CO(6), LK.CO.logphi(2), '.b', 'MarkerSize', 12,...
        'DisplayName',legend_text_PrCO{P1_length+2})
    % legend('Lee Kessler: Pr = 0.2950; Tr = 2.1391',...
    % 'Lee Kessler: Pr = 1.4611; Tr = 0.5666')
    end
    
    xlim([0.5 4]); xlabel('Tr'); ylabel('log(phi)'); ylim([-2.6 0.7]);
    title('Carbon Monoxide Isobars (y_C_O = 0.6939)')
    legend('show');

    % fourth tile of changing temperature; plots the change in the fugacity
    % coefficient of H2 as a function of temp at y2_H2 = 0.3061 across
    % 10 different Pressures. 
    % Plotting log_phi(H2) against: T1(P1, y2_H2), T2(P2, y2_H2), 
    % T3(P3, y2_H2), ... , T10(P10, y2_H2)
    hold on;
    subplot(2,2,4)
    plot(Tr_H2,values.y1(35).log_phi.H2(:,n),'DisplayName',legend_text_PrH2{n})
    
    if n == 50
    % Lee-Kesler value at Pr = 0.1664; Tr = 3.84;
    plot(Tr_H2(10), LK.H2.logphi(1), '.r', 'MarkerSize', 12,...
        'DisplayName',legend_text_PrH2{P1_length+1})
    % Lee-Kesler value at Pr = 3.9383; Tr = 2.2683;    
    plot(Tr_H2(6), LK.H2.logphi(2), '.b', 'MarkerSize', 12,...
        'DisplayName',legend_text_PrH2{P1_length+2})
    % legend('Lee Kessler: Pr = 0.1664; Tr = 3.84',...
    % 'Lee Kessler: Pr = 3.9383; Tr = 2.2683')
    end
   
    xlim([1 10]); xlabel('Tr'); ylabel('log(phi)'); ylim([-0.5 0.1]);
    title('Hydrogen Isobars (y_H_2 = 0.3061)')
    legend('show');
end
end
hold off;


%% Values of log(phi) Sorted by Composition


% Making separate variable for log_phi at different compositions and the
% above specified T & P combinations for ease of graphing
log_phi_composition_P12_CO = zeros(y1_length, P1_length);
log_phi_composition_P26_CO = zeros(y1_length, P1_length);
log_phi_composition_P12_H2 = zeros(y1_length, P1_length);
log_phi_composition_P26_H2 = zeros(y1_length, P1_length);

for j = 1:T1_length
    for i = 1:y1_length
        % creates an array where each value is log(phi): the 
        % row corresponds to one specific isobaric, isothermic value (T(j),
        % P(12); and T(j), P(26)). The column corresponds to different
        % composition (y1) values at these consistent P, T combinations.
    log_phi_composition_P12_T_CO(j,i) = values.y1(i).log_phi.CO(12,j);
    log_phi_composition_P26_T_CO(j,i) = values.y1(i).log_phi.CO(26,j);

    log_phi_composition_P12_T_H2(j,i) = values.y1(i).log_phi.H2(12,j);
    log_phi_composition_P26_T_H2(j,i) = values.y1(i).log_phi.H2(26,j);
    end
end

%% Plotting log(phi) vs Composition


% Plotting log(phi) vs changing compositions (constant T and P)
% using 13 isobaric, isothermic lines at: P(12) = 2.05e+06 Pa; 
% and P(26) = 5.11e+06 Pa. This is spanned over all 50 temperatures.
figure(3)
hold on;
sgtitle('Changing Composition (Constant Pressure and Temperature)')

for n = 2:4:P1_length
    
    % first tile of changing composition; plots the change in the fugacity
    % coefficient of CO as a function of comp at P1 = 2.05e+06 Pa across
    % 10 different Temperatures. 
    % Plotting log_phi(CO) against: y1_CO(T1, P1), y2_CO(T2, P1), 
    % y3_CO(T3, P1), ... , y10_CO(T10, P1)
    hold on;
    subplot(2,2,1)
    plot(y1_glob,log_phi_composition_P12_T_CO(n,:),'DisplayName', ...
    strcat(legend_text_TrCO{n}))
    xlim([0 1]); xlabel('y1'); ylabel('log(phi)'); ylim([-0.6 0.25]);
    title('Carbon Monoxide (Isobaric and Isothermic)')
    legend('show');

    % second tile of changing composition; plots the change in the fugacity
    % coefficient of H2 as a function of comp at P1 = 2.05e+06 Pa across
    % 10 different Temperatures. 
    % Plotting log_phi(H2) against: y1_H2(T1, P1), y2_H2(T2, P1), 
    % y3_H2(T3, P1), ... , y10_H2(T10, P1)
    hold on;
    subplot(2,2,2)
    plot(y2_glob,log_phi_composition_P26_T_CO(n,:),'DisplayName', ...
    strcat(legend_text_TrH2{n}))
    xlim([0 1]); xlabel('y2'); ylabel('log(phi)'); ylim([-1 0.3]);
    title('Hydrogen (Isobaric and Isothermic)')
    legend('show');
    
    % third tile of changing composition; plots the change in the fugacity
    % coefficient of CO as a function of comp at P2 = 5.11e+06 Pa across
    % 10 different Temperatures. 
    % Plotting log_phi(CO) against: y1_CO(T1, P2), y2_CO(T2, P2), 
    % y3_CO(T3, P2), ... , y10_CO(T10, P2)
    hold on;
    subplot(2,2,3)
    plot(y1_glob,log_phi_composition_P12_T_H2(n,:),'DisplayName', ...
    strcat(legend_text_TrCO{n}))
    xlim([0 1]); xlabel('y1'); ylabel('log(phi)'); ylim([-0.6 0.1]);
    title('Carbon Monoxide (Isobaric and Isothermic)')
    legend('show'); hold on;
    
    % fourth tile of changing composition; plots the change in the fugacity
    % coefficient of H2 as a function of comp at P2 = 5.11e+06 Pa across
    % 10 different Temperatures. 
    % Plotting log_phi(H2) against: y1_H2(T1, P2), y2_H2(T2, P2), 
    % y3_H2(T3, P2), ... , y10_H2(T10, P2)
    hold on;
    subplot(2,2,4)
    plot(y2_glob,log_phi_composition_P26_T_H2(n,:),'DisplayName', ...
    strcat(legend_text_TrH2{n}))
    xlim([0 1]); xlabel('y2'); ylabel('log(phi)'); ylim([-0.6 0.25])
    title('Hydrogen (Isobaric and Isothermic)')
    legend('show'); hold on;

end
hold off;


%% Function for Fugacity, Fugacity Coefficient (phi), ln(phi), and log(phi)


function [out] = EoS (l,m,q)
% function solving for ln phi, log phi, phi, and fugacity of a binary mix

% l is the input for the number of Pressure values
% m is the input for the number of Temperature values
% q is the input for the number of composition values

% making variables global so they can be accessed from within the function
global P1 T1 structure amix bmix R a1 a2 a12 b1 b2 y1 y2

for z = 1:q

    for k = 1:l
       
        for s = 1:m
    
    % fugacity coefficient for CO
    out.y1(z).ln_phi.CO(k,s) = ((log((R*T1(s))/((structure.vol_CO(k,s)-structure.bmix(z))*(P1(k)))))+ ...
        (structure.bmix(z)/((structure.vol_CO(k,s)-structure.bmix(z)))) ...
        -((2*structure.amix(z))/(R*T1(s)*structure.vol_CO(k,s))));
    out.y1(z).phi.CO(k,s) = exp(out.y1(z).ln_phi.CO(k,s));
    out.y1(z).log_phi.CO(k,s) = log10(out.y1(z).phi.CO(k,s));
    out.y1(z).fugacity.CO(k,s) = out.y1(z).phi.CO(k,s)*(P1(k));
    
    % fugacity coefficient for H2
    out.y1(z).ln_phi.H2(k,s) = ((log((R*T1(s))/((structure.vol_H2(k,s)-structure.bmix(z))*(P1(k)))))+ ...
        (structure.bmix(z)/(structure.vol_H2(k,s)-structure.bmix(z))) ...
        -((2*structure.amix(z))/(R*T1(s)*structure.vol_H2(k,s))));
    out.y1(z).phi.H2(k,s) = exp(out.y1(z).ln_phi.H2(k,s));
    out.y1(z).log_phi.H2(k,s) = log10(out.y1(z).phi.H2(k,s));
    out.y1(z).fugacity.H2(k,s) = out.y1(z).phi.H2(k,s)*(P1(k));

        end

    end

end

end


%% General Interpolation Formula


function y = interp(yl, yh, xl, x, xh)

global LK

% Interpolation formula: y = y_low + (y_high - y_low)*((x-x_low)/(x_high-x_low))

y = yl + (((yh - yl))*((x - xl)/(xh - xl)));

end


%% Solves for log(phi) of a Species Based on Lee-Kesler Equation of State

function [log_phi_CO, log_phi_H2] = Lee_Kesler(y1, y2)

global LK

% takes interpolated data from inputs (solved for earlier in the script)
% to give log_phi at (CO_ideal), (CO_second), (H2_ideal), and (H2_second)

w_CO = 0.049; % acentric factor of CO
w_H2 = -0.22; % acentric factor of H2

% This will give the log(phi)^0 and log(phi)^1 of (Pr, Tr)
% interpolate from xl = P_low, x = Pr, xh = P_high
log_phi0_CO = [interp(y1(1,1),y2(1,1),LK.CO.Pr(1,1),LK.CO.Pr(1,2), ...
               LK.CO.Pr(1,3)); % ideal
               interp(y1(1,2),y2(1,2),LK.CO.Pr(2,1),LK.CO.Pr(2,2), ...
               LK.CO.Pr(2,3))]; % second

log_phi1_CO = [interp(y1(2,1),y2(2,1),LK.CO.Pr(1,1),LK.CO.Pr(1,2), ...
               LK.CO.Pr(1,3)); % ideal
               interp(y1(2,2),y2(2,2),LK.CO.Pr(2,1),LK.CO.Pr(2,2), ...
               LK.CO.Pr(2,3))]; % second

log_phi0_H2 = [interp(y1(1,3),y2(1,3),LK.H2.Pr(1,1),LK.H2.Pr(1,2), ...
               LK.H2.Pr(1,3)); % ideal
               interp(y1(1,4),y2(1,4),LK.H2.Pr(2,1),LK.H2.Pr(2,2), ...
               LK.H2.Pr(2,3))]; % second

log_phi1_H2 = [interp(y1(2,3),y2(2,3),LK.H2.Pr(1,1),LK.H2.Pr(1,2), ...
               LK.H2.Pr(1,3)); % ideal
               interp(y1(2,4),y2(2,4),LK.H2.Pr(2,1),LK.H2.Pr(2,2), ...
               LK.H2.Pr(2,3))]; % second

% Finally, we can solve the departure function with log(phi)^0, log(phi)^1,
% and acentric factor

% returns the log(phi) of CO and H2 at the ideal and second points
log_phi_CO = [log_phi0_CO(1,1) + w_CO*log_phi1_CO(1,1); % ideal
              log_phi0_CO(2,1) + w_CO*log_phi1_CO(2,1)]; % second
log_phi_H2 = [log_phi0_H2(1,1) + w_H2*log_phi1_H2(1,1); % ideal
              log_phi0_H2(2,1) + w_H2*log_phi1_H2(2,1)]; % second

end
