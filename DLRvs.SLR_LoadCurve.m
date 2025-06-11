clc; clear;

%% SLR IN DLR (ACRS 110kV) “Drake” 26/7 ACSR 
hours = 1:23;
n = numel(hours);

% Prazni vektorji
I_DLR = zeros(1,n);
I_SLR  = zeros(1,n);
V_arr  = zeros(1,n);
Ta_arr = zeros(1,n);
IT_arr = zeros(1,n);
delta_deg_arr = zeros(1,n);

% Geometrijski in konstantni parametri
D = 30.6e-3;           % premer vodnika [m]
D_1 = 10.19e-3;        % notranji premer jeklenega jedra [m]
d = 3.4e-3;            % zunanji premer vodnika [m]
T_s = 80;              % maksimalna dopustna T vodnika [°C]
eps = 0.75;            % emisijski faktor vodnika
alpha_s = 0.75;        % absorpcijski (sončni) faktor
R20 = 7.283e-5;        % ac upornost Al pri 25°C [Ω/m]
alfa_20 = 0.0045;      % linearni temperaturni koeficient upornosti  
T1 = 20;               % referencna tem.
%delta_deg = 20;       % kot med smerjo vetra in vodnikom [°]
beta = 0;              % nagib vodnika [°]
y = 300;               % nadmorska visina [m]
azim_v_deg= 90;        % azimut vodnika [°] 90=proti vzhodu
phi_deg = 46;          % geogr. širina lokacije [°]
F = 0.1;               % odbito sevanje od tal(albedo)
Z_deg = -15;           % hour angle of the sun
c_f = 1005;            % specificna toplota pri kontanktne tlaku
%ro_20 = 69.0e-8;      % specificna upornost jekla pri 20°C (Operato)
ro_20 = 9.71e-8;       % specificna upornost jekla pri 20°C (splosno)
k_ac = 1.0172;         % faktor koznega pojava
beta_20 = 2.0e-6;      % kvadratni temperaturni koeficient upornosti

sigma = 5.670374e-8;   % Stefan‑Boltzmann [W/m²·K^4]
g = 9.81;              % gravitacijski pospešek [m/s²]
R_air = 287;           % plinska konst. zraka [J/kg·K]
N_i = 117;             % dan v letu

%% 1) IZRACUN DLR
for k = 1:n
    t = hours(k);

    % Ssprememba smeri vetra v stopinjah [°]
    smer_vetra_deg = 130 + 70 * sin(2 * pi * (t - 14) / 24);


    % Izracun beta kota
    delta = abs(smer_vetra_deg - azim_v_deg);
    delta_deg = min(delta, 360 - delta);  % vedno med 0° in 180°


    % Spremenljivke cez dan
    V = 2.5 + 1.0 * sin(2*pi*(t - 14)/24);                                  % hitrost vetra [m/s]
    T_a = 20  + 5.0 * cos(2*pi*(t - 15)/24);                                % temperatura okolice[°C]
    N_s = max(min((t>=5 & t<20) .* (0.8 + 0.2*cos(2*pi*(t-14)/24)),1),0);   % koeficient jasnosti

 
    T_f = 0.5 * (T_s + T_a);                                                % temperatura tanke plasti zraka v neposrednem stiku s povrsino vodnika
    lambda_f = 2.368e-2 + 7.23e-5 * T_f - 2.763e-8 * T_f^2;                 % [W/Km] toplotna prevodnost zraka  pri temperaturi T_f
    gamma = (1.293 - 1.525e-4 * y + 6.379e-9 * y^2) / (1 + 0.00367 * T_f);  % gostoto zraka 
    mu_f = (17.239 + 4.635e-2 * T_f - 2.03e-5 * T_f^2) * 1e-6;              % dinamicno viskoznost zraka 
    nu_f = mu_f / gamma;                                                    % enacba 2.24

    % Soncno obsevanje
    if t>=5 && t<20
        delta_s = deg2rad(23.4 * sin(2*pi*(284+N_i)/365));                  % enacba 2.14
        phi = deg2rad(phi_deg);                                             % zemljepisna sirina
        Z = deg2rad(Z_deg);
        H_s = asin( sin(phi)*sin(delta_s) + cos(phi)*cos(delta_s)*cos(Z) ); % enacba 2.13
        gamma_s = abs(asin( cos(delta_s)*sin(Z) / cos(H_s) ));              % enacba 2.17
        eta = acos( cos(H_s)*cos(gamma_s - deg2rad(azim_v_deg)) );          % enacba 2.16
        I_B0 = N_s * 1280 * sin(H_s) / (sin(H_s) + 0.314);                  % enacba 2.11
        I_B = I_B0 * (1 + 1.4e-4 * y * (1367/I_B0 - 1));                    % enacba 2.12
        I_d = (430.5 - 0.3288 * I_B) * sin(H_s);                            % enacba 2.15
        I_T = I_B * (sin(eta) + (pi/2)*F*sin(H_s)) + I_d*(1+(pi/2)*F);      % globalno sončno sevanje 
    else
        I_T = 0;  % ponoči
    end

    % Prisilna konvekcija(za primer prepletenega vodnika R_s>0.05 / tabela)
    Re = V * D / nu_f;                                                      % Reynoldsovim številom 2.23
    if Re < 2650    
        B = 0.641; n_Re = 0.471;
    else
        B = 0.048; n_Re = 0.8;
    end
    Nu_90 = B * Re^n_Re;                                                    % Nusseltovega števila 2.27
    delta_rad = deg2rad(delta_deg);      
    Nu_delta = Nu_90 * (0.42 + 0.58 * sin(delta_rad)^0.90);                 %(za primer prepletenega vodnika delta>24 / tabela)
    P_c_forced = pi * lambda_f * (T_s - T_a) * Nu_delta; 

    % Naravna konvekcija
    Gr = D^3 * (T_s - T_a) * g / ((T_f + 273) * nu_f^2);    % Grashofovo stevilo 2.32
    Pr = c_f * mu_f / lambda_f;                             % Pradtlovo stevilo 2.33
    Ra  = Gr * Pr;                                          % glede na produkt pogledamo v tabelo in dolocimo koeficienta A in m
    if Ra < 1e2
        A = 1.02; m = 0.148;
    elseif Ra < 1e4
        A = 0.850; m = 0.188;
    elseif Ra < 1e7
        A = 0.480; m = 0.250;
    else      
        A = 0.125; m = 0.333;
    end
    Nu_nat = A * Ra^m;                                      % enacba 2.34
    Nu_beta = Nu_nat * (1 - 6.76e-6 * beta^2.5);            % enacba 2.36
    P_c_nat = pi * lambda_f * (T_s - T_a) * Nu_beta;        % enacba 2.20

    P_c = max(P_c_forced, P_c_nat);                         % vzamemo visjo vrednost

    % Izguba zaradi radiacije
    P_r = pi * D * sigma * eps * ((T_s + 273)^4 - (T_a + 273)^4);  % enacba 2.38

    % Absorpcija soncne energije
    P_s = alpha_s * I_T * D;                                       % enacba 2.10

    % Upornost (razlicne definicije)

    if abs(T_a - T1) < 1e-6                                        % popravek, da preprecimo deljenje z nic
        T_a = T1 + 0.01;
    end
    % Izracun iz cigre - podan R20
    %R_tem     = R20 * (1 + alfa_20*(T_a - T1)); 
    %Rdc       = R20 + (T_s - T1)*(R_tem - R20)/(T_a - T1)
    %Rac = k_ac * Rdc;

    % V splosni dokumentaciji ponavadi podan ro
    ro = ro_20 * ( 1 + alfa_20*(T_a - 20));                         % specifična ohmska upornost 2.5
        %beta_20 * (T_a - 20)^2);                                   % drugi kvadratni del lahko zanemarimo
    A = ((pi * (D)^2) /4);                                          % presek
    Rdc = ro / A;                                                   % enacba 2.4
    Rac = k_ac * Rdc;                                               % enacba 2.7 upostevan kozni efekt
   
    I_DLR(k) = sqrt((P_r + P_c - P_s) / Rac);                       % koncni izracun maksimalnega toka

    V_arr(k) = V;
    Ta_arr(k) = T_a;
    IT_arr(k) = I_T;
    delta_deg_arr(k) = delta_deg;
end

%% 2) IZRACUN SLR ( izracun enak kot DLR samo da vzamem konservativne predpostavke V,T_a in I_T)
for k = 1:n
    % Konstantne spremenljivke
    V = 0.61;                % hitrost vetra
    T_a = 35;                % okoljska temperatura
    I_T = 1000;              % globalno obsevanje
    delta_deg = 60;          % smer vetra glede na vodnik
    
    T_f = 0.5 * (T_s + T_a);
    lambda_f = 2.368e-2 + 7.23e-5 * T_f - 2.763e-8 * T_f^2;
    gamma = (1.293 - 1.525e-4 * y + 6.379e-9 * y^2) / (1 + 0.00367 * T_f);
    mu_f = (17.239 + 4.635e-2 * T_f - 2.03e-5 * T_f^2) * 1e-6;
    nu_f = mu_f / gamma;

    % Prisila in naravna konvekcija (enako kot zgoraj)
    Re = V * D / nu_f;
    if Re < 2650
        B = 0.641; n_Re = 0.471;
    else
        B = 0.048; n_Re = 0.8;
    end
    Nu_90 = B * Re^n_Re;
    Nu_delta = Nu_90 * (0.42 + 0.58 * sin(deg2rad(delta_deg))^0.90);
    P_c_forced= pi * lambda_f * (T_s - T_a) * Nu_delta;

    Gr = D^3 * (T_s - T_a) * g / ((T_f + 273) * nu_f^2);
    Pr = c_f * mu_f / lambda_f;
    Ra = Gr * Pr;
    if Ra < 1e2
        A = 1.02; m = 0.148;
    elseif Ra < 1e4
        A = 0.850; m = 0.188;
    elseif Ra < 1e7
        A = 0.480; m = 0.250;
    else
        A = 0.125; m = 0.333;
    end
    Nu_nat = A * Ra^m;
    Nu_beta = Nu_nat * (1 - 6.76e-6 * beta^2.5);
    P_c_nat = pi * lambda_f * (T_s - T_a) * Nu_beta;
    P_c = max(P_c_forced, P_c_nat);

    % Izguba zaradi sevanja
    P_r = pi * D * sigma * eps * ((T_s + 273)^4 - (T_a + 273)^4);
    P_s = alpha_s * I_T * D;

    % Upornost

    %R_tem = R20 * (1 + alfa_20*(T_a - T1));
    %Rac = R20 + (T_s - T1)*(R_tem - R20)/(T_a - T1);

    ro = ro_20 * ( 1 + alfa_20*(T_a - 20));
    A = ((pi * (D)^2) /4);
    Rdc = ro / A;
    Rac = k_ac * Rdc;

    I_SLR(k) = sqrt((P_r + P_c - P_s) / Rac);
end

%% 3) GRAFI
% Graf: Prikaz treh spremenljivk za izracun DLR
 figure;
 yyaxis left
 plot(hours, V_arr,  'r-', 'LineWidth',1.5); hold on;
 plot(hours, Ta_arr, 'b-', 'LineWidth',1.5);
 ylabel('V (m/s), T_a (°C)');
 
 yyaxis right
 plot(hours, IT_arr, 'g-', 'LineWidth',1.5);
 ylabel('I_T (W/m^2)');
 
 xlabel('Ura [h]');
 title('24-urni potek hitrosti vetra, temperature okolice in soncnega obsevanja');
 legend({'V_{vetra}','T_a','I_T'}, 'Location','best');
 grid on;
 xticks(0:2:24); labels = arrayfun(@(h) sprintf('%02d:00',h), 0:2:24, 'UniformOutput', false);
 xticklabels(labels);
 
 % Prikaz smeri vetra glede na vodnik
 figure;
 plot(hours, delta_deg_arr, 'g-', 'LineWidth', 2);
 xlabel('Ura [h]');
 ylabel('Kot \delta [°]');
 title('24-urni potek smeri vetra glede na vodnik (kot \delta)');
 grid on;
 xlim([0 23]);
 ylim([0 180]);
 xticks(0:2:24);
 labels = arrayfun(@(h) sprintf('%02d:00',h), 0:2:24, 'UniformOutput', false);
 xticklabels(labels);
 
% Graf: potek DLR in SLR cez dan
 figure;
 plot(hours, I_DLR, 'r-', 'LineWidth',1.5); hold on;
 plot(hours, I_SLR, 'b-', 'LineWidth',1.5);
 xlabel('Ura [h]');
 ylabel('Tokovna zmogljivost [A]');
 title('Potek DLR in SLR');
 legend({'I_{DLR}','I_{SLR}'}, 'Location','best');
 grid on;
 xticks(0:2:24);
 xticklabels(labels);
%% 4) Loadability Duration Curve za DLR (24 ur)

P_DLR = zeros(1,n);
P_SLR  = zeros(1,n);

U = 110e3;
P_SLR(k) = sqrt(3) * U * I_SLR(k) / 1e6;
P_DLR(k) = sqrt(3) * U * I_DLR(k) / 1e6;

% DLR tokovna zmogljivost od najvecje proti najmanjsi
I_sorted = sort(I_DLR, 'descend');

% Casa v dnevu (ker imamo 23 vrednosti, vsak predstavlja ~4.35 %)
n = length(I_sorted);
percent_time = ((1:n) / n) * 100;

% Staticni tok
I_static = mean(I_SLR); 

% Izris Loadability Duration Curve
figure;
plot(percent_time, I_sorted, 'g-', 'LineWidth', 2); hold on;
xlabel('Procent časa v dnevu [%]');
ylabel('Tokovna zmogljivost [A]');
title('Loadability Duration Curve – 24h DLR');

% Staticne meje
yline(I_static, '--', 'static limit', 'Color', 'k');
yline(1.1 * I_static, '--', '110% static limit', 'Color', [0 0 1]);
yline(1.2 * I_static, '--', '120% static limit', 'Color', [0.2 0.6 1]);

% Navpicne crte za vizualno interpretacijo
xline(60, '--', 'Color', [0.3 0.3 0.8]);
xline(80, '--', 'Color', [0.3 0.3 0.8]);
xline(90, '--', 'Color', [0.3 0.3 0.8]);

legend('I_{DLR}', 'Static limit', '110% limit', '120% limit', ...
       '60% časa', '80% časa', '90% časa', 'Location','northeast');
grid on;

% Razlika med SLR in DLR
% figure;
% plot(hours, I_DLR - I_SLR, 'k-o', 'LineWidth', 1.5);
% xlabel('Ura [h]');
% ylabel('Razlika I_{DLR} - I_{SLR} [A]');
% title('Razlika med DLR in SLR čez dan');
% grid on;

%% 5) Loadability Duration Curve za moč (MW)

% Napetostni nivo v V
%U = 110e3;

% Pretvorba toka v moč [MW]
%P_DLR = sqrt(3) * U * I_DLR / 1e6;
%P_SLR = sqrt(3) * U * I_SLR / 1e6;

% DLR moc od najvecje proti najmanjsi
%P_sorted = sort(P_DLR, 'descend');
%n = length(P_sorted);
%percent_time = (1:n) / n * 100;

% Statična moc (konzervativno npr. P_SLR(1) ali povprečje)
%P_static = mean(P_SLR);

% Izris grafa
%figure;
%plot(percent_time, P_sorted, 'm-', 'LineWidth', 2); hold on;
%xlabel('Procent časa v dnevu [%]');
%ylabel('Moč [MW]');
%title('Loadability Duration Curve – 24h DLR (moč)');

% Staticne meje
%yline(P_static, '--', 'Static limit', 'Color', 'k');
%yline(1.1 * P_static, '--', '110% static limit', 'Color', [0 0 1]);
%yline(1.2 * P_static, '--', '120% static limit', 'Color', [0.2 0.6 1]);

% Navpicne crte za casovno interpretacijo
%xline(60, '--', 'Color', [0.3 0.3 0.8]);
%xline(80, '--', 'Color', [0.3 0.3 0.8]);
%xline(90, '--', 'Color', [0.3 0.3 0.8]);

%legend('P_{DLR}', 'Static limit', '110% limit', '120% limit', ...
%       '60% časa', '80% časa', '90% časa', 'Location','northeast');
%grid on;
