close all;

Fe = 48000;         % Fréquence d'échantillonnage
Te = 1/Fe;          % Période d'échantillonnage
Ts = 1/300;         % Période de NRZ
Ns = fix(Ts/Te);    % Nombre d'échantillons
N_bit = 1000;
Nb_bit_secondes = 300;
bit = randi([0,1],1,N_bit);
NRZ = repelem(bit, 1, Ns);

Temps = linspace(0, N_bit/Nb_bit_secondes, N_bit*Ns);

%% Calcul de la DSP de NRZ en échelle fréquentielle.
[DSP, F]=pwelch(NRZ,[],[],[],Fe,'twosided');


%% Calcul de SNRZ
SNRZ = (1/4)*Ts*(sinc(Ts*F).^2) + (1/4)*dirac(F);

%% Construction d'un signal constitué de deux cosinus de fréquence F0 et F1 echantillonnés à Te
F0 = 1180;
F1 = 980;
phi0 = rand*2*pi;
phi1 = rand*2*pi;
x = (1-NRZ) .* cos (2*pi*F0*Temps + phi0) + NRZ .* cos (2*pi*F1*Temps + phi1);

%% Perturbations

Px = mean(abs(x).^2);
SNR = 50;
Sigma = sqrt(Px / 10^(SNR/10));
bruit = Sigma*randn(1,N_bit*Ns);
%bruit = 0;  % A décommenter si on ne veut plus de bruit
x_bruit = x + bruit;

%% 6 Démodulateur de fréquence adapté à la norme V21

% Démodulateur avec gestion du problème de phases 

theta_0 = rand*pi - 50;
theta_1 = rand*5*pi;
x00 = cos(2*pi*F0*reshape(Temps,Ns,N_bit)+theta_0);
x01 = sin(2*pi*F0*reshape(Temps,Ns,N_bit)+theta_0);

x10 = cos(2*pi*F1*reshape(Temps,Ns,N_bit)+theta_1);
x11 = sin(2*pi*F1*reshape(Temps,Ns,N_bit)+theta_1);

inte_0 = sum(reshape(x_bruit,Ns,N_bit).* x00);
inte_1 = sum(reshape(x_bruit,Ns,N_bit).* x01);
inte_2 = sum(reshape(x_bruit,Ns,N_bit).* x10);
inte_3 = sum(reshape(x_bruit,Ns,N_bit).* x11);

inte_0 = inte_0.^2;
inte_1 = inte_1.^2;
inte_2 = inte_2.^2;
inte_3 = inte_3.^2;
terme_moins = inte_0 + inte_1;
terme_plus = inte_2 + inte_3;

Matrice_bits_fin = (terme_plus-terme_moins)>0;
erreur_fin = (Matrice_bits_fin == bit);
taux_erreur_fin = 100-100*mean(erreur_fin)
