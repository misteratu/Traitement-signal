close all;

Fe = 48000;         % Fréquence d'échantillonnage
Te = 1/Fe;          % Période d'échantillonnage
Ts = 1/300;         % Période de NRZ
Ns = fix(Ts/Te);    % Nombre d'échantillons
N_bit = 1000;
Nb_bit_secondes = 300;
bit = randi([0,1],1,N_bit);
NRZ = repelem(bit, 1, Ns);

%% Tracé de NRZ en échelle temporelle.
Temps = linspace(0, N_bit/Nb_bit_secondes, N_bit*Ns);

%% Tracé de la DSP de NRZ en échelle fréquentielle.
[DSP, F]=pwelch(NRZ,[],[],[],Fe,'twosided');

%% Tracé de SNRZ
SNRZ = (1/4)*Ts*(sinc(Ts*F).^2) + (1/4)*dirac(F);

%% Tracé des deux cosinus de fréquence F0 et F1 echantillonnés à Te
F0 = 1180;
F1 = 980;
phi0 = rand*2*pi;
phi1 = rand*2*pi;
x = (1-NRZ) .* cos (2*pi*F0*Temps + phi0) + NRZ .* cos (2*pi*F1*Temps + phi1);

% Densité spectrale de puissance du signal modulé en fréquence.
[S_x_estime, F_x] = pwelch (x, [], [], [], Fe, 'twosided');

%% Perturbations

Px = mean(abs(x).^2);
SNR = 50;
Sigma = sqrt(Px / 10^(SNR/10));
bruit = Sigma*randn(1,N_bit*Ns);
%bruit = 0;  % A décommenter si on ne veut plus de bruit
x_bruit = x + bruit;
% Densité spectrale de puissance du signal modulé en fréquence.
[S_x_estime_bruit, F_x_bruit] = pwelch (x_bruit, [], [], [], Fe, 'twosided');


%% 6 Démodulateur de fréquence adapté à la norme V21

% 6.1) Contexte de synchronisation idéale

x0 = cos(2*pi*F0*reshape(Temps,Ns,N_bit)+phi0);
x1 = cos(2*pi*F1*reshape(Temps,Ns,N_bit)+phi1);

inte_0 = sum(reshape(x_bruit,Ns,N_bit).* x0);
inte_1 = sum(reshape(x_bruit,Ns,N_bit).* x1);

Matrice_bits2 = (inte_1-inte_0)>0;
erreur_s = (Matrice_bits2 == bit);
taux_erreur_synchronise = 100-100*mean(erreur_s)





