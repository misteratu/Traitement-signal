close all;

Fe = 48000;         % Fréquence d'échantillonnage
Te = 1/Fe;          % Période d'échantillonnage
Ts = 1/300;         % Période de NRZ
Ns = fix(Ts/Te);    % Nombre d'échantillons
N_bit = 1000;       % Nombre de bits
Nb_bit_secondes = 300;      % Débits de bits
bit = randi([0,1],1,N_bit); % Vecteur de bits de référence à faire passer dans le modem

%% Construction du signal NRZ

NRZ = repelem(bit, 1, Ns);  

%% Tracé de NRZ en échelle temporelle.

Temps = linspace(0, N_bit/Nb_bit_secondes, N_bit*Ns);
figure('Name','Figure 1 : NRZ','NumberTitle','off'); 
p = plot(Temps,NRZ);
p.LineWidth = 3;
xlabel("Temps en secondes (s)");
ylabel("NRZ");
title("Tracé de NRZ en fonction du temps");

%% Tracé de la DSP de NRZ en échelle fréquentielle.
[DSP, F]=pwelch(NRZ,[],[],[],Fe,'twosided');

figure('Name','Figure 2 : Densité spectrale de NRZ','NumberTitle','off');
Freq = linspace(-Fe/2, Fe/2, length(F));
semilogy(Freq,fftshift(DSP));
xlabel("f(Hz)");
ylabel("DSP de NRZ");
title("Tracé de DSP(NRZ) = f(f)");

figure('Name','Figure 3 : Densité spectrale de NRZ et SNRZ','NumberTitle','off'); 
semilogy(Freq,fftshift(DSP));
xlabel("f(Hz)");
ylabel("DSP de NRZ");
title("Tracé de DSP(NRZ) = f(f)");

%% Construction de SNRZ

SNRZ = (1/4)*Ts*(sinc(Ts*Freq).^2) + (1/4)*dirac(Freq);


%% Tracé de SNRZ

hold on
semilogy(Freq,SNRZ);
xlabel("f(Hz)");
ylabel("DSP de NRZ et SNRZ");
title("Tracé de DSP(NRZ) = f(f) et SNRZ");
legend("DSP de NRZ", "SNRZ");

%% Construction d'un signal constitué de deux cosinus de fréquence F0 et F1 echantillonnés à Te
F0 = 1180;  %1180 ou 6000
F1 = 980;   %980 ou 2000
phi0 = rand*2*pi;
phi1 = rand*2*pi;

x = (1-NRZ) .* cos (2*pi*F0*Temps + phi0) + NRZ .* cos (2*pi*F1*Temps + phi1);

% Densité spectrale de puissance du signal modulé en fréquence.
[S_x_estime, F_x] = pwelch (x, [], [], [], Fe, 'twosided');

figure('Name','Figure 4 : Densité spectrale de x','NumberTitle','off');
semilogy(linspace(-Fe/2, Fe/2, length(F_x)), fftshift(S_x_estime));
xlabel("f(Hz)");
ylabel("DSP de x");
title("Tracé de DSP(x) = f(f)");

%% Perturbations et construction d'un signal bruité

Px = mean(abs(x).^2);
SNR = 50;
Sigma = sqrt(Px / 10^(SNR/10));
bruit = Sigma*randn(1,N_bit*Ns);
bruit = 0;  % A décommenter si on ne veut plus de bruit
x_bruit = x + bruit;

% Densité spectrale de puissance du signal modulé en fréquence.
[S_x_estime_bruit, F_x_bruit] = pwelch (x_bruit, [], [], [], Fe, 'twosided');

%% Affichage de la comparaison du signal bruité avec le non bruité 

figure('Name','Figure 5 : Tracé du signal x','NumberTitle','off')
plot(Temps, x);
xlabel("Temps en secondes");
ylabel("x");
title("Tracé de x en fonction du temps");

figure('Name','Figure 6 : Tracé du signal avant et après le bruitage','NumberTitle','off')
subplot(2,1,1);
plot(Temps, x);
xlabel("Temps en secondes");
ylabel("x");
title("Tracé de x non bruité");
subplot(2,1,2);
plot(Temps, x_bruit);
xlabel("Temps en secondes");
ylabel("x bruité");
title("Tracé de x bruité");

%% Démodulation par filtrage (Passe bas)

% Construction du filtre passe bas 

B = (F1 + F0)/2;
taille = 30; % Mettre 30 si on veut tester de nouveau
Taille_filtre = [-taille:1:taille];
Passe_bas_i = 2*B/Fe*sinc(2*B/Fe*Taille_filtre);

% Affichage de la réponse impulsionnelle du filtre passe bas 

figure('Name','Figure 7 : Réponses fréquentielles et impulsionnelles des filtres','NumberTitle','off');
subplot(2,2,1);
plot(Taille_filtre/Fe, Passe_bas_i);
xlabel("Temps en secondes");
ylabel("h(t)");
title("Tracé de la réponse impulsionnelle du passe bas");

% Affichage de la réponse en fréquence du filtre passe bas 

Passe_bas_f = fftshift(abs(fft(Passe_bas_i)));
subplot(2,2,2);
plot(linspace(-Fe/2, Fe/2, length(Passe_bas_f)), Passe_bas_f);
xlabel("Fréquence en Hertz");
ylabel("H(f)");
title("Tracé de la réponse fréquentielle du passe bas");

%% Démodulation par filtrage (Passe haut)

% Construction filtre passe haut

Passe_haut_f = 1 - Passe_bas_f;
Passe_haut_i = - Passe_bas_i;
Passe_haut_i(taille + 1) = 1 + Passe_haut_i(taille + 1);

% Affichage de la réponse impulsionnelle du filtre passe haut

subplot(2,2,3);
plot(Taille_filtre/Fe, real(Passe_haut_i));
xlabel("Temps en secondes");
ylabel("h(t)");
title("Tracé de la réponse impulsionnelle du passe haut");

% Affichage de la réponse en fréquence du filtre passe haut 

subplot(2,2,4);
plot(linspace(-Fe/2, Fe/2, length(Passe_haut_f)), Passe_haut_f);
xlabel("Fréquence en Hertz");
ylabel("H(f)");
title("Tracé de la réponse fréquentielle du passe haut");

%% Sortie des filtres

Y_h = filter(Passe_haut_i, 1, x_bruit); % Signal en sortie du filtre passe haut
Y_b = filter(Passe_bas_i, 1, x_bruit);  % Signal en sortie du filtre passe bas

% Affichage des réponses fréquentielles et densité spectrale de l'entrée et de la sortie des filtres


%% Tracé de la DSP de la sortie du filtre en échelle fréquentielle.

% Densité spectrale de puissance du signal bruité modulé en fréquence.
[S_x_estime_bruit, F_x_bruit] = pwelch (x_bruit, [], [], [], Fe, 'onesided');

% DSP de la sortie du passe bas 
[DSP_b, F_b]=pwelch(Y_b,[],[],[],Fe,'onesided');


figure('Name','Figure 8 : Réponses fréquentielles et densité spectrale de l`entrée et de la sortie des filtres','NumberTitle','off');

% Affichage de la réponse fréquentielle et densité spectrale de l'entrée du
% filtre passe bas

subplot(2,2,1);
plot(linspace(-Fe/2, Fe/2, length(Passe_bas_f)), Passe_bas_f);
hold on
semilogy(F_x_bruit,S_x_estime_bruit/max(S_x_estime_bruit));
xlabel("fréquence en Hz");
ylabel("DSP");
title("Tracé de la DSP de l'entrée du filtre passe bas");
legend("Réponse fréquentielle du passe bas", "Densité spectrale");

% DSP de la sortie du passe haut
[DSP_h, F_h]=pwelch(Y_h,[],[],[],Fe,'onesided');

% Affichage de la réponse fréquentielle et densité spectrale de l'entrée du
% filtre passe haut

subplot(2,2,2);
plot(linspace(-Fe/2, Fe/2, length(Passe_haut_f)), Passe_haut_f);
hold on
semilogy(F_x_bruit,S_x_estime_bruit/max(S_x_estime_bruit));
xlabel("fréquence en Hz");
ylabel("DSP");
title("Tracé de la DSP de l'entrée du filtre passe haut");
legend("Réponse fréquentielle du passe haut", "Densité spectrale");

% Affichage de la réponse fréquentielle et densité spectrale de la sortie du
% filtre passe bas

subplot(2,2,3);
plot(linspace(-Fe/2, Fe/2, length(Passe_bas_f)), Passe_bas_f);
hold on
semilogy(F_b,DSP_b/max(DSP_b));
xlabel("fréquence en Hz");
ylabel("DSP");
title("Tracé de la DSP de la sortie du filtre passe bas");
legend("Réponse fréquentielle du passe bas", "Densité spectrale");

% Affichage de la réponse fréquentielle et densité spectrale de la sortie du
% filtre passe haut

subplot(2,2,4);
plot(linspace(-Fe/2, Fe/2, length(Passe_haut_f)), Passe_haut_f);
hold on
semilogy(F_h,DSP_h/max(DSP_h));
xlabel("fréquence en Hz");
ylabel("DSP");
title("Tracé de la DSP de la sortie du filtre passe haut");
legend("Réponse fréquentielle du passe haut", "Densité spectrale");

%% Détection d'énergie

% Signal filtré passe bas 
signal_demodule_b = reshape(Y_b,Ns,N_bit);

figure('Name','Figure 9 : Sortie des filtres','NumberTitle','off');
subplot(2,1,1);
plot(Y_b)
title("Signal démodulé après filtrage passe bas")
xlabel("Temps (t)")
ylabel("y(t)")

% Signal filtré passe haut
signal_demodule_h = reshape(Y_h,Ns,N_bit);

subplot(2,1,2);
plot(Y_h)
title("Signal démodulé après filtrage passe haut")
xlabel("Temps (t)")
ylabel("y(t)")

%% Taux d'erreur pour les filtres 

% Taux d'erreur pour le passe bas

energie_b = sum(signal_demodule_b.^2);
K = (max(energie_b) + min(energie_b))/2;

Matrice_bit_b = energie_b > K;
erreur = (Matrice_bit_b == bit);
Taux_erreur_b = 100-100*mean(erreur);

% Taux d'erreur pour le passe haut

energie_h = sum(signal_demodule_h.^2);
K = (max(energie_h) + min(energie_h))/2;
Matrice_bit_h = energie_h < K;
erreur = (Matrice_bit_h == bit);
Taux_erreur_h = 100-100*mean(erreur);

%% 6 Démodulateur de fréquence adapté à la norme V21

%% 6.1) Contexte de synchronisation idéale

x0 = cos(2*pi*F0*reshape(Temps,Ns,N_bit)+phi0);
x1 = cos(2*pi*F1*reshape(Temps,Ns,N_bit)+phi1);

inte_0 = sum(reshape(x_bruit,Ns,N_bit).* x0);
inte_1 = sum(reshape(x_bruit,Ns,N_bit).* x1);

Matrice_bits2 = (inte_1-inte_0)>0;
erreur_s = (Matrice_bits2 == bit);
taux_erreur_synchronise = 100-100*mean(erreur_s);

%% 6.2) Présence d'une erreur de synchronisation de phase porteuse
phi0_prime = 3000;
phi1_prime = 7000;
x0 = cos(2*pi*F0*reshape(Temps,Ns,N_bit)+phi0_prime);
x1 = cos(2*pi*F1*reshape(Temps,Ns,N_bit)+phi1_prime);

inte_0 = sum(reshape(x_bruit,Ns,N_bit).* x0);
inte_1 = sum(reshape(x_bruit,Ns,N_bit).* x1);

Matrice_bits2 = (inte_1-inte_0)>0;
erreur_e = (Matrice_bits2 == bit);
taux_erreur_erreur_synchronise = 100-100*mean(erreur_e);

%% Récepteur avec gestion d'erreur de synchronisation des phases 

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
taux_erreur_fin = 100-100*mean(erreur_fin);

