function sortie = demoduler(signal)

Fe = 48000;                     % Fréquence d'échantillonnage
Te = 1/Fe;                      % Période d'échantillonnage
Nb_bit_secondes = 300;          % Débits de bits
Ts = 1/Nb_bit_secondes;         % Période de NRZ
Ns = fix(Ts/Te);                % Nombre d'échantillons
N_bit = length(signal)/Ns;      % Nombre de bits

F0 = 1180;
F1 = 980;
theta_0 = rand*2*pi;
theta_1 = rand*2*pi;
t = 0:Te:(length(signal)-1)*Te;
x00 = cos(2*pi*F0*t+theta_0);
x01 = sin(2*pi*F0*t+theta_0);

x10 = cos(2*pi*F1*t+theta_1);
x11 = sin(2*pi*F1*t+theta_1);

inte_0 = reshape(signal.*x00,Ns,N_bit);
inte_1 = reshape(signal.*x01,Ns,N_bit);
inte_2 = reshape(signal.*x10,Ns,N_bit);
inte_3 = reshape(signal.*x11,Ns,N_bit);

inte_0 = sum(inte_0);
inte_1 = sum(inte_1);
inte_2 = sum(inte_2);
inte_3 = sum(inte_3);

inte_0 = inte_0.^2;
inte_1 = inte_1.^2;
inte_2 = inte_2.^2;
inte_3 = inte_3.^2;
terme_moins = inte_0 + inte_1;
terme_plus = inte_2 + inte_3;

sortie = (terme_plus-terme_moins)>0;
end