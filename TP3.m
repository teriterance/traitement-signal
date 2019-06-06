%% TP3 - filtrage anti-repliement par la méthode de la transposition bilinéaire

close all
clear all
clc


%% I Sous-échantillonnage et aliasing
% 1.a Chargement et représentation temporelle
% A COMPLETER
%chargement du son qui sera exploite
[signal, fe] = audioread('glockenspiel.wav');
N = length(signal);
%la duree de son enregistrement est:
duree =  (N-1)/fe;
%representation temporelle 
t = 0:1:N-1;
t = t/fe;
plot(t,signal)
xlabel('temps en(s)');
ylabel('amplitude');
title('Amplitude du signal de depart');
%ecoute du son
soundsc(signal, fe);

% 1.b Calcul de la TFD et spectre en amplitude
% A COMPLETER
%calcul de la fft
Nfft = 3*N;% on definit un grand nombre de valeur
Fft_signal = fftshift(fft(signal, Nfft));
%trace de son module
figure;
f = -fe/2:fe/Nfft:fe/2-fe/Nfft;
plot(f,abs(Fft_signal));
xlabel('frequence en(Hz)');
ylabel('amplitude');
title('Module de la fft du signal de depart');

% 1.c Sous-échantillonnage du signal
% A COMPLETER
%le sous echantillonage
r = 4; % facteur de sous echantillonage
signal_sub = zeros(floor(N/r),2);
for i = 1:r:floor(N/r)
    signal_sub(i,1) = signal(2*i,1);
    signal_sub(i,2) = signal(2*i, 2);
end 
% la duree du signal ne devrai pas etre considere comme modifier car notre
% sous echantillonage s'acompagne d'une nouvelle frequence qui est la
% moitie de la frequence originale soit une nouvelle periode 2 fois plus
% longue
fer = fe/r; %definition de la frequence du signal sous echantillonee(au sous echantillonage)

% 1.d Calcul de la TFD et spectre en amplitude du signal sous-échantillonné
% A COMPLETER
%calcul de la fft
Fft_signal_sub = fftshift(fft(signal_sub, Nfft));
%trace de son module
figure;
f = -fer/2:fer/Nfft:fer/2-fer/Nfft;
plot(f,abs(Fft_signal_sub));
xlabel('frequence en(Hz)');
ylabel('amplitude');
title('module de la fft du signal sous echantillone');
%ecoute du son
soundsc(signal_sub, fe);

%plus on augmente r, plus on a de la peine a distinguer les deux signaux pl, ceci est du au fait qu'en
%bassant la frequence d'echantillonage, on suprime certaine partie du
%signal en temps.

%% II Gabarit
% COMPLETER les paramètres du gabarit ci-dessous
fc = fer/2;
fp = fer/4;
Gmin = 0.95;
Gmax = 0.001;
Gabf = [0 fp fp; fc fc fe/2];
GabH = [Gmin Gmin Gmax; Gmin Gmax Gmax];
figure
plot(Gabf',GabH');
xlabel('frequence (Hz)');
ylabel('Gain');
title('Gabarit du filtre')
xlim([0 fe/2])


%% III Filtre de Butterworth
% 3.a Calcul des paramètres du filtre
% - l'ordre N_b
% - la fréquence de coupure (à -3dB) f0_b
% COMPLETER les paramètres du filtre
N_b  = (log10(10^(Gmax/10)-1) - log10(10^(Gmin/10)-1))/(2*log10(0.5)); 
N_b = ceil(N_b);%passage a la partie entiere
f0_b = fp/( (10^(Gmax/10)-1)^(1/(2*N_b)));

%-------------------------------------PARTIE DONNEE : PAS TOUCHE !
% Création et visualisation de la fonction de transfert
[b_b,a_b] = butter(N_b,f0_b/(fe/2)); 
[H_butter,f_butter] = freqz(b_b,a_b,'whole',N);
figure;
plot(f_butter*(fe/2)/pi,abs(H_butter));
xlabel('frequence (Hz)');
ylabel('gain')
title('Fonction de transfert en f du filtre de Butterworth : gain');
hold on
plot(Gabf',GabH','r');
xlim([0 fe/2])

A = H_butter./abs(H_butter);
phi_butter = [0; cumsum(angle(A(2:end)./A(1:end-1))/pi*180)];
figure;
plot(f_butter*(fe/2)/pi,phi_butter);
xlabel('frequence (Hz)')
ylabel('phase (degres)')
title('Fonction de transfert en f du filtre de Butterworth : phase');
xlim([0 fe/2])
ylim([-1200 0])
%-------------------------------------


% 3.b Filtrage
% A COMPLETER
rep_imp_filtre = ifft(H_butter);
signal_filtre = zeros(N*2-1,2);
signal_filtre(:,1) = conv(rep_imp_filtre,signal(:,1));
signal_filtre(:,2) = conv(rep_imp_filtre,signal(:,2));
% 3.c Sous-échantillonnage
% A COMPLETER
N2 = length(signal_filtre);
r2 = 4; % facteur de sous echantillonage
signal_sub_2 = zeros(floor(N2/r2),2);
for i = 1:r2:floor(N/r2)
    signal_sub_2(i,1) = signal_filtre(2*i,1);
    signal_sub_2(i,2) = signal_filtre(2*i, 2);
end 
% 3.d Calcul de la TFD et spectre en amplitude
% A COMPLETER

Fft_signal_filtre = fftshift(fft(signal_filtre, Nfft));
%trace de son module
figure;
fer2 = fe/r2;
f = -fer2/2:fer2/Nfft:fer2/2-fer2/Nfft;
plot(f,abs(Fft_signal_filtre));
xlabel('frequence en(Hz)');
ylabel('amplitude');
title('module de la fft du signal filtre Butterworth sous echantillone');
%ecoute du son
soundsc(signal_filtre, fer2);

%% IV Filtre de Tchebychev
% 4.a et 4.b Calcul des parmatères du filtre
% - la fréquence de coupure f0_t
% - le taux d'ondulation epsilon2
% - l'ordre du filtre N_t
% COMPLETER les paramètres du filtre
f0_t     = fp;
epsilon2 = 10^(Gmax/10)-1;
N_t      = acos((10^(Gmin/10)-1)*epsilon2^2*f0_t/fc)/acos(f0_t/fc);
N_t = ceil(N_t);
%-------------------------------------PARTIE DONNEE : PAS TOUCHE !
% Création du filtre et visualistion de la fonction de transfert
[b_t,a_t] = cheby1(N_t,10*log10(1/Gmin),f0_t/(fe/2));
[H_cheby,f_cheby] = freqz(b_t,a_t,'whole',N);
figure;
plot(f_cheby*(fe/2)/pi,abs(H_cheby));
xlabel('frequence (Hz)');
ylabel('gain')
title('Fonction de transfert en f du filtre de Tchebychev : gain');
hold on
plot(Gabf',GabH','r');
xlim([0 fe/2])

A   = H_cheby./abs(H_cheby);
phi_cheby = [0; cumsum(angle(A(2:end)./A(1:end-1))/pi*180)];
figure;
plot(f_cheby*(fe/2)/pi,phi_cheby);
xlabel('frequence (Hz)')
ylabel('phase (degres)')
title('Fonction de transfert en f du filtre de Tchebychev : phase');
xlim([0 fe/2])
ylim([-1000 0])
%-------------------------------------


% 4.c Filtrage
% A COMPLETER
rep_imp_filtre2 = ifft(H_cheby);%reponce impulsionnele du filtre de Tchebychev
signal_filtre2 = zeros(N*2-1,2);
signal_filtre2(:,1) = conv(rep_imp_filtre2,signal(:,1));
signal_filtre2(:,2) = conv(rep_imp_filtre2,signal(:,2));

% 4.d Sous-échantillonnage
% A COMPLETER
N2 = length(signal_filtre2);
r2 = 4; % facteur de sous echantillonage
signal_sub_3 = zeros(floor(N2/r2),2);
for i = 1:r2:floor(N/r2)
    signal_sub_3(i,1) = signal_filtre2(2*i,1);
    signal_sub_3(i,2) = signal_filtre2(2*i, 2);
end 

% 4.e Calcul de la TFD et spectre en amplitude
% A COMPLETER
Fft_signal_filtre2 = fftshift(fft(signal_filtre2, Nfft));
%trace de son module
figure;
fer2 = fe/r2;
f = -fer2/2:fer2/Nfft:fer2/2-fer2/Nfft;
plot(f,abs(Fft_signal_filtre2));
xlabel('frequence en(Hz)');
ylabel('amplitude');
title('module de la Fft du signal filtre selon Tchebychev sous echantillone');
%ecoute du son
soundsc(signal_filtre, fer2);

