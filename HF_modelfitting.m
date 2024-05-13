%% Clean up
    clear;
    close all;
    clc;

%% Paramétetek
    mu_A        = 1;    % V
    mu_B        = 0;    % V
    mu_C        = 2;    % V
    o_a         = .1;   % V
    o_w         = .1;   % V
    ro          = .1;   % 
    t_0         = 10;   % ms
    f_0         = 50;   % Hz
    rendszer    = 'E';  % 
    r           = .88;  % 
    P           = 25;   % 

%% Paraméterek koherensé tétele
    f_0 = f_0 / 1000; % kHz

%% Multisin
    % Parameterek
        M = 100;
        N_minta = 1000;
        rng('default');  % Visszaallitjuk a generator alapertelmezett allapotat a reprodukalhatosag erdekeben
    % Kezdofazis: 0 vagy veletlenszeru (0 es 2pi kozott)
        phi = 2*pi*rand(M, 1);
        phi = phi - pi;
    
    % Sulyozo egyutthatok generalasa
        x = 1/2 * exp(1i*phi);
    % DC tag, x es konjugaltja
        x = [1 ; x ; flip(conj(x))];
    % Mintaszam
        N = 2*M + 1;
        n = 0 : 1 : N_minta;
        m = 0 : 1 : N-1;
    % Bazis generalasa
        c = exp(1i * 2*pi/N * m' * n);
    % Kimenet generalasa
        y = c.'*x;
        Signal = real(y);
%% System
    numerator = conv([0, 1-r], [1, 0 -1]);
    denominator = 2*[1, 0, 0, 0, -r];
    answear = filter(numerator,denominator,Signal);

    r_2= r/2;
    numerator_2 = conv([0, 1-r_2], [1, 0 -1]);
    denominator_2 = 2*[1, 0, 0, 0, -r_2];
    answear_2 = filter(numerator_2,denominator_2,Signal);

%% 9. feladat
    M = P;
    mu = .001/M; 
    w = zeros(M,1);    % szoro egyutthatok
    [e,w,yk,wn] = myLMS(M, mu, Signal, answear, w);
    [e_2,w_2,yk_2,wn_2] = myLMS(M, mu, Signal, answear_2, w);

    w_m = [wn(:,1),wn(:,2),wn(:,3),w,wn_2(:,1),wn_2(:,2),wn_2(:,3),w_2];

    Imp = [1, zeros(1,M-1)];
    og = filter(numerator, denominator,Imp);
    og_2 = filter(numerator_2, denominator_2,Imp);

    jel_og =20*log10(abs(fft(og)));
    jel_og_2 =20*log10(abs(fft(og_2)));
    jel_1 = 20*log10(abs(fft(w_m(:,1))));
    jel_2 = 20*log10(abs(fft(w_m(:,2))));
    jel_3 = 20*log10(abs(fft(w_m(:,3))));
    jel_4 = 20*log10(abs(fft(w_m(:,4))));
    jel_5 = 20*log10(abs(fft(w_m(:,5))));
    jel_6 = 20*log10(abs(fft(w_m(:,6))));
    jel_7 = 20*log10(abs(fft(w_m(:,7))));
    jel_8 = 20*log10(abs(fft(w_m(:,8))));

    figure(1);
        hold on;
        h=1;
        plot(jel_og(2:ceil(M/h)), LineWidth=3);
        plot(jel_og_2(2:ceil(M/h)), LineWidth=3);
        plot(jel_1(2:ceil(M/h)), LineWidth=1.2);
        plot(jel_2(2:ceil(M/h)), LineWidth=1.2);
        plot(jel_3(2:ceil(M/h)), LineWidth=1.2);
        plot(jel_4(2:ceil(M/h)), LineWidth=1.2);
        plot(jel_5(2:ceil(M/h)), LineWidth=1.2);
        plot(jel_6(2:ceil(M/h)), LineWidth=1.2);
        plot(jel_7(2:ceil(M/h)), LineWidth=1.2);
        plot(jel_8(2:ceil(M/h)), LineWidth=1.2);
        title('Átvitelek');
        xlabel('Minta [1]');
        ylabel('Erősítés [dB]');
        legend('r', 'r/2', '1', '2', '3', '4', '5', '6', '7', '8');
        grid on;
        hold off;

    indicator = zeros(1,2*N_minta);
    indicator(0+100) = 5;
    indicator(0+200) = 5;
    indicator(0+500) = 5;  
    indicator(0+N_minta) = 5;
    indicator(N_minta+100) = 5;
    indicator(N_minta+200) = 5;
    indicator(N_minta+500) = 5;
    indicator(N_minta+N_minta) = 5;
%}
    figure(2);
        hold on
        plot(abs([e;e_2]));
        plot(indicator, LineWidth=2);
        title('LMS hibajele');
        xlabel('Minta [1]');
        ylabel('Amplitúdó [1]');
        grid on;
        hold off;

    figure(3);
        hold on;
        plot(og, LineWidth=3);
        plot(w_m(:,1), LineWidth=1.2);
        plot(w_m(:,2), LineWidth=1.2);
        plot(w_m(:,3), LineWidth=1.2);
        plot(w_m(:,4), LineWidth=1.2);
        title('Konvergencia diagramm (r)');
        xlabel('Minta [1]');
        ylabel('Amplitúdó [1]');
        legend('r', '1', '2', '3', '4');
        grid on;
        hold off;
    figure(4);
        hold on;
        plot(og_2, LineWidth=3);
        plot(w_m(:,5), LineWidth=1.2);
        plot(w_m(:,6), LineWidth=1.2);
        plot(w_m(:,7), LineWidth=1.2);
        plot(w_m(:,8), LineWidth=1.2);
        title('Konvergencia diagramm (r/2)');
        xlabel('Minta [1]');
        ylabel('Amplitúdó [1]');
        legend('r/2', '1', '2', '3', '4');
        grid on;
        hold off;
%}
% Nem mukodik
N_sor = 100000;
u = [1; zeros(N_sor-1,1)];
y = zeros(N_sor,1);

y(1) =  0;
y(2) =  0.06;
y(3) = -0.06;
y(4) =  0;
d = @(n) (1-r)/2 *u(n-1) - (1-r)/2 *u(n-2) + r*y(n-4);

for k=5:N_sor
    y(k) = r/2*y(k-4);
end
%figure(3);
%plot(20*log10(abs(fft(y))));


