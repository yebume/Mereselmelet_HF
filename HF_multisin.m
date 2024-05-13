%% Clean up
    clear;
    close all;
    clc;

%% Parametetek
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

%% Parameterek koherense telete
    t_0 = t_0/1000;

%% 6. Multiszinusz generator
% Parameterek
    M = 100;
    N_minta = 2000;
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
    u = real(y);    % a komplex kinjugaltak kiejtik egymast, csak a szamitasi hiba miatt nem teljesen 0 az im resz

    figure(1);
        plot(u, LineWidth=1.5);
        grid on;
        title('Generált multisinus jel');
        xlabel('Idő [minta]');
        ylabel('Amlitudó [1]');

%% 7. Multiszinusz analizator
[amplitudes, phases, y_error, checksum] = multisinusoid_analyzer(u, M);

% Konvergencia ellenorzese
    figure(2);
        plot(real(y_error), LineWidth=1.5);
        grid on;
        title('A hibajel');
        xlabel('Idő [minta]');
        ylabel('Amlitudó [1]');
    figure(3)
        plot(ones(1,M)-amplitudes);
         grid on;
        title('eredeti - rekreált amplitúdók');
        xlabel('f [1/minta]');
        ylabel('Amlitudó [1]');
    figure(4)
        plot(phi-phases);
        grid on;
        title('eredeti - rekreált fázisok');
        xlabel('f [1/minta]');
        ylabel('Fázis [rad]');


%% Nyolcadik feladat feladat
    % Modellezendo rendszer atviteli fuggvenye
        numerator = conv([0, 1-r], [1, 0 -1]);
        denominator = 2*[1, 0, 0, 0, -r];
        sys = filter(numerator, denominator, [1, zeros(1, 2*M-1)]);

        answear_sys = filter(numerator, denominator, u);
        [sys_amplitudes, sys_phases, sys_y_error, sys_checksum] = multisinusoid_analyzer(answear_sys, M);
    
    figure(5);
        plot(real(sys_y_error), LineWidth=1.5);
        grid on;
        title('A hibajel');
        xlabel('Idő [minta]');
        ylabel('Amlitudó [1]');

    figure(6);
        hold on;
        original_sys = abs(fft(sys));
        plot(20*log10(original_sys(1:M)), LineWidth=1.5);
        plot(20*log10(sys_amplitudes./[0;amplitudes]), LineWidth=1);
        title('Amplitúdó karakterisztika');
        xlabel('f [1/minta]');
        ylabel('Amlitudó [dB]');
        legend("Eredeti 'E' rendszer", "Amplitúdokból számolt");
        grid on;
        
    % Az abrazolashoz szukseges, hogy a fuggveny ertek ne 'ugraljon'
        d_phases = sys_phases-[0;phases];
        for k=1:M
            if (d_phases(k) >= pi)
                d_phases(k) = d_phases(k) - 2*pi;
            end
            if (d_phases(k) < -pi)
                d_phases(k) = d_phases(k) + 2*pi;
            end
        end
    figure(7);
        hold on;
        original_sys = angle(fft(sys));
        plot(original_sys(1:M), LineWidth=1.5);
        plot(d_phases, LineWidth=1);
        title('Fázis karakterisztika');
        xlabel('f [1/minta]');
        ylabel('Fázis [rad]');
        legend("Eredeti 'E' rendszer", "Fázisokból számolt");
        grid on;

