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

%% Elso feladat - Jelgenerator
    % Parameterek
    N_5   =   5;
    N_10  =  10;
    N_100 = 100;        % minta darabszama
   
    Periods = .1;    % vett periodusok
    
    % Valtozok
    t_5   = linspace(0+t_0, ((Periods/f_0)+t_0), N_5);
    t_10  = linspace(0+t_0, ((Periods/f_0)+t_0), N_10);
    t_100 = linspace(0+t_0, ((Periods/f_0)+t_0), N_100);

    % A, B, C parameterek meghatarozasa
    C_aa = o_a^2 .* [1, ro, ro^2; ro, 1, ro; ro^2, ro, 1];
    rng('default');  % Visszaallitjuk a generator alapertelmezett allapotat a reprodukalhatosag erdekeben
    X = mvnrnd([mu_A, mu_B, mu_C], C_aa, 1);
    A = X(1);
    B = X(2);
    C = X(3);

    % Zajos jel generalasa
    sn = @(t) A*sin(2*pi*f_0*t) + B*cos(2*pi*f_0*t) + C + o_w*randn;
    u_5   = zeros(N_5,1);
    u_10  = zeros(N_10,1);
    u_100 = zeros(N_100,1);

    for k=1:N_5
        u_5(k)   = sn(t_5(k));
    end
    for k=1:N_10
        u_10(k)  = sn(t_10(k));
    end
    for k=1:N_100
        u_100(k) = sn(t_100(k));
    end
    figure(1);
        plot(t_5,u_5,t_10,u_10,t_100,u_100, LineWidth=1.5);
        grid on;
        title('Generált zajos jel');
        xlabel('Idő [s]');
        ylabel('Amlitudó [V]');
        legend('5 minta', '10 minta', '100 minta');

%% Masodik feladat - MS becslo
    % Megfigyelesi matrix
    U     = @(n,t) [sin(2*pi*f_0*t'),cos(2*pi*f_0*t'),ones(n,1)];
    C_nn  = @(n) eye(n)*o_w^2;
    C_aaz = @(n,t) U(n,t).'*inv(C_nn(n))*U(n,t)+inv(C_aa);

    MS = @(u,n,t) [mu_A; mu_B; mu_C] + inv(C_aaz(n,t)) * U(n,t).'*inv(C_nn(n))*(u-(U(n,t)*[mu_A;mu_B;mu_C]));

    % a - a tenyleges es a becsult erteket es az elteresuket
    MS_5   = MS(u_5,   N_5,   t_5);
    MS_10  = MS(u_10,  N_10,  t_10);
    MS_100 = MS(u_100, N_100, t_100);
    
    % b - a becsles torzitasat
    b = @(n, a) ([mu_A; mu_B; mu_C] - a)/(1+n*(o_a^2/o_w^2));
    b_5   = b(5, MS_5);
    b_10  = b(10, MS_10);
    b_100 = b(100, MS_100);

    % c - a becsult ertekek kovarianciajat (kovariancia matrixat)
    C_aaz_5   = C_aaz(N_5,   t_5);
    C_aaz_10  = C_aaz(N_10,  t_10);
    C_aaz_100 = C_aaz(N_100, t_100);
   
%% Harmadik feladat - ML becslo (Gauss - Markov)
    GM = @(u, n, t) inv(U(n,t)'*inv(C_nn(n))*U(n,t))*U(n,t)'*inv(C_nn(n))*u;
    
    ML_5   = GM(u_5,   N_5,   t_5);
    ML_10  = GM(u_10,  N_10,  t_10);
    ML_100 = GM(u_100, N_100, t_100);
    
%% Negyedik feladat - LS becslo + minosegjelzo
    LS = @(u, n, t) inv(U(n,t).'*U(n,t))*U(n,t).'*u;
    J = @(u, a_k, n, t) u'*(u- U(n,t)*a_k);

    LS_5   = LS(u_5,   N_5,   t_5);
    LS_10  = LS(u_10,  N_10,  t_10);
    LS_100 = LS(u_100, N_100, t_100);

    J_5    = J(u_5,   LS_5,   N_5,   t_5);
    J_10   = J(u_10,  LS_10,  N_10,  t_10);
    J_100  = J(u_100, LS_100, N_100, t_100);
    
    bar([J_5 J_10, J_100]);

%% Otodik feladat - u(t) = D sin(2pi * f_0 * t + fi) + C
    % Addicios telellel: D*cos(fi)*sin(2*pi*f_0*t) + D*sin(fi)*cos(2*pi*f_0*t) + C
    %                    A*sin(2*pi*f_0*t) + B*cos(2*pi*f_0*t) + C
    % D = sqrt (A^2 + B^2)
    % phi = atan(B/A)
    U_d  = @(n,t) [sin(2*pi*f_0*t'), cos(2*pi*f_0*t'), ones(n,1)];
    LS_d = @(u,n,t) inv(U_d(n,t).'*U_d(n,t))*U_d(n,t).'*u;

    LS_d_5   = LS_d(u_5,   N_5,   t_5);
    LS_d_10  = LS_d(u_10,  N_10,  t_10);
    LS_d_100 = LS_d(u_100, N_100, t_100);

    % Egyenletrendszer megoldasai
    p = LS_d_5;
    D_5  = sqrt(p(1)^2 + p(2)^2);
    fi_5 = atan(p(2)./p(1));
    dC_5 = p(3);

    p = LS_d_10;
    D_10  = sqrt(p(1)^2 + p(2)^2);
    fi_10 = atan(p(2)./p(1));
    dC_10 = p(3);

    p = LS_d_100;
    D_100  = sqrt(p(1)^2 + p(2)^2);
    fi_100 = atan(p(2)./p(1));
    dC_100 = p(3);

    % Minosegjelzo
    J_d = @(u, a_k, n, t) u'*(u- (U_d(n,t)*a_k));
    J_d_5    = J_d(u_5,   [D_5;  fi_5;   dC_5],   N_5,   t_5);
    J_d_10   = J_d(u_10,  [D_10; fi_10;  dC_10],  N_10,  t_10);
    J_d_100  = J_d(u_100, [D_10; fi_100; dC_100], N_100, t_100);

    % CRLB
    % Derivaltak meghatarozasa
        % u = D sin(2pi * f_0 * t + fi) + C
        udD  = @(t,D,fi,C) 1*sin(2*pi*f_0*t+fi); 
        udfi = @(t,D,fi,C) D*cos(2*pi*f_0*t+fi) * 1;
        udC  = @(t,D,fi,C) 1;

    % Hasznalhato ez a keplet ,mert:
    % u = s(a) + w
    % zaj: gauss, mu_zaj = 0, szorasa o_w
        CRLB = @(udx) o_w^2 / ((udx*ones(length(udx),1))^2);  % ahol udx egy sorvektor
    
    % Derivalt jelek meghatarozasa
        % 5 elemu
            s_udD_5  = zeros(1,N_5);
            s_udfi_5 = zeros(1,N_5);
            s_udC_5  = zeros(1,N_5);
            for k = 1:N_5
                s_udD_5(1,k)  = udD( t_5(k), D_5, fi_5, dC_5);
                s_udfi_5(1,k) = udfi(t_5(k), D_5, fi_5, dC_5);
                s_udC_5(1,k)  = udC( t_5(k), D_5, fi_5, dC_5);
            end
        % 10 elemu
            s_udD_10  = zeros(1,N_10);
            s_udfi_10 = zeros(1,N_10);
            s_udC_10  = zeros(1,N_10);
            for k = 1:N_10
                s_udD_10(1,k) =  udD( t_10(k), D_10, fi_10, dC_10);
                s_udfi_10(1,k) = udfi(t_10(k), D_10, fi_10, dC_10);
                s_udC_10(1,k) =  udC( t_10(k), D_10, fi_10, dC_10);
            end
        % 100 elemu
        s_udD_100 = zeros(1,N_100);
        s_udfi_100 = zeros(1,N_100);
        s_udC_100 = zeros(1,N_100);
        for k = 1:N_100
            s_udD_100(1,k) =  udD( t_100(k), D_100, fi_100, dC_100);
            s_udfi_100(1,k) = udfi(t_100(k), D_100, fi_100, dC_100);
            s_udC_100(1,k) =  udC( t_100(k), D_100, fi_100, dC_100);
        end

    % CRLB-k kiszamolasa
        CRLB_D_5   =  CRLB(s_udD_5);
        CRLB_D_10  =  CRLB(s_udD_10);
        CRLB_D_100 =  CRLB(s_udD_100);

        CRLB_fi_5   = CRLB(s_udfi_5);
        CRLB_fi_10  = CRLB(s_udfi_10);
        CRLB_fi_100 = CRLB(s_udfi_100);

        CRLB_C_5   =  CRLB(s_udC_5);
        CRLB_C_10  =  CRLB(s_udC_10);
        CRLB_C_100 =  CRLB(s_udC_100);
        
    %% Kiértékelés
    %{
    k = 1;
    figure(2);
        og = [A, B, C];
        mu = [mu_A,mu_B,mu_C];
        diag_title = ['A paraméter becslői';'B paraméter becslői';'C paraméter becslői'];
        show = [og(k), mu(k), MS_5(k), ML_5(k), LS_5(k); og(k), mu(k),MS_10(k), ML_10(k), LS_10(k); og(k), mu(k), MS_100(k), ML_100(k), LS_100(k); 0, 0, 0, 0, 0];
        bar(show);
        title(diag_title(k,:));
        xticklabels({'5 minta', '10 minta', '100 minta'});
        legend('Eredeti','mu','MS becslő','ML becslő','LS becslő');
        ylabel('Amplitudó [V]');
        grid on;
    %}
    k = 3;
    figure(3);
        if (k == 1)
        show = [A, mu_A, MS_5(1), MS_10(1), MS_100(1); B, mu_B, MS_5(2), MS_10(2), MS_100(2); C, mu_C, ML_5(3), ML_10(3), ML_100(3); 0, 0, 0, 0, 0];
        end
        if (k == 2)
        show = [A, mu_A, ML_5(1), ML_10(1), ML_100(1); B, mu_B, ML_5(2), ML_10(2), ML_100(2); C, mu_C, MS_5(3), MS_10(3), MS_100(3); 0, 0, 0, 0, 0];
        end
        if (k == 3)
        show = [A, mu_A, LS_5(1), LS_10(1), LS_100(1); B, mu_B, LS_5(2), LS_10(2), LS_100(2); C, mu_C, LS_5(3), LS_10(3), LS_100(3); 0, 0, 0, 0, 0];
        end
        bar(show);
        xticklabels({'A', 'B', 'C'});
        ylabel('Amplitudó [V]');
        grid on;
        diag_title = ['MS becslők'; 'ML becslők'; 'LS becslők'];
        title(diag_title(k,:));
        legend('Eredeti', 'mu', '5 minta', '10 minta', '100 minta');

    %{
    k = 1;
    figure(4);
        og = [A, B, C];
        mu = [mu_A,mu_B,mu_C];
        diag_title = ['eredeti érték - A paraméter becslői';'eredeti érték - B paraméter becslői';'eredeti érték - C paraméter becslői'];
        show = [og(k)-MS_5(k), og(k)-ML_5(k), og(k)-LS_5(k); og(k)-MS_10(k), og(k)-ML_10(k), og(k)-LS_10(k); og(k)-MS_100(k), og(k)-ML_100(k), og(k)-LS_100(k); 0, 0, 0];
        bar(show);
        title(diag_title(k,:));
        xticklabels({'5 minta', '10 minta', '100 minta'});
        legend('MS becslő','ML becslő','LS becslő');
        ylabel('Amplitudó [V]');
        grid on;
    %}
    signal = @(t,p) p(1)*sin(2*pi*f_0*t) + p(2)*cos(2*pi*f_0*t) + p(3);

    s_MS_5   = signal(t_5,  MS_5);
    s_MS_10  = signal(t_10, MS_10);
    s_MS_100 = signal(t_100,MS_100);

    s_ML_5   = signal(t_5,  ML_5);
    s_ML_10  = signal(t_10,  ML_10);
    s_ML_100 = signal(t_100,  ML_100);

    s_LS_5   = signal(t_5,  LS_5);
    s_LS_10  = signal(t_10, LS_10);
    s_LS_100 = signal(t_100,LS_100);

    figure(5);
        plot(t_100, u_100 ,t_5, s_MS_5, t_10, s_MS_10, t_100, s_MS_100, LineWidth=1.5);
        title('MS becslők');
        legend('generált zajos jel', 'becslő 5 mintával','becslő 10 mintával','becslő 100 mintával');
        ylabel('Amplitudó [V]');
        xlabel('Idő [s]');
        grid on;
        
    figure(6);
        plot(t_100, u_100 ,t_5, s_ML_5, t_10, s_ML_10, t_100, s_ML_100, LineWidth=1.5);
        title('ML becslők');
        legend('generált zajos jel', 'becslő 5 mintával','becslő 10 mintával','becslő 100 mintával');
        ylabel('Amplitudó [V]');
        xlabel('Idő [s]');
        grid on;

    figure(7);
        plot(t_100, u_100 ,t_5, s_LS_5, t_10, s_LS_10, t_100, s_LS_100, LineWidth=1.5);
        title('LS becslők');
        legend('generált zajos jel', 'becslő 5 mintával','becslő 10 mintával','becslő 100 mintával');
        ylabel('Amplitudó [V]');
        xlabel('Idő [s]');
        grid on;

    
    D_signal = @(t,D,fi,C) D*sin(2*pi*f_0 *t+fi) + C;
    D_LS_5   = D_signal(t_5,   D_5,   fi_5,   dC_5);
    D_LS_10  = D_signal(t_10,  D_10,  fi_10,  dC_10);
    D_LS_100 = D_signal(t_100, D_100, fi_100, dC_100);
    figure(8);
        plot(t_100, u_100 ,t_5, D_LS_5, t_10, D_LS_10, t_100, D_LS_100, LineWidth=1.5);
        title('LS becslők');
        legend('generált zajos jel', 'becslő 5 mintával','becslő 10 mintával','becslő 100 mintával');
        ylabel('Amplitudó [V]');
        xlabel('Idő [s]');
        grid on;
    %{
    figure(9);
        plot(t_100, u_100 ,t_100, s_MS_100, t_100, s_ML_100, t_100, s_LS_100, LineWidth=1);
        title('Becslők');
        legend('generált zajos jel', 'MS', 'ML', 'LS');
        ylabel('Amplitudó [V]');
        xlabel('Idő [s]');
        grid on;
    %}
    figure(10);
        show = [D_5, D_10, D_100; fi_5, fi_10, fi_100; dC_5, dC_10, dC_100];
        bar(show);
        title('D, fi és C LS becslői');
        xticklabels({'D', 'fi', 'C'});
        legend('5 minta', '10 minta', '100 minta');
        ylabel('Amplitudó [V] & Fázis [rad]');
        grid on;

    close all;
