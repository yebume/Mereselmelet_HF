function [amplitudes, phases, y_err, check_sum] = multisinusoid_analyzer(y, M)

    N = (2*M)+1;
    n = 0 : 1 : length(y)-1;
    m = 0 : 1 : N-1;

    amplitudes = zeros(M,1);
    phases = zeros(M,1);

    % Bazis generalasa
        c = exp(1i * 2*pi/N * m' * n);
    % Reciprok bazis generalasa
        g = 1/N * exp(-1i*2*pi/N*m'*n); % 1/N * conj(c)

    % Becslo szamitasa rekurziv modon
        x_k   = zeros(N,   1);
        y_k   = zeros(N+1, 1);
        y_err = zeros(N+1, 1);

        for k = 1 : length(y)
            % Kimenet becsloje
                y_k(k) = c(:, k).' * x_k;
            % Kimenet becslojenek hibaja
                y_err(k) = y(k) - y_k(k);
            % Sulytenyezok becsloje
                x_k = x_k + g(:,k) * y_err(k);
        end
    % N lepeses konvergencia ellenorzese
        check_sum = zeros(N, N);
        for k = 1:N
            check_sum = check_sum + g(:,k) * conj(c(:,k)');
        end
    % Vegertekek beallitasa
        for k = 1:M+1
            amplitudes(k) = 2*abs(x_k(k));
            phases(k)     = angle(x_k(k));
        end
end