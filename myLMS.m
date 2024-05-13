function [e,w,yk,w_n] = myLMS(M, mu, s, d, w)
    %OUTPUT: [e,w]
    %M                  % szuro egyutthatok szama
    %mu                 % batorsagi tenyezo ~(1/M) (beallitas tapasztalati uton)
    %s                  % reprezentalni kivant rendszer bemenete
    %d                  % reprezentalni kivant rendszer kimenete
    
    
    L = length(s);      % gerj hossza
    %w = zeros(M,1);    % szoro egyutthatok
    x = zeros(M,1);     % forgo vektor, kesleltetok erteket tartalmazza
    e = zeros(L,1);     % hiba vektor
    % y                 % szuro altal becsult kimenet
    yk = zeros(L,1);
    
    l = 1;
    w_n = zeros(M,3);

    for k = 1:L     % k: futo idovaltozo
        x(1)   = s(k);                  % FIR tartalmazza
        y      = w' * x;                % FIR tartalmazza
        e(k)   = d(k)  - y;
        w      = w + 2*mu*e(k)*x;      % modositott coef
        if (k == 100 || k == 200 || k == 500)
            w_n(:,l) = w;
            l =l+1;
        end
        x(2:M) = x(1:(M-1));            % FIR tartalmazza
        yk(k)  = y;
    end
end