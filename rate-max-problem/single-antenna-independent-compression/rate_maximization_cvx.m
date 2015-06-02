function [Omega R objective feasible] = rate_maximization_cvx(params)
    L = params.L;
    K = params.K;
    Nl = params.Nl;
    sigma = params.sigma;
    h = params.h;
    E = params.E;
    r = params.r;
    C = params.C;
    P = params.P;
    omega = params.omega;
    
    cvx_begin
        variable OmegaI(Nl, Nl, L) hermitian semidefinite
        variable R(Nl * L, Nl * L, K) hermitian semidefinite
        
        Omega = blkdiag(OmegaI(:, :, 1), OmegaI(:, :, 2), OmegaI(:, :, 3));
        
        % objective: maximize rate
        sumR = 0;
        sumr = 0;
        for i = 1 : K
            sumR = sumR + R(:, :, i);
            sumr = sumr + r(:, :, i);
        end
        rate = 0;
        for i = 1 : K
            sumW = sumR - R(:, :, i);
            sumw = sumr - r(:, :, i);

            X = sigma(i) + h(:, i)' * (sumW + Omega) * h(:, i);
            Y = sigma(i) + h(:, i)' * (sumw + omega) * h(:, i);
            rate = rate + log(sigma(i) + h(:, i)' * (sumR + Omega) * h(:, i)) / log(2) ...
                   - log(Y) / log(2) - trace(inv(Y) * (X - Y)) / log(2);
        end
        maximize rate;
        
        subject to
            % capacity contraints & pwer constraints
            for l = 1 : L
                X = E(:, :, l) * (sumR + Omega) * E(:, :, l)';
                Y = E(:, :, l) * (sumr + omega) * E(:, :, l)';

                log(Y) / log(2) + trace(inv(Y) * (X - Y)) / log(2) ...
                - log(E(:, :, l) * Omega * E(:, :, l)') / log(2) <= C(l);
            
                trace(X) <= P(l);
            end
    cvx_end
    
    objective = cvx_optval;
    feasible = false;
    if strfind(cvx_status, 'Solved')
        feasible = true;
    end
end