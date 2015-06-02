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
        variable Omega(Nl * L, Nl * L) hermitian semidefinite
        variable R(Nl * L, Nl * L, K) hermitian semidefinite
        
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
            Y = real(sigma(i) + h(:, i)' * (sumw + omega) * h(:, i));
            rate = rate + log(sigma(i) + h(:, i)' * (sumR + Omega) * h(:, i)) / log(2) ...
                   - log(Y) / log(2) - trace(inv(Y) * (X - Y)) / log(2);
        end
        maximize rate;
        
        subject to
		    %capacity constraints
            number = 0;
            while (number < 7)
                number = number + 1;
                a = DecBinConverter(number, L);
                capacity = 0;
                Es = [];
                for i = 1 : L
                    if (0 ~= a(i))
                        X = E(:, :, i) * (sumR + Omega) * E(:, :, i)';
                        Y = E(:, :, i) * (sumr + omega) * E(:, :, i)';
                        capacity = capacity + log(Y) / log(2) + trace(inv(Y) * (X - Y)) / log(2);
                        Es = [Es; E(:, :, i)];

                        % power constraint
                        if (1 == sum(a))
                            trace(X) <= P(i);
                        end
                    end
                end
                capacity = capacity - log_det(Es * Omega * Es') / log(2);
                capacity <= sum(a) * C(1);
            end
    cvx_end
    
    objective = cvx_optval;
    feasible = false;
    if strfind(cvx_status, 'Solved')
        feasible = true;
    end
end