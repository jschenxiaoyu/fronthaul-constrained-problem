function [Omega R objective feasible] = power_minization_cvx(params)
    L = params.L;
    K = params.K;
    Nl = params.Nl;
    gama = params.gama;
    sigma = params.sigma;
    h = params.h;
    E = params.E;
    r = params.r;
    C = params.C;
    omega = params.omega;
    
    cvx_begin
        variable Omega(Nl * L, Nl * L) hermitian semidefinite
        variable R(Nl * L, Nl * L, K) hermitian semidefinite
        
        % objective: minimize power
        sumR = 0;
        sumr = 0;
        for i = 1 : K
            sumR = sumR + R(:, :, i);
            sumr = sumr + r(:, :, i);
        end
        power = trace(sumR + Omega);
        minimize power;
        
        subject to
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
                    end
                end
                if (Nl == size(Es, 1))
                    capacity = capacity - log(Es * Omega * Es') / log(2);
                else
                    capacity = capacity - log_det(Es * Omega * Es') / log(2);
                end
%                 capacity = capacity - log_det(Es * Omega * Es') / log(2);
                capacity <= sum(a) * C(1);
            end

            % SINR constraints
            for i = 1 : K
                h(:, i)' * R(:, :, i) * h(:, i) >= ...
                gama * (h(:, i)' * (sumR - R(:, :, i) + Omega) * h(:, i) + sigma(i));
            end
    cvx_end
    
    objective = cvx_optval;
    feasible = false;
    if strfind(cvx_status, 'Solved')
        feasible = true;
    end
end