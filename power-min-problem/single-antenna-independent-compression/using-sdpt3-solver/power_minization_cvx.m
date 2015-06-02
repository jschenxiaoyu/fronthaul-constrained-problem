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

    cvx_begin quiet
        cvx_solver('sdpt3');
        cvx_solver_settings('max_iters', 50000, 'eps', 1e-6);
        variable OmegaI(Nl, Nl, L) hermitian semidefinite
        variable R(Nl * L, Nl * L, K) hermitian semidefinite

        Omega = blkdiag(OmegaI(:, :, 1), OmegaI(:, :, 2), OmegaI(:, :, 3));

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
            % capacity contraint:
            for l = 1 : L
                X = E(:, :, l) * (sumR + Omega) * E(:, :, l)';
                Y = E(:, :, l) * (sumr + omega) * E(:, :, l)';

                log(Y) / log(2) + trace(inv(Y) * (X - Y)) / log(2) ...
                - log(E(:, :, l) * Omega * E(:, :, l)') / log(2) <= C(l);
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