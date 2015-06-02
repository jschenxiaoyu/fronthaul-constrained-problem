function [Omega R objective feasible solving_time] = power_minization_cvx(params)
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

    time_pat_cvx_scs = 'Timing: Total solve time: (?<total>[\d\.e[+-]\d]+)';

    cvx_begin %sdp
        cvx_solver('scs');
        cvx_solver_settings('max_iters', 50000, 'eps', 1e-6);
        %cvx_solver_settings('max_iters', 50000, 'eps', 1e-6, 'use_indirect', 1);
        variable OmegaI(Nl, Nl, L) hermitian semidefinite
        variable R(Nl * L, Nl * L, K) hermitian semidefinite

        Omega = blkdiag(OmegaI(:, :, 1), OmegaI(:, :, 2), OmegaI(:, :, 3));
%         for i = 1 : K
%            R(:, :, i) >= 0;
%         end
%         for i = 1 : L
%             OmegaI(:, :, i) >= 0;
%         end

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
    output = evalc('cvx_end');

    timing = regexp(output, time_pat_cvx_scs, 'names');
    solving_time = str2num(timing.total);
    
    objective = cvx_optval;
    feasible = false;
    if strfind(cvx_status, 'Solved')
        feasible = true;
    end
end