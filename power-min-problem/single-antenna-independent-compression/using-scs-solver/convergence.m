function [col results] = convergence(realizetion_times)
    L = 3;  % number of base stations
    K = 3;  % number of users
    Nl = 1; % number of antennas of BS
    gama = 2; % SNR threshold for each user
    sigma = ones(K, 1); % variance of noise

    % CSI matrix h
    for i = 1 : 3
        for j = 1 : 3
            if (i == j)
                h(i, j) = 1;
            else
                h(i, j) = 0.5;
            end
        end
    end

    E = zeros(Nl, Nl * L, L);
    for i = 1 : L
        E(:, (i - 1) * Nl + 1 : i * Nl ,i) = eye(Nl);
    end

    params = [];
    params.L = L;
    params.K = K;
    params.Nl = Nl;
    params.gama = gama;
    params.sigma = sigma;
    params.h = h;
    params.E = E;

    point = 0;
    results = [];
%     results.mini = [];
%     results.iter_times = [];
%     results.totaltime_cvx = [];
    x = 2 : 1 : 20;
    for z = x  % range of fronthaul capacity
        tmark = tic;
        % initialize omega and r
        r = rand(Nl * L, Nl * L, K);
        for i = 1 : K
            r(:, :, i) = r(:, :, i) * r(:, :, i)';
        end
%         clear omega;
%         for i = 1 : L
%             omega(:, :, i) = rand(Nl, Nl);
%             omega(:, :, i) = omega(:, :, i) * omega(:, :, i)';
%         end
%         omega = blkdiag(omega(:, :, 1), omega(:, :, 2), omega(:, :, 3));
        omega = rand(Nl * L, Nl * L);
        omega = omega * omega';

        C = z * ones(L, 1);
        params.C = C;

        point = point + 1;
        iteration = 0;
        obj = [];
        time_cvx = [];
        totaltime_cvx = 0;
        totaltime_solve = 0;
        totaltime_trans = 0;
        while 1
            iteration = iteration + 1;
            params.r = r;
            params.omega = omega;
            tic;
            [Omega R objective feasible solving_time] = power_minization_cvx(params);
            time_cvx(iteration) = toc;
            totaltime_cvx = totaltime_cvx + time_cvx(iteration);
            totaltime_solve = totaltime_solve + solving_time;
            if (false == feasible)
                break;
            end
            omega = Omega;
            r = R;

            obj(iteration) = objective;
            if ((iteration > 1) && ...
                (abs(obj(iteration) - obj(iteration - 1)) / obj(iteration) < 1e-6))
                break;
            end
        end
        totaltime_trans = totaltime_cvx - totaltime_solve;

		results.totaltime_cvx(point) = totaltime_cvx;
        results.totaltime_solve(point) = totaltime_solve;
        results.totaltime_trans(point) = totaltime_trans;
        results.iter_times(point) = iteration;
        if (true == feasible)
            results.mini(point) = obj(iteration);
        else
            results.mini(point) = 0;
        end
        disp(sprintf('%dth realization, the capacity is: %d', realizetion_times, z));
        toc(tmark);
    end
    results.x = x;
    col = point;
end