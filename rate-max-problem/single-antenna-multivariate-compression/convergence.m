function [col results] = convergence()
    clc;
    clear;

    L = 3;  % number of base stations
    K = 3;  % number of users
    Nl = 1; % number of antennas of BS
    P = 100 * ones(L, 1); % power constraint for each BS
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

    params.L = L;
    params.K = K;
    params.Nl = Nl;
    params.P = P;
    params.sigma = sigma;
    params.h = h;
    params.E = E;

    point = 0;
    results = [];
    for z = 2 : 1 : 20  % range of fronthaul capacity
        % initialize omega and r
        r = zeros(Nl * L, Nl * L, K);
        clear omega;
        for i = 1 : L
            omega = rand(Nl * L, Nl * L);
            omega = omega * omega';
        end

        C = z * ones(L, 1);
        params.C = C;

        point = point + 1;
        iteration = 0;
        obj = [];
        while 1
            iteration = iteration + 1;
            params.r = r;
            params.omega = omega;
            [Omega R objective feasible] = rate_maximization_cvx(params);
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
        if (true == feasible)
            results(point) = obj(iteration);
        else
            results(point) = 0;
        end
    end
    col = point;
end