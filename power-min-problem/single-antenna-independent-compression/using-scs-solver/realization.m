clc;
clear all;
close all;

realization_times = 10;

for i = 1 : realization_times
    [col results] = convergence(i);
    if (1 == i)
        mini = ones(1, col) * 1e4;
        iter = zeros(1, col);
        totaltime_cvx = zeros(1, col);
        totaltime_solve = zeros(1, col);
        totaltime_trans = zeros(1, col);
        num = zeros(1, col);
    end

    for j = 1 : col
        if (0 ~= results.mini(j))
            if (mini(j) > results.mini(j))
                mini(j) = results.mini(j);
            end
            num(j) = num(j) + 1;
        end
        iter(j) = iter(j) + results.iter_times(j);
        totaltime_cvx(j) = totaltime_cvx(j) + results.totaltime_cvx(j);
        totaltime_solve(j) = totaltime_solve(j) + results.totaltime_solve(j);
        totaltime_trans(j) = totaltime_trans(j) + results.totaltime_trans(j);
    end
end

% compute average cvx time for each iteration
one_iter_time_cvx = [];
one_iter_time_solve = [];
one_iter_time_trans = [];
for i = 1 : col
    one_iter_time_cvx(i) = totaltime_cvx(i) / iter(i);
    one_iter_time_solve(i) = totaltime_solve(i) / iter(i);
    one_iter_time_trans(i) = totaltime_trans(i) / iter(i);
end

for i = 1 : col
    if (0 == num(i))
        mini(i) = NaN;
    end

    iter(i) = iter(i) / realization_times;
    totaltime_cvx(i) = totaltime_cvx(i) / realization_times;
    totaltime_solve(i) = totaltime_solve(i) / realization_times;
    totaltime_trans(i) = totaltime_trans(i) / realization_times;
end

x = results.x;
plot(x, mini, 'bo-');
grid on;

clear i j col results
save('result.mat');
