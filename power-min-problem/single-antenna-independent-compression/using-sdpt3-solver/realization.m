clc;
clear all;
close all;

realizatoin_times = 10;

for i = 1 : realizatoin_times
    [col results] = convergence(i);
    if (1 == i)
        mini = ones(1, col) * 1e4;
        iter = zeros(1, col);
        totaltime = zeros(1, col);
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
        totaltime(j) = totaltime(j) + results.totaltime_cvx(j);
    end
end

% compute average cvx time for each iretation
one_iter_time = [];
for i = 1 : col
    one_iter_time(i) = totaltime(i) / iter(i);
end

for i = 1 : col
    if (0 == num(i))
        mini(i) = NaN;
    end

    iter(i) = iter(i) / realizatoin_times;
    totaltime(i) = totaltime(i) / realizatoin_times;
end

x = results.x;
plot(x, mini, 'bo-');
grid on;

clear i j col results
save('result.mat');
