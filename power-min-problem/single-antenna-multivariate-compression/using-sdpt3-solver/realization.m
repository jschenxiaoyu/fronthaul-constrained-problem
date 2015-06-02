clc;
clear;

for i = 1 : 100
    [col results] = convergence;
    if (1 == i)
        sum = zeros(1, col);
        num = zeros(1, col);
    end
    
    for j = 1 : col
        if (0 ~= results(j))
            sum(j) = sum(j) + results(j);
            num(j) = num(j) + 1;
        end
    end
end

for i = 1 : col
    if (0 ~= num(j))
        sum(i) = sum(i) / num(i);
    else
        sum(i) = 0;
    end
end

save('result.mat');
