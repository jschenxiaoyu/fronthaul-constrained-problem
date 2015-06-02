function a = DecBinConverter(number, L)
    i = 0;
    a = zeros(1, L);
    while (0 ~= number)
        i = i + 1;
        a(i) = mod(number, 2);
        number = floor(number / 2);
    end
%     k = size(a, 2);
%     b = a;
%     for i = 1 : k
%         a(i) = b(k - i + 1);
%     end
end