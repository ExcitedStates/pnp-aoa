function str = int0str(num, lg)
if num <0
    str = num2str(num);
else
    if num == 0
        str = repmat('0', 1, lg);
    else
        str = [repmat('0', 1, (lg - floor(log10(num)) - 1) ) num2str(num)];
    end
end
end