function y = Fbox(q, p)
    y = (0 < q & q < 1) & (0 < p & p < 1);
    y = double(y);
end