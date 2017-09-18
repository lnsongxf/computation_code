function uc = u(c,gamma)
if gamma==1
    uc = log(c);
else
    uc = c.^(1-gamma) / (1-gamma);
end
end