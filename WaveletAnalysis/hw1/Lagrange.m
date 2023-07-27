function [fx] = Lagrange(x, x0, f)
% 根据插值点列x0计算Lagrange多项式函数fx
% 直接根据基函数计算
S = size(x);
fx = zeros(S);
for xi = x0
    tmp = ones(S);
    for xj = x0
        if xj == xi
            continue;
        end
        tmp = tmp .* (x-xj)./(xi-xj);
    end
    fx = fx + tmp .* f(xi);
end
end