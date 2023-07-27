function [x] = smooth_loop(x, Ts, threshold)
% Ts 是loop细分算法中使用的数据结构，蕴含细分网格的邻接信息
n =length(Ts);
model_size = max(max(x)-min(x));
x = x + 2*model_size; % 将平均值带离坐标轴附近

for i=1:n
    x = loop_divide(x, Ts{i});
end
x(abs(x)<threshold) = 0;
for i=n:-1:1
    x = loop_reconstruct(x, Ts{i});
end

x = x - 2*model_size;
end

function x = loop_divide(x, T)
[n_e,~] = size(T);
n = T(n_e,5); n_v = n - n_e;

n_neighbor = zeros(n_v, 1);
for Ti=transpose(T)
    n_neighbor(Ti(1)) = n_neighbor(Ti(1))+1;
    n_neighbor(Ti(2)) = n_neighbor(Ti(2))+1;
end
gamma_v = 8/5*(3/8+1/4*cos(2*pi./n_neighbor).^2);
delta_v = (1 - 8/5*(3/8+1/4*cos(2*pi./n_neighbor).^2))./n_neighbor;
ve_mat = sparse([T(:,1) T(:,2)], [1:n_e 1:n_e], 1, n_v, n_e);

x(1:n_v,:) = x(1:n_v,:) - delta_v.* (ve_mat * x(n_v+1:n,:));
x(1:n_v,:) = x(1:n_v,:) ./ gamma_v;
x(n_v+1:n,1:3) = x(n_v+1:n,1:3) - 3/8* (x(T(:,1),1:3) + x(T(:,2),1:3)) ...
    - 1/8* (x(T(:,3),1:3) + x(T(:,4),1:3)); % 不允许边界情况出现
end

function x = loop_reconstruct(x, T)
[n_e,~] = size(T);
n = T(n_e,5); n_v = n - n_e;

n_neighbor = zeros(n_v, 1);
for Ti=transpose(T)
    n_neighbor(Ti(1)) = n_neighbor(Ti(1))+1;
    n_neighbor(Ti(2)) = n_neighbor(Ti(2))+1;
end
gamma_v = 8/5*(3/8+1/4*cos(2*pi./n_neighbor).^2);
delta_v = (1 - 8/5*(3/8+1/4*cos(2*pi./n_neighbor).^2))./n_neighbor;
ve_mat = sparse([T(:,1) T(:,2)], [1:n_e 1:n_e], 1, n_v, n_e);

x(n_v+1:n,1:3) = x(n_v+1:n,1:3) + 3/8* (x(T(:,1),1:3) + x(T(:,2),1:3)) ...
    + 1/8* (x(T(:,3),1:3) + x(T(:,4),1:3)); % 不允许边界情况出现
x(1:n_v,:) = x(1:n_v,:) .* gamma_v;
x(1:n_v,:) = x(1:n_v,:) + delta_v.* (ve_mat * x(n_v+1:n,:));
end