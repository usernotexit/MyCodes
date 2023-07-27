%% 验证 定理2
f = @(x) abs(sin(6*x)).^3 - cos(5 * exp(x));
%% calculate V and get baseline
h = 0.0001;
X = (-1+3e-6):h:1;
Y = diff(f(X))/h;   % 一阶导数
Z = diff(Y)/h;      % 二阶导数
W = diff(Z)/h;      % 三阶导数
T = diff(W);      % 四阶 差分

%V = inte(abs(T), h); % 复化梯形方法求W的总变差 V
%V = sum(abs(T));
V = 37684 + 1296*2*3; % mathematica 分段积分 + 间断点
baseLine_2 = @(n) 4/pi*V /2 ./(n-2).^2;
baseLine_3 = @(n) 4/pi*V /3 ./(n-3).^3;

%% plot the errs
% 比较系数]
figure;
M = 10000;    % number of grids
m = 3:10;
for N = 2.^m
    args_Cheb = getChebshevArgs(f, N); % 计算切比雪夫插值系数
    x_all = (-1+1e-5):(2/M):1;
    x_all_ = 2 * rand([1,M]) - 1;
    
    yi = f(x_all);
    yi_Chev = calChebshev(args_Cheb, x_all);
    yi_ = f(x_all_);
    yi_Chev_ = calChebshev(args_Cheb, x_all_);
    
    err_Chev = max(abs(yi - yi_Chev));
    err_Chev_ = max(abs(yi_ - yi_Chev_));

    plot(log2(N), log10(err_Chev), '.', 'Color', 'g','MarkerSize',10);
    plot(log2(N), log10(err_Chev_), '.', 'Color', 'b','MarkerSize',10);
    hold on
end
plot(m, log10(baseLine_2(2.^m)))
plot(m, log10(baseLine_3(2.^m)))
xlabel('log2 N')
ylabel('log10 error')

function F = inte(W, h)
% 复化梯形方法
    n = length(W);
    weight = [1, repmat([2],1, n-2), 1];
    F = sum(W.*weight) * h/2;
end