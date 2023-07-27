%% 验证 定理3
f = @(x) 1./(1 + 25*x.^2) - sin(20*x);
b = 0.1987;
rho = b + sqrt(1+b^2); %(1+sqrt(26))/5; % ???
M = abs(f(b*1i)); % f在椭圆内的界，显然真正的上界比M大
%% baseline
baseLine = @(n) 4*M/(rho-1) *rho.^(-n);

%% plot the errs
% 比较系数]
figure;
M = 10000;    % number of grids
m = 3:8;
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

    plot(N, log10(err_Chev), '.', 'Color', 'g','MarkerSize',10);
    plot(N, log10(err_Chev_), '.', 'Color', 'b','MarkerSize',10);
    hold on
end

plot(2.^m, log10(baseLine(2.^m)))
xlabel('N')
ylabel('log10 error')