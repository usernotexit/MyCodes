% 切比雪夫插值
f1 = @(x) abs(sin(6*x)).^3 - cos(5 * exp(x));
f2 = @(x) 1./(1 + 25*x.^2) - sin(20*x);

for f = {f1, f2}
    f = f{1,1};
%% 绘图
M = 10000;    % number of grids
id = 0;
fprintf('***************** new function ********************\n');
figure;
for N = [8, 16, 32, 64, 128]
    args_Cheb = getChebshevArgs(f, N); % 计算切比雪夫插值系数
    x_all = (-1+1e-5):(2/M):1;
    x_all_ = 2 * rand([1,M]) - 1;
    
    yi = f(x_all);
    yi_Chev = calChebshev(args_Cheb, x_all);
    yi_ = f(x_all_);
    yi_Chev_ = calChebshev(args_Cheb, x_all_);
    
    err_Chev = max(abs(yi - yi_Chev));
    err_Chev_ = max(abs(yi_ - yi_Chev_));
    fprintf("N = %d, ", N);
    %fprintf("Max Error of grid: %e\n", err_Chev);      % 切比雪夫点
    fprintf("Max Error of grid: %e\t%e\n", err_Chev, err_Chev_);      % 切比雪夫点
    
    id = id + 1;
    subplot(3,2,id);
    plot(x_all, yi)
    hold on
    plot(x_all, yi_Chev)
    hold on
    %legend("Origin","Chebyshev")
    xlabel("X-Axis");
    ylabel('Y-Axis');
    title(['N = ' num2str(N)])
end
end