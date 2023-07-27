function args = getChebshevArgs(f, n, varargin)
% 对于函数f，做切比雪夫插值，次数为n。
% 默认定义域为 [-1, 1] ，通过varargin修改采样范围
% 4/29 更新：采用自己写的二次幂长度的fft，请确保n是2的幂次。
% inputs:
%   f    -- 待插值函数，如 @(x) sin(x)
%   n    -- 插值多项式次数，注意插值点数为 n+1
%   varargin -- 其他输入，如插值区间（未完成）
xi = cos((0:(1/n):(2-1/n))*pi);
yi = f(xi);
%ci = real(fft(yi));
ci = real(myFFT(yi));

args = ci(1:n+1)./n;
args(1) = args(1)/2;
args(n+1) = args(n+1)/2;
end