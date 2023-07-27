function px = calChebshev(args, x, varargin)
% 根据算好的插值系数args计算切比雪夫插值在x处的值
% 默认定义域为 [-1, 1] ，通过varargin修改采样范围
% inputs:
%   args    --  插值系数，程序将根据其长度确定插值次数 n
%   x       --  待求值点，可以是矩阵，逐个计算并返回大小的矩阵
%   varargin-- 其他输入，如插值区间（未完成）
assert(length(size(args))<3, 'Invaild args!')

size_x = size(x);% 输入向量化
x = reshape(x, [], 1);

n = length(args) - 1;
Tk = cos(acos(x) * (0:n));

px = args * transpose(Tk);
px = reshape(px, size_x);
end