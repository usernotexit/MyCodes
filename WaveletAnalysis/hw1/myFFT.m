function Y = myFFT(y)
% 长度为2的幂次的（均匀）数列的离散傅里叶变换
% make sure that the length y is 2^N, where N is a positive integer
N = floor(log2(length(y))); assert(length(y)==2^N);
y = complex(y);
y_ = y;

weights = exp(-1i*pi*2^(-N+1)*(0:(2^N-1)));
for n=0:N-1
    for k = 0:(2^(N-1-n)-1)
        
        for j = 0:(2^n-1)
            idx = 2^n*k+j;
            u = y(idx +1);
            v = y(idx+2^(N-1) +1) * weights(j*2^(N-n-1) +1);
            y_(2^(n+1)*k+j +1) = (u+v);
            y_(2^(n+1)*k+j+2^(n) +1) = (u-v);
        end
    end
    y = y_;
end
Y = y;
end