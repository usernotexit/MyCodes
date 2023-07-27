clear;
X = imread('figures/02.jpg');
subplot(131);image(X); axis off         %画出原始图像
title('原始图像');
%用小波函数sym4对X进行2层小波分解
%c为各层分解系数,s为各层分解系数长度,也就是大小
[c,s]=wavedec2(X,2,'sym4');
sizec=size(c);      %sizec为[1xc]的行向量。

%% 分段线性函数
c1 = c;
%对分解系数进行处理，
for i=1:sizec(2)    %从1到第二个值c
   if(c(i)>300)     %对低频系数分段增强
      %增强系数
      c(i)=2*c(i);
   else            %对高频系数进行弱化
      %弱化系数
      c(i)=0.5*c(i);
   end
end
XX=uint8(waverec2(c,s,'sym4'));        %对处理后的系数进行重构
%画出重构后的图像
subplot(132);image(XX); axis off
title('增强图像（分段线性）');

%% 光滑函数
%对分解系数进行处理，
theta = 300;
c1 = 4*c1.*c1./(c1+theta);

XX1 = uint8(waverec2(c1,s,'sym4'));
subplot(133); image(XX1);axis off
title('增强图像（光滑函数）');