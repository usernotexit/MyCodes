clear;
X = imread('figures/01.jpg');
subplot(131);image(X); axis off         %画出原始图像
title('原始图像');
%% 两层分解
%用小波函数sym4对X进行2层小波分解
[c,s]=wavedec2(X,2,'sym4');
sizec=size(c);      %sizec为[1xc]的行向量。
for i=1:sizec(2)    %从1到第二个值c
   if(c(i)>200)     %对低频系数分段增强
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
title('增强图像（两层分解）');

%% 三层分解
%用小波函数sym4对X进行2层小波分解
[c,s]=wavedec2(X,3,'sym4');
sizec=size(c);      %sizec为[1xc]的行向量。
for i=1:sizec(2)    %从1到第二个值c
   if(c(i)>200)     %对低频系数分段增强
      %增强系数
      c(i)=2*c(i);
   else            %对高频系数进行弱化
      %弱化系数
      c(i)=0.5*c(i);
   end
end
XX=uint8(waverec2(c,s,'sym4'));        %对处理后的系数进行重构
%画出重构后的图像
subplot(133);image(XX); axis off
title('增强图像（三层分解）');