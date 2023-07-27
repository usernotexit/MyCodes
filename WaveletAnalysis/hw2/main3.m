clear;
X = imread('figures/02.jpg');
subplot(131);image(X); axis off         %画出原始图像
title('原始图像');
%用小波函数sym4对X进行2层小波分解
%c为各层分解系数,s为各层分解系数长度,也就是大小
[c,s]=wavedec2(X,2,'sym4');
sizec=size(c);      %sizec为[1xc]的行向量。

%对分解系数进行处理以突出轮廓部分，弱化细节部分
for i=1:sizec(2)    %从1到第二个值c
   if(c(i)>500)     %对低频系数分段增强
       %增强系数
      c(i)=1.6*c(i);
   else            
       if (c(i)>200)
           c(i) = 1.3*c(i);
       else %对高频系数进行弱化
       %弱化系数
            c(i)=0.8*c(i);
       end
   end
end
XX=uint8(waverec2(c,s,'sym4'));        %对处理后的系数进行重构
%画出重构后的图像
subplot(132);image(XX); axis off
title('增强图像');

subplot(133); image(uint8(1.6*X));axis off
title('直接提高亮度');