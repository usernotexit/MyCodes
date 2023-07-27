clear;
X = imread('figures/01.jpg');
subplot(141);image(X); axis off         %画出原始图像
title('原始图像');

%% haar
[c,s]=wavedec2(X,2,'haar');
sizec=size(c);
%对分解系数进行处理以突出轮廓部分，弱化细节部分
for i=1:sizec(2)    %从1到第二个值c
   if(c(i)>400)     %对低频系数分段增强
       %增强系数
      c(i)=2*c(i);
   else             %对高频系数进行弱化
       %弱化系数
      c(i)=0.5*c(i);
   end
end
XX=uint8(waverec2(c,s,'haar'));        %对处理后的系数进行重构
%画出重构后的图像
subplot(142);image(XX); axis off
title('Haar');

%% sym4
[c,s]=wavedec2(X,2,'sym4');
sizec=size(c);
%对分解系数进行处理以突出轮廓部分，弱化细节部分
for i=1:sizec(2)    %从1到第二个值c
   if(c(i)>300)     %对低频系数分段增强
       %增强系数
      c(i)=2*c(i);
   else             %对高频系数进行弱化
       %弱化系数
      c(i)=0.5*c(i);
   end
end
XX=uint8(waverec2(c,s,'sym4'));        %对处理后的系数进行重构
%画出重构后的图像
subplot(143);image(XX); axis off
title('sym4');


%% db3
[c,s]=wavedec2(X,2,'db3');
sizec=size(c);
%对分解系数进行处理以突出轮廓部分，弱化细节部分
for i=1:sizec(2)    %从1到第二个值c
   if(c(i)>300)     %对低频系数分段增强
       %增强系数
      c(i)=2*c(i);
   else             %对高频系数进行弱化
       %弱化系数
      c(i)=0.5*c(i);
   end
end
XX=uint8(waverec2(c,s,'db3'));        %对处理后的系数进行重构
%画出重构后的图像
subplot(144);image(XX); axis off
title('db3');
