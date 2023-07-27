clear
im = imread('figures/01.jpg');
im = im2gray(im);
% im = imnoise(im, );
draw_decomp(im, int8(3), 'sym4');

function draw_decomp(im, k, type)
% 画出灰度图im的 k 次二维小波分解结果
% k一般设置为 1~3
assert(isinteger(k));
if k<=0
    return
end
if nargin<3
    type = 'haar';
end
[cA1, cH1, cV1, cD1] = dwt2(im, type);
draw_decomp(cA1, k-1, type); %递归调用

k=double(k);
cA1 = uint8(cA1);
cH1 = uint8(cH1*2.^k);
cV1 = uint8(cV1*2.^k);
cD1 = uint8(cD1*2.^k);

figure
sgtitle(['第' num2str(k) '层分解'])
subplot('Position',[0.1 0.1 0.4 0.4]), imshow(cA1);
subplot('Position',[0.52 0.1 0.4 0.4]), imshow(cH1);
subplot('Position',[0.1 0.52 0.4 0.4]), imshow(cV1);
subplot('Position',[0.52 0.52 0.4 0.4]), imshow(cD1);
end

