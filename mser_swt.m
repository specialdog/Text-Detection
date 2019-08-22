%%
%提取MSER区域并保存
%mser.tif为二值化的mser区域图像
im=imread('000078.jpg');
im=im2double(rgb2gray(im));
[m,n]=size(im);
%检测MSER区域
[mserregions,mserconcomp]=detectMSERFeatures(im,'ThresholdDelta',4);
figure,imshow(im);
hold on
plot(mserregions, 'showPixelList', true,'showEllipses',false)
mserstats=regionprops(mserconcomp,'all');
box=vertcat(mserstats.BoundingBox);
region=mserconcomp.PixelIdxList;
region=vertcat(region{:});
im2=zeros(size(im));
im2(region)=im(region);
figure,imshow(~im2);
imwrite(~im2,'mser.tif');
%%
%在mser区域的基础上计算笔画宽度图像
%保存为swtmap.tif
str='mser.tif';
swtmap2=SwtTransform(str,1);
imwrite(swtmap2,'swtmap.tif');
%%
%计算原图的笔画宽度图像
%保存为swtmap2.tif
str2='000078.jpg';
swtmap2=SwtTransform(str2,1);
imwrite(swtmap2,'swtmap2.tif');

