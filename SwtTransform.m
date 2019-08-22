%函数功能：计算输入图像的笔画宽度图像
%包含两次计算，第二次为笔画拐角处的计算值的修正
function [swtmap2]=SwtTransform(str,searchDirection)
 im=imread(str);
 if size(im,3)==3
     im=rgb2gray(im);
 end
im=im2double((im));
%figure,imshow(im);
edgeMap=edge(im,'canny');
%figure,imshow(edgeMap,[]);

sobelMask = fspecial('sobel');
dx = imfilter(im,sobelMask);
dy = imfilter(im,sobelMask');
[edgePointRows, edgePointCols] = find(edgeMap);

%初始化梯度角度矩阵
theta=zeros(size(edgeMap,1),size(edgeMap,2));
%计算边缘图像中每个像素的梯度方向
for i=1:size(edgeMap,1)
    for j=1:size(edgeMap,2)
        if edgeMap(i,j)==1
            theta(i,j)=atan2(dy(i,j),dx(i,j));
        end
    end
end
[m,n] = size(edgeMap);
% 初始化笔画宽度
swtMap = zeros(m,n);
for i=1:m
    for j=1:n
        swtMap(i,j) =1000;
    end
end

% 笔划宽度最大值
maxStrokeWidth = 350;

% 初始化矩阵储存所有可能的笔画点
strokePointsX = zeros(size(edgePointCols));
strokePointsY = zeros(size(strokePointsX));
sizeOfStrokePoints = 0;
widthvalue=zeros(1,10000);
widthvalue(1)=inf;
% Iterate through all edge points and compute stoke widths
for i=1:size(edgePointRows)
    step = 1;
    initialX = edgePointRows(i);
    initialY = edgePointCols(i);
    isStroke = 0;
    initialTheta = theta(initialX,initialY);
    sizeOfRay = 0;
    pointOfRayX = zeros(maxStrokeWidth,1);
    pointOfRayY = zeros(maxStrokeWidth,1);
    
    % 记录射线的第一个点
    pointOfRayX(sizeOfRay+1) = initialX;
    pointOfRayY(sizeOfRay+1) = initialY;
    
    % 增加射线长度
    sizeOfRay = sizeOfRay + 1;

    while step < maxStrokeWidth
        nextX = round(initialX + cos(initialTheta) * searchDirection * step);
        nextY = round(initialY + sin(initialTheta) * searchDirection * step);
        
        step = step + 1;
        
        
        if nextX < 1 || nextY < 1 || nextX > m || nextY > n
            break
        end
        
        % 记录射线的下一个点
        pointOfRayX(sizeOfRay+1) = nextX;
        pointOfRayY(sizeOfRay+1) = nextY;
        
        %增加射线长度
        sizeOfRay = sizeOfRay + 1;
        
        % 找到对应的边缘像素点
        if edgeMap(nextX,nextY)
            
            oppositeTheta = theta(nextX,nextY);
            
            % 判断梯度方向角度
            if abs(abs(initialTheta - oppositeTheta) - pi) < pi/2
                isStroke = 1;
                strokePointsX(sizeOfStrokePoints+1) = initialX;
                strokePointsY(sizeOfStrokePoints+1) = initialY;
                sizeOfStrokePoints = sizeOfStrokePoints + 1;
            end                   
            
  break           
        end
    end
    
    if isStroke
        
        % 计算笔划宽度
        strokeWidth = sqrt((nextX - initialX)^2 + (nextY - initialY)^2);
        widthvalue(i)=strokeWidth;
        % Iterate all ray points and populate with the minimum stroke width
        for j=1:sizeOfRay
            swtMap(pointOfRayX(j),pointOfRayY(j)) = min(swtMap(pointOfRayX(j),pointOfRayY(j)),strokeWidth);
        end
    end
    
end

a=max(widthvalue,[],2);
[a,b]=size(swtMap);
for i=1:m
    for j=1:n
        if swtMap(i,j)>=0.1*a
            swtMap(i,j)=1000;
        end
    end
end

for i=1:sizeOfStrokePoints
    step = 1;
    initialX = strokePointsX(i);
    initialY = strokePointsY(i);
    initialTheta = theta(initialX,initialY);
    sizeOfRay = 0;
    pointOfRayX = zeros(maxStrokeWidth,1);
    pointOfRayY = zeros(maxStrokeWidth,1);
    swtValues = zeros(maxStrokeWidth,1);
    sizeOfSWTValues = 0;
    
    % 记录射线第一个点
    pointOfRayX(sizeOfRay+1) = initialX;
    pointOfRayY(sizeOfRay+1) = initialY;
   % 增加射线长度
    sizeOfRay = sizeOfRay + 1;
    
    % Record the swt value of first stoke point
    swtValues(sizeOfSWTValues+1) = swtMap(initialX,initialY);
    sizeOfSWTValues = sizeOfSWTValues + 1;
    
    while step < maxStrokeWidth
        nextX = round(initialX + cos(initialTheta) * searchDirection * step);
        nextY = round(initialY + sin(initialTheta) * searchDirection * step);
        
        step = step + 1;
        
        % 记录射线的下一个点
        pointOfRayX(sizeOfRay+1) = nextX;
        pointOfRayY(sizeOfRay+1) = nextY;
        
        % 增加射线长度
        sizeOfRay = sizeOfRay + 1;
        
        % Record the swt value of next stoke point
        swtValues(sizeOfSWTValues+1) = swtMap(nextX,nextY);
        sizeOfSWTValues = sizeOfSWTValues + 1;
        
        %找到对应的像素点
        if edgeMap(nextX,nextY)
            break
        end
    end
    
    % 计算笔画宽度中值
    strokeWidth = median(swtValues(1:sizeOfSWTValues));
    
    % 取最小值
    for j=1:sizeOfRay
        swtMap(pointOfRayX(j),pointOfRayY(j)) = min(swtMap(pointOfRayX(j),pointOfRayY(j)),strokeWidth);
    end
    
end

se=strel('square',2);
swtmap=imopen(swtMap,se);
swtmap2=imclose(swtmap,se);
figure,imshow(swtmap2,[]);

swtmap2=swtmap2/1000;
