%函数功能：提取样本集的HOG特征并训练相应的SVM分类器
%分类器保存为svm_model文件夹内的model
%本项目已经训练好了一个model
imdsTrain = imageDatastore('E:\百度云\图像处理\text-detection\HOG_SVM',...  
    'IncludeSubfolders',true,...  
    'LabelSource','foldernames'); 
Train_disp = countEachLabel(imdsTrain);
disp(Train_disp);
imageSize = [256,256];% 对所有图像进行此尺寸的缩放 ?
image1 = readimage(imdsTrain,1); 
scaleImage = imresize(image1,imageSize); 
[features, visualization] = extractHOGFeatures(scaleImage); 
imshow(scaleImage);hold on; plot(visualization) 
numImages = length(imdsTrain.Files); 
featuresTrain = zeros(numImages,size(features,2),'single'); % featuresTrain为单精度 ?
for i = 1:numImages 
 imageTrain = readimage(imdsTrain,i); 
 imageTrain = imresize(imageTrain,imageSize); 
 featuresTrain(i,:) = extractHOGFeatures(imageTrain); 
end 
trainLabels = imdsTrain.Labels;
classifer = fitcsvm(featuresTrain,trainLabels,'KernelFunction','linear'); 
save('./svm_model/model','classifer');
