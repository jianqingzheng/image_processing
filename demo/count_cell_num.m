close all 
clear
%% Directory

subfolder_names=["20220831 irradiation BM chimera C57 HTP"];

folder_path=fullfile('data');
paths={'D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\cells\123.jpg'}; 
% path={'/Users/LengHoufu/Desktop/baf 1_T001 copy.tif'};
%     '/Users/LengHoufu/Desktop/baf 1_T001 copy.tif',
%     '/Users/LengHoufu/Desktop/baf 1_T001 copy.tif'
%     };

% paths=get_dirs(folder_path,'.jpg');

%% Parameter setting
thresh0=0.08;         %
thresh1=0.4;         %
erode_rate=6;    %
sigma0=0.4;            %
sigma1=1.2;            %
min_area0=50;       %
min_area1=400;       %

k_size=701;


%% For loop calculation
k=-ones(k_size,k_size);
m=-sum(sum(k))-1;
k((k_size-1)/2,(k_size-1)/2)=m;
kernel=k./m;
areas=cell(numel(paths),3);
for i=1:numel(paths)
img=imread(paths{i});
% figure(1);
% imshow(img,[]);
R=img(:,:,1);
B=img(:,:,3);
Z=zeros(size(R));
% imshow(cat(3,Z,Z,B),[]);
B=double(B)./256;

figure(i)
BM=imgaussfilt(B,sigma0);
BM=conv2(BM,kernel,'same');
% B=reshape(normalize(reshape(B,[1024*1024,1]),'scale','std'),[1024,1024]);
BM=double(BM>thresh0);
% se = strel('ball',dilation_rate,dilation_rate);
% B=imdilate(double(B>thresh),se);
BM=imgaussfilt(BM,sigma1)>thresh1;
BM=double(BM);
% BM=conv2(BM,kernel,'same');
BM=bwareaopen(BM,min_area0);
se = strel('cube',erode_rate);
BM=imerode(BM,se);

BM=bwareaopen(BM,min_area1);
imshow(cat(3,BM./4,Z,B),[]);

lb_bw=bwlabel(BM,8);
num=max(max(max(lb_bw)))

end
