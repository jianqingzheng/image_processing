clear
%%
folder_path='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data';
path={'D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\red_area.jpg'}; 
% path={'/Users/LengHoufu/Desktop/baf 1_T001 copy.tif'};
%     '/Users/LengHoufu/Desktop/baf 1_T001 copy.tif',
%     '/Users/LengHoufu/Desktop/baf 1_T001 copy.tif'
%     };
thresh0=0.08;         %
thresh1=0.4;         %
erode_rate=6;    %
sigma0=0.4;            %
sigma1=1.2;            %
min_area0=50;       %
min_area1=400;       %

k_size=701;
k=-ones(k_size,k_size);
m=-sum(sum(k))-1;
k((k_size-1)/2,(k_size-1)/2)=m;
kernel=k./m;
%%
close all 
paths=get_dirs(folder_path,'.jpg');
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


% BM=imgaussfilt(B,sigma)>tg;
% se = strel('line',erode_rate,90);
% se = offsetstrel('ball',erode_rate,erode_rate);
se = strel('cube',erode_rate);
BM=imerode(BM,se);

BM=bwareaopen(BM,min_area1);
imshow(cat(3,BM./4,Z,B),[]);

lb_bw=bwlabel(BM,8);
num=max(max(max(lb_bw)))

end

 
 
 
 
 
function [dirs_out]=get_dirs(folder,dirPattern,commonstr,dirs_in)
% get_dirs - Get the directories for required files in one folder and its subfolders
%--------------------------------------------------------------------------
%   get_dirs(folder,dirPattern,commonstr,dirs_in) is a recursive function,
%   which returns a cell including the directories of all the required
%   files.
%--------------------------------------------------------------------------
%   [dirs_out] = get_dirs(folder,dirPattern,commonstr,dirs_in)
%   'dirs_out'  - output cell of the directories of all the required files.
%   'folder'    - a dir string of the folder to search.
%   'dirPattern'- a string of required files' pattern.
%   'commonstr' - a string required in the dirs of those selected files.
%   'dirs_in'   - a cell of dirs for recursion
%--------------------------------------------------------------------------
%   Examples:
%      % Before starting recursion:
%      >> folder='/data/someone/Aneurysm_Seg'
%      >> dirPat='.jpg'
%      >> commonstr={'epoch_1','_2'}
%      >> dirs_out=get_dirs(folder,dirPat,commonstr)
%  
%      % after completing recursion
%      dirs_out =
%      
%      12?ó1 cell array
% 
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_10_2.jpg'
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_11_2.jpg'
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_12_1.jpg'
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_12_2.jpg'
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_13_2.jpg'
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_14_2.jpg'
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_15_2.jpg'
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_16_2.jpg'
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_17_2.jpg'
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_18_2.jpg'
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_19_2.jpg'
%         '/data/someone/Aneurysm_Seg//Trained_1/prediction/epoch_1_2.jpg'
%--------------------------------------------------------------------------
%   MATLAB Ver >= R2016b
%--------------------------------------------------------------------------
%   $ Author: Jachin $
%   $ Revision: 2.1 $  $ Date: 2018/03/13 22:24 $
%--------------------------------------------------------------------------
%   See also: 'fileparts', 'contains'
 
%--------------------------------------------------------------------------
 
%% parameter transfer
if nargin<4
    dirs_in={};
    if nargin<3
        commonstr='';
        if nargin<2
            dirPattern='';
        end
    end
end
%% verify and append
[~,~,pattern]=fileparts(folder);
if (strcmp(pattern,dirPattern)||isempty(dirPattern))
    if ~isempty(commonstr)
        if iscell(commonstr)
            result=0;
            for i=1:size(commonstr,1)
                temp=1;
                for j=1:size(commonstr,2)
                    temp=temp&contains(folder,commonstr{i,j});
                end
                result=result|temp;
            end
        else
            result=contains(folder,commonstr);
        end
    else
        result=1;
    end
    if result
        dirs_out=[dirs_in;folder];
    else
        dirs_out=dirs_in;
    end
else
    dirs_out=dirs_in;
end
%% recursion
dirs_recur=dirs_out;
files=dir(folder);
if contains(folder,'/')
    slash='/';
else
    slash='\';
end
for i=3:size(files,1)
    fn=files(i).name;
    folder_recur=sprintf('%s%s%s', folder,slash, fn);
    dirs_out=get_dirs(folder_recur,dirPattern,commonstr,dirs_recur);
    dirs_recur=dirs_out;
end
end
