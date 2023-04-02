close all
clear
% vox_folder='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\vols\';
% subfolder_names=["C1 L_1_20220303_102301"];
% vox_folder='E:\houfu\';
% subfolder_names=[   
%                 "20220619 TRAF3 KO experiment TRAF3 KO group",...
%                 "20220525 aged mice micro-CT",...
%                 "20220303 C57 HTP MM bone new",...
%                 "20220215 micro-CT Ana Houfu",...
%                 "20220210 C57 NZB HTP n=15",...
%                 "20211124 NZB mice bone microCT",...
%                 "20211004 Ana UZC aged mice tibia microCT"
%                 ];
vox_folder='F:\collab\houfu\';
% subfolder_names=["20221118 High fat diet bone scan"];
subfolder_names=["20220831 irradiation BM chimera C57 HTP"];

rot_axis=   [1,0,0;     ];
angs=       [0;         ];
sec_locs=   [250,350;   ];

vol_id=1;
group_id=1;
thresh_bone0=1200+4000;
thresh_bone1=2300+4000;
thresh_bone2=1100+4000;

% thresh_bone0=1200;
% thresh_bone1=2300;
% thresh_bone2=1100;

min_area1=10000;
% remove_small_vol=1;
remove_small_vol=0;
remove_vol_rank=2;
remove_vol_sz=[1000000,1500000];

sigma=4;
sigma_vol=6;
erode_iter=2;
erode_sz=1;
dilate_iter=12;
dilate_sz=4;

% subfolder_name=subfolder_names(vol_id);
axis=rot_axis(vol_id,:);
ang=angs(vol_id);
% sec_id=[1:sec_locs(vol_id,1),sec_locs(vol_id,1):end];

% dcmfiles = dir();

dirinfo = dir(vox_folder+subfolder_names(group_id));
dirinfo= dirinfo(3:end);
dcmfiles = cell(length(dirinfo),1);
for K = 1 : length(dirinfo)
  thisdir = dirinfo(K).name;
  dcmfiles{K} = vox_folder+subfolder_names(group_id)+'\'+thisdir;
end

% vox_files=get_dirs(vox_folder,'.dcm');
% sz=size(dicomread(vox_files{vn}));
% num_files=numel(vox_files);
% img=zeros(num_files,sz(1),sz(2));
% for vn=1:num_files
% [img(vn,:,:), map] = dicomread(vox_files{vn});
% % info = dicominfo(vox_files{vn})
% end
for pth_id =1:numel(dcmfiles)

[img,spatial,dim] = dicomreadVolume(dcmfiles{pth_id});
img=img(:,:,:);
img=imrotate3(img,ang,axis);

%% thresholding
vol_shape0=(img>thresh_bone0);
vol_shape0=bwareaopen(vol_shape0,min_area1./10);
vol_shape1=(imgaussfilt3(img,sigma)>thresh_bone1);
vol_shape1=bwareaopen(vol_shape1,min_area1);
vol_shape2=(imgaussfilt3(img,sigma)>thresh_bone2);
vol_shape2=bwareaopen(vol_shape2,min_area1);
%% filter0
% vol_shape=bwareaopen(vol_shape,min_area0);
if remove_small_vol
    L = bwconncomp(vol_shape0,26);% 
    stats = regionprops(L,'Area');
    Ar = cat(1, stats.Area);
    [~,ind]=sort(Ar,'descend');
    LM = labelmatrix(L);
%     vol_shape(find(LM~=ind(2)))=0;%
    vol_shape0(find(LM==ind(remove_vol_rank)))=0;%
end
if remove_small_vol
    L = bwconncomp(vol_shape2,26);% 
    stats = regionprops(L,'Area');
    Ar = cat(1, stats.Area);
    [~,ind]=sort(Ar,'descend');
    LM = labelmatrix(L);
%     vol_shape(find(LM~=ind(2)))=0;%
    vol_shape2(find(LM==ind(remove_vol_rank)))=0;%
end


%% filter1
se=strel('sphere',erode_sz);
for i=1:erode_iter
    vol_shape1=imerode(vol_shape1,se);
end
%% section
% section
img(:,:,sec_locs(vol_id,2):size(img,1))=0;
img(:,:,1:sec_locs(vol_id,1))=0;
vol_shape1(:,:,sec_locs(vol_id,2):size(vol_shape1,1))=0;
vol_shape1(:,:,1:sec_locs(vol_id,1))=0;
vol_shape0(:,:,sec_locs(vol_id,2):size(vol_shape0,1))=0;
vol_shape0(:,:,1:sec_locs(vol_id,1))=0;
vol_shape2(:,:,sec_locs(vol_id,2):size(vol_shape2,1))=0;
vol_shape2(:,:,1:sec_locs(vol_id,1))=0;


vol_shape1=bwareaopen(vol_shape1,min_area1);
for i=1:erode_iter
    vol_shape1=imdilate(vol_shape1,se).*vol_shape0;
end
vol_shape1=imgaussfilt3(double(vol_shape1),sigma_vol)>0.3;
% vol_shape=vol_shape1;
for i=1:7
    vol_shape1=imdilate(vol_shape1,se).*vol_shape0;
end
vol_shape1=bwareaopen(vol_shape1,min_area1);
% =================================
% vol2
% vol_shape2(:,:,sec_locs(vol_id,2):size(vol_shape1,1))=0;
% vol_shape2(:,:,1:sec_locs(vol_id,1))=0;
% vol_shape2=bwareaopen(vol_shape2,min_area1);
se=strel('sphere',dilate_sz);
for i=1:dilate_iter
    vol_shape2=imdilate(vol_shape2,se);%.*vol_shape0;
end
vol_shape2=imgaussfilt3(double(vol_shape2),sigma_vol)>0.1;
for i=1:dilate_iter+erode_iter+10
    vol_shape2=imerode(vol_shape2,se);
end

% for i=1:3
%     vol_shape2=imdilate(vol_shape1,se);%.*vol_shape0;
% end
vol_shape2=bwareaopen(vol_shape2,min_area1);


%% select volume
% vol_shape=vol_shape0.*vol_shape1;
vol_shape=vol_shape0.*vol_shape2;
% vol_shape=vol_shape0.*~vol_shape2;
% vol_shape=vol_shape2;
% vol_shape=vol_shape0.*~vol_shape1;
% vol_shape=bwareaopen(vol_shape,min_area1);
vol_low_density=sum(sum(sum(vol_shape0.*vol_shape2)));
vol_total=sum(sum(sum(vol_shape2)));

ratio=vol_low_density./vol_total
%% plot
intensity = [0 20 40 120 220 1024];
alpha = [0 0.1 0.15 0.3 0.38 0.5];
% alpha = [0 0 0.15 0.3 0.38 0.5].*.7;
color = ([0 0 0; 43 0 0; 103 37 20; 199 155 97; 216 213 201; 255 255 255])/ 255;
queryPoints = linspace(min(intensity),max(intensity),256);
alphamap = interp1(intensity,alpha,queryPoints)';
colormap = interp1(intensity,color,queryPoints);
ViewPnl = uipanel(figure,"Title","4-D Dicom Volume");

% figure(pth_id);
volshow(vol_shape);
% volshow(img.*uint16(vol_shape),"Colormap",colormap,"Alphamap",alphamap,"Parent",ViewPnl,'BackgroundColor',[0 0.4470 0.7410]);
% volshow(img.*uint16(vol_shape0.*~vol_shape1),"Colormap",colormap,"Alphamap",alphamap,"Parent",ViewPnl,'BackgroundColor',[0 0.4470 0.7410]);
% volshow(img.*uint16(vol_shape0.*vol_shape1),"Colormap",colormap,"Alphamap",alphamap,"Parent",ViewPnl,'BackgroundColor',[0 0.4470 0.7410]);
end
% volshow(vol_shape);
% lighting none
% colormap copper
% niftiwrite(int8(vol_shape),[vox_folder,subfolder_name,'_vol_shape.nii']);
%%
function [predirs,foldernames]=get_predirs(dirs_in,layernumber)
predirs=dirs_in;
predirs_temp=predirs;
for i=0:layernumber
    if i==layernumber
        foldernames_temp=cell(size(dirs_in));
    end
    for j=1:numel(predirs)
        [predir,foldername,~]=fileparts(predirs_temp{j});
        if i==layernumber
            foldernames_temp{j}=foldername;
        else
            predirs_temp{j}=predir;
        end
    end    
end
if layernumber>0
    predirs=predirs_temp;
    
end
foldernames=foldernames_temp;
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
%      12×1 cell array
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
% if contains(folder,'/')
%     slash='/';
% else
    slash='\';
% end
for i=3:size(files,1)
    fn=files(i).name;
    folder_recur=sprintf('%s%s%s', folder,slash, fn);
    dirs_out=get_dirs(folder_recur,dirPattern,commonstr,dirs_recur);
    dirs_recur=dirs_out;
end
    
end