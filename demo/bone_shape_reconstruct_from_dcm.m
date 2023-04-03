close all
clear
% vox_folder='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\vols\';
% subfolder_names=["C1 L_1_20220303_102301"];
% vox_folder='E:\houfu\';
vox_folder='F:\collab\houfu\';

% subfolder_names=[
%                 "20220620 WT",...
%                 "20220619 TRAF3 KO experiment TRAF3 KO group",...
%                 "20220525 aged mice micro-CT",...
%                 "20220303 C57 HTP MM bone new",...
%                 "20220215 micro-CT Ana Houfu",...
%                 "20220210 C57 NZB HTP n=15",...
%                 "20211124 NZB mice bone microCT",...
%                 "20211004 Ana UZC aged mice tibia microCT"
%                 ];
% subfolder_names=["20221118 High fat diet bone scan"];
subfolder_names=["20220831 irradiation BM chimera C57 HTP"];


rot_axis=   [1,0,0;     ];
rotx_axis=   [1,0,0     ];
roty_axis=   [0,1,0     ];
% angs=      [15,0;
%             0,0;
%             ];
angs=      [0,0;
            0,0;
            ];
% sec_locs=   [250,350;
%              200,300;
%              200,300;
%              240,340;
%              230,330;
%             ];
sec_locs=   [0,350;
            ];
start_id=21;%%%
vol_id=1;
group_id=1;
filter_flag=1;
intensity_bias=[0,0,-100,0,400,100,0,0];
thresh_bone0=800+4000+intensity_bias(group_id);
thresh_bone1=1900+4000+intensity_bias(group_id);
thresh_bone2=800+4000+intensity_bias(group_id);
% sec_bias=245;
% sec_bias=240;
% sec_bias=250;
sec_bias=80;
min_area1=500;
% remove_small_vol=1;
remove_small_vol=0;
remove_vol_rank=2;
remove_vol_sz=[1000000,1500000];

sigma=4;
sigma_vol=6;
erode_iter=2;
erode_sz=1;
dilate_iter=12;
% dilate_iter=20;
dilate_sz=4;
% dilate_sz=3;

% subfolder_name=subfolder_names(vol_id);
axis=rot_axis(vol_id,:);
ang=angs(vol_id,:);
% sec_id=[1:sec_locs(vol_id,1),sec_locs(vol_id,1):end];

% dcmfiles = dir();

dirinfo = dir(vox_folder+subfolder_names(group_id));
dirinfo= dirinfo(3:end);
% dcmfiles = cell(length(dirinfo),1);
dcmfiles = cell(0);
for K = 1 : length(dirinfo)
  thisdir = dirinfo(K).name;
  if ~ismember('.tif',char(thisdir))
  dcmfiles =[dcmfiles;vox_folder+subfolder_names(group_id)+'\'+thisdir];
  end
end

% vox_files=get_dirs(vox_folder,'.dcm');
% sz=size(dicomread(vox_files{vn}));
% num_files=numel(vox_files);
% img=zeros(num_files,sz(1),sz(2));
% for vn=1:num_files
% [img(vn,:,:), map] = dicomread(vox_files{vn});
% % info = dicominfo(vox_files{vn})
% end

for pth_id =start_id:2:numel(dcmfiles)
    pth_id
    close all;
[img,spatial,dim] = dicomreadVolume(dcmfiles{pth_id});
img=img(:,:,:);
if pth_id>size(sec_locs,1)
    ang=angs(end,:);
    sec_loc=sec_locs(end,:);
else
    sec_loc=sec_locs(pth_id,:);
    ang=angs(pth_id,:);
end

img=imrotate3(img,ang(1),rotx_axis);
img=imrotate3(img,ang(2),roty_axis);
%% section
if pth_id>size(sec_locs,1)
    sec_loc=sec_locs(end,:);
else
    sec_loc=sec_locs(pth_id,:);
end
[~,~,z_coord]=ind2sub(size(img),find(img>thresh_bone1));
sec_loc(2)=max(z_coord)-sec_bias;
sec_loc(1)=sec_loc(2)-100;
% sec_loc(1)=0;

img(:,:,sec_loc(2):size(img,1))=[];
img(:,:,1:sec_loc(1))=[];
%% thresholding
vol_shape0=(img>thresh_bone0);
vol_shape0=bwareaopen(vol_shape0,min_area1./10);
if filter_flag
    vol_shape1=(imgaussfilt3(img,sigma)>thresh_bone1);
    vol_shape1=bwareaopen(vol_shape1,min_area1);
    % vol_shape2=(imgaussfilt3(img,sigma)>thresh_bone2);
    vol_shape2=(img>thresh_bone2);
    vol_shape2=bwareaopen(vol_shape2,min_area1./5);
end
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

if filter_flag %%%%%%
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
    % section
    % vol_shape1(:,:,sec_locs(vol_id,2):size(vol_shape1,1))=0;
    % vol_shape1(:,:,1:sec_locs(vol_id,1))=0;
    
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
    % se=strel('sphere',dilate_sz);
    vol_shape2=fill_slide_centroid(vol_shape2);
    se=strel('cube',dilate_sz);
    for i=1:dilate_iter
        vol_shape2=imdilate(vol_shape2,se);%.*vol_shape0;
    end
    vol_shape2=imgaussfilt3(double(vol_shape2),sigma_vol)>0.1;
    for i=1:dilate_iter%+erode_iter+12
    % for i=1:dilate_iter+20
        vol_shape2=imerode(vol_shape2,se);
    end
    
    % for i=1:3
    %     vol_shape2=imdilate(vol_shape1,se);%.*vol_shape0;
    % end
    vol_shape2=bwareaopen(vol_shape2,min_area1*200);

end%%%%%%%%%%
%% select volume
% vol_shape=vol_shape0.*vol_shape1;
vol_shape=vol_shape0;
% vol_shape=vol_shape0.*vol_shape2.*(~vol_shape1);
% vol_shape=vol_shape0.*~vol_shape2;
% vol_shape=vol_shape2;
% vol_shape=vol_shape0.*~vol_shape1;
% vol_shape=bwareaopen(vol_shape,min_area1);
if filter_flag
    vol_shape=vol_shape.*(vol_shape2);
vol_low_density=sum(sum(sum(vol_shape0.*vol_shape2)));
vol_total=sum(sum(sum(vol_shape2)));
ratio=vol_low_density./vol_total
end
%% plot
intensity = [0 20 40 120 220 1024];
alpha = [0 0.1 0.15 0.3 0.38 0.5];
% alpha = [0 0 0.15 0.3 0.38 0.5].*.7;
% color = ([0 0 0; 43 0 0; 103 37 20; 199 155 97; 216 213 201; 255 255 255])/ 255;
color = ([0 0 0; 43 43 43; 103 103 103; 199 199 199; 216 216 216; 255 255 255])/ 255;
queryPoints = linspace(min(intensity),max(intensity),256);
alphamap = interp1(intensity,alpha,queryPoints)';
colormap = interp1(intensity,color,queryPoints);
ViewPnl = uipanel(figure,"Title","4-D Dicom Volume");

figure(pth_id);
% vol_shape = squeeze(vol_shape);
vol_shape(:,:,1)=0;
% vol_shape(:,:,end)=0;
vs=volshow(double(vol_shape));
set(vs,'Renderer', 'Isosurface');
% set(vs,'Renderer', 'VolumeRendering');
% set(vs, 'colormap', colormap,'alphamap',alphamap)
vs.BackgroundColor='w';
vs.IsosurfaceColor=[1.,1.,1.];
vs.Isovalue=.5;
% vs.Lighting=1;
% vs.CameraViewAngle
% vs.CameraPosition
% vs.CameraUpVector
% vs.CameraTarget
vs.CameraTarget=[-0.1,-0.1,0];
vs.CameraViewAngle=10;
% vs.CameraPosition = [-0.3509   -1.2034    3.6760];
% vs.CameraUpVector= [0.0850    0.2791    0.9565];

vs.CameraPosition = [-1.5098    0.4677    3.5471];
vs.CameraUpVector= [-0.1328   -0.0914    0.9869];

% vs.CameraPosition = [0.5,0.5,.5];
% set(vs, 'CameraPosition', [0.,0.,5.]);

print(dcmfiles(pth_id)+'_view0'+'.png', '-dpng', '-r900');
% print(dcmfiles(pth_id)+'_view0'+'.tif', '-depsc','-tiff', '-r900');
% I = getframe(gcf); 
% imwrite(I.cdata, dcmfiles{pth_id}+'_view0'+'.tif', 'TIFF','Resolution',900);
% [indI,cm] = rgb2ind(I.cdata,256);
% imwrite(indI, cm, dcmfiles(pth_id)+'_view0'+'.tif', 'tif','Resolution',30000);
% volshow(img.*uint16(vol_shape),"Colormap",colormap,"Alphamap",alphamap,"Parent",ViewPnl,'BackgroundColor',[0 0.4470 0.7410]);
% volshow(img.*uint16(vol_shape0.*~vol_shape1),"Colormap",colormap,"Alphamap",alphamap,"Parent",ViewPnl,'BackgroundColor',[0 0.4470 0.7410]);
% volshow(img.*uint16(vol_shape0.*vol_shape1),"Colormap",colormap,"Alphamap",alphamap,"Parent",ViewPnl,'BackgroundColor',[0 0.4470 0.7410]);
end
% volshow(vol_shape);
% lighting none
% colormap copper
% niftiwrite(int8(vol_shape),[vox_folder,subfolder_name,'_vol_shape.nii']);
%%

function [out_vol]=fill_slide_centroid(in_vol)
vol_sum3=sum(sum(in_vol,2),1);
[Y,X,Z]=meshgrid(1:size(in_vol,2),1:size(in_vol,1),1:size(in_vol,3));
cX=sum(sum(X.*in_vol,1),2)./vol_sum3;
cY=sum(sum(Y.*in_vol,1),2)./vol_sum3;
cZ=sum(sum(Z.*in_vol,1),2)./vol_sum3;
out_vol=in_vol;
idx=round(permute(sub2ind(size(out_vol),cX,cY,cZ),[3,1,2]));
out_vol(idx)=1;
end


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