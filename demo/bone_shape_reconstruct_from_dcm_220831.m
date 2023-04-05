close all
clear
% set directory
vox_folder=fullfile('data');
subfolder_names=["20220831 irradiation BM chimera C57 HTP"];

% set parameters
% rot_axis =   [1,0,0;     ];
rotx_axis =   [1,0,0     ];
roty_axis =   [0,1,0     ];

angs=      [0,0;
            0,0;
            ];
sec_locs=   [0,350;];

vol_id=1;
subfolder_id=1;
filter_flag=1;
intensity_bias=[0,0,-100,0,400,100,0,0];
thresh_bone0=800+4000+intensity_bias(subfolder_id);
thresh_bone1=1900+4000+intensity_bias(subfolder_id);
thresh_bone2=800+4000+intensity_bias(subfolder_id);
% sec_bias=245;
sec_bias=180;
% sec_bias=250;
% sec_bias=80;

% min_area1=500;
min_area1=100;

remove_small_vol=1;
remove_vol_rank=2;

sigma=4;
sigma_vol=6;
erode_iter=2;
erode_sz=1;
dilate_iter=12;
% dilate_iter=20;
dilate_sz=4;
% dilate_sz=3;

% subfolder_name=subfolder_names(vol_id);
% axis=rot_axis(vol_id,:);
ang=angs(vol_id,:);
% sec_id=[1:sec_locs(vol_id,1),sec_locs(vol_id,1):end];

% dcmfiles = dir();

dirinfo = dir(fullfile(vox_folder,subfolder_names(subfolder_id)));
dirinfo= dirinfo(3:end);
% dcmfiles = cell(length(dirinfo),1);
dcmfiles = cell(0);
for K = 1 : length(dirinfo)
  thisdir = dirinfo(K).name;
  if ~all(ismember('.png',char(thisdir)))
  dcmfiles =[dcmfiles;fullfile(vox_folder,subfolder_names(subfolder_id),thisdir)];
  end
end

for pth_id =1:1:numel(dcmfiles)
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
    vol_shape2=(img>thresh_bone2);
    vol_shape2=bwareaopen(vol_shape2,min_area1./5);
end

%% filter0
if remove_small_vol
    vol_shape0=keep_remove_ranked_vol(vol_shape0,[-remove_vol_rank]);
end

if filter_flag %%%%%%
    if remove_small_vol
        vol_shape2=keep_remove_ranked_vol(vol_shape2,[-remove_vol_rank]);
    end

    %% filter1
    vol_shape1 = serial_dilate_erode_vol(vol_shape1,kron([-erode_sz,erode_sz],ones(1,erode_iter)),vol_shape0,min_area1);
    vol_shape1=imgaussfilt3(double(vol_shape1),sigma_vol)>0.3;
    vol_shape1 = serial_dilate_erode_vol(vol_shape1,kron([erode_sz],ones(1,7)),vol_shape0,min_area1);
    vol_shape2=fill_slide_centroid(vol_shape2);
    vol_shape2 = serial_dilate_erode_vol(vol_shape2,kron([dilate_sz],ones(1,dilate_iter)),1,0,'cube');
    vol_shape2=imgaussfilt3(double(vol_shape2),sigma_vol)>0.1;
    vol_shape2 = serial_dilate_erode_vol(vol_shape2,kron([-dilate_sz],ones(1,dilate_iter)),1,min_area1*200,'cube');
end%%%%%%%%%%
%% select volume
% calculate the ratio of cancellous bone divided by the volume inside the
% cortical bone.
vol_shape=vol_shape0;
if filter_flag
    vol_shape=vol_shape.*(vol_shape2);
    vol_low_density=sum(sum(sum(vol_shape0.*vol_shape2)));
    vol_total=sum(sum(sum(vol_shape2)));
    ratio=vol_low_density./vol_total
end

%% plot
figure(pth_id);
% vol_shape(:,:,1)=0;
% vol_shape(:,:,end)=0;
vs=volshow(double(vol_shape));
set(vs,'Renderer', 'Isosurface');
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
