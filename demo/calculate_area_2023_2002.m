clear
close all 
%% Parameter setting
sigma_smooth=20;
sigma_low=2;
% thresh_low=0.2;
thresh_low=0.1;
min_conn_area=10;

% folder_path='G:\Nastia_test dataset';
folder_path='G:\AA37_test';

% suffix='.tif';
suffix='.TIF';

files=cell(0);

files =[files;folder_path];
%%
for f =1:numel(files)
    img_pths=get_dirs(files{f},suffix);
    img=0;
    for i=1:numel(img_pths)%
    img_tmp=double(read(Tiff(img_pths{i},'r')));
    img_tmp(img_tmp<1)=max(max(max(max(img_tmp(img_tmp>1)))));

    % Filter the image using a high passing filter
    img_tmp=img_normalize(img_tmp);
    img_low_pass=img_tmp;
    img_low_pass=imgaussfilt(img_low_pass,sigma_smooth);
    img_high_pass=img_normalize(img_tmp-img_low_pass);
    img_tmp=img_high_pass;
    img_low_pass=img_tmp;
    img_low_pass=imgaussfilt(img_low_pass,sigma_low);
    img_high_pass=img_normalize(abs(img_tmp-img_low_pass));

    % Plot the image filtered by a high passing filter
    figure(1)
    imshow(img_high_pass,[]);
    
    % plot the mask of the cells
    figure(2)
    msk_tmp = img_high_pass>thresh_low;
    msk_cells = any(msk_tmp,3);
    imshow(double(msk_cells),[]);
    
    % Fill the small black background area to find the gap
    figure(3)
    dila_erod_sz=[10,-20,35,-50,-30,30,10];
    msk_interest=serial_dilate_erode_vol(msk_cells,dila_erod_sz);

    % Find all the connected white components in the image
    CC = bwconncomp(~msk_interest);
    
    % Get the number of pixels in each component
    numPixels = cellfun(@numel, CC.PixelIdxList);
    
    % Find the component with the largest number of pixels
    [~, maxIdx] = max(numPixels);
    
    % Get the center of the largest component
    center = regionprops(CC, 'Centroid');
    center = center(maxIdx).Centroid;
    
    % Define the size of the cropped area
    org_img_sz=size(img_tmp);
    img_size = [org_img_sz(1), org_img_sz(2)/2]; % [height, width]
    
    % Calculate the coordinates of the top-left corner of the cropped area
    x1 = center(1) - img_size(2)/2;
    y1 = center(2) - img_size(1)/2;
    
    % Crop the area of interest
    img_msked=zeros(size(img_tmp));
    img_msked(:,int16(x1):int16(x1+img_size(2)),:)=img_high_pass(:,int16(x1):int16(x1+img_size(2)),:);
    
    % Remove the small connected area    
    msk_cells = any(img_msked>thresh_low./1,3);
    msk_cells = serial_dilate_erode_vol(msk_cells,[4,-3]);
    msk_cells = bwareaopen(msk_cells,min_conn_area,8);
    msk_cells = ~bwareaopen(~msk_cells,min_conn_area,8);
    
    % Display the cropped image    
    imshow(msk_cells,[])
    
    % Write the cropped binary mask of the cells
    imwrite(msk_cells,[img_pths{i}(1:end-4),'.png'],'png');
%     areas(i)=sum(sum(msk_cells));
%     ratio(i)=sum(sum(msk_cells))/numel(msk_cells).*3;
    areas(i)=numel(msk_cells)./2-sum(sum(msk_cells));
    ratio(i)=sum(sum(msk_cells))/numel(msk_cells).*2;

    end

end
table=[img_pths,num2cell(areas'),num2cell(ratio')];
writecell(table,[folder_path,'_stat.xlsx'],'Sheet',1)%,'Range','D1')
 




%% Function list

function out_vol=serial_dilate_erode_vol(vol,dilate_erode_rates)
for r= dilate_erode_rates
%     se=strel('sphere',abs(r));
    se=strel('cube',abs(r));
    if r>0
        vol=imdilate(vol,se);%.*vol_shape0;
    else
        vol=imerode(vol,se);%.*vol_shape0;
    end
end
out_vol=vol;
end

function [nm_img]=img_normalize(img)
img_min=min(min(min(img)));
img_max=max(max(max(img)));
nm_img=(img-img_min)./(img_max-img_min+1e-5);
end

function [mask]= get_purple(img,thresh)
mask=(img(:,:,1)>thresh(1))&(img(:,:,2)./img(:,:,1)<thresh(2))&(img(:,:,3)./img(:,:,1)<thresh(3));
end
function [mask]= get_red(img,thresh)
mask=(img(:,:,1)>thresh(1))&(img(:,:,2)./img(:,:,1)<thresh(2))&(img(:,:,3)./img(:,:,1)<thresh(3));
end
function [mask]= get_blue(img,thresh)
mask=(img(:,:,1)./img(:,:,3)<thresh(1))&(img(:,:,2)./img(:,:,3)<thresh(2))&(img(:,:,3)>thresh(3));
end
function [mask]= get_green(img,thresh)
mask=(img(:,:,1)./img(:,:,2)<thresh(1))&(img(:,:,2)>thresh(2))&(img(:,:,3)./img(:,:,2)<thresh(3));
end
function [mask]=get_white(img,thresh)
img_intensity_sum=sum(img,3);
mask=(img(:,:,1)>thresh(4))&(img(:,:,2)>thresh(4))&(img(:,:,3)>thresh(4))&(img(:,:,1)*3./img_intensity_sum>thresh(1))&(img(:,:,2)*3./img_intensity_sum>thresh(2))&(img(:,:,3)*3./img_intensity_sum>thresh(3));
end

function [mask]=circ_mask(img_size,mask_size)
mask=zeros(img_size);
% [xr,yr]=mask_size;
[idx,idy]=meshgrid([1:img_size(1)],[1:img_size(2)]);
mask(((idx-img_size(1)/2)./mask_size(1)).^2+((idy-img_size(2)/2)./mask_size(1)).^2<=1)=1;
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

function proc_img=white_bal(img,weights)
R=img(:,:,1);G=img(:,:,2);B=img(:,:,3);
RAVG=mean(mean(R))*weights(1);
GAVG=mean(mean(G))*weights(2);
BAVG=mean(mean(B))*weights(3);
KAVG=(RAVG+GAVG+BAVG)/3;
proc_img=cat(3,KAVG*R/RAVG,KAVG*G/GAVG,KAVG*B/BAVG);
end



