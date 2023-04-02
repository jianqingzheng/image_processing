clear
%%
% folder_path='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data';
% folder_path='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\houfu\20211013 BM OC TRAP aged HTP mice copy\';
% folder_path='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\houfu\20211013 BM OC TRAP aged HTP mice copy\plate 2';

% folder_path='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\20220228 Expt1 and Expt2\20220309 Expt1-1';
% folder_path='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\20220228 Expt1 and Expt2\20220309 Expt1-2';
% folder_path='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\20220228 Expt1 and Expt2\20220309 Expt1-3';
% folder_path='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\20220228 Expt1 and Expt2\20220309 Expt2-1';
% folder_path='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\20220228 Expt1 and Expt2\20220309 Expt2-2';
% folder_path='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\202207\tiff20220705';
folder_path='D:\Users\JachinZ\Documents\MATLAB\bio_image_process\data\20220725\QuPath_staining July';%\20220210 n=15 histology22-62';


% dirinfo = dir(folder_path);
% dirinfo= dirinfo(3:end);
files=cell(0);
% for K = 1 : numel(dirinfo)
%   thisdir = dirinfo(K).name;
% %   if ~ismember('.tif',char(thisdir))
%   files =[files;folder_path,'\',thisdir];
% %   end
% end
% 
files =[files;folder_path];

% path={'/Users/LengHoufu/Desktop/baf 1_T001 copy.tif',
%     '/Users/LengHoufu/Desktop/baf 1_T001 copy.tif',
%     '/Users/LengHoufu/Desktop/baf 1_T001 copy.tif'
%     };
thresh0=0.1;         %
thresh1=0.7;         %
thresh2=0.7;
ch_w=[0.9,1.,1.3];
circ_mask_size=[490,490];
%%
close all 
% img_pths=glob([folder_path,'*.jpg']);
% areas=cell(numel(img_pths),3);
% areas=zeros(numel(img_pths),1);
% ratio=zeros(numel(img_pths),1);
sigma=3;
thresh_white=[0.8,0.8,0.8,0.5];
thresh_red=[0.3,0.5,0.5];
thresh_green=[];
thresh_max=50;
for f =1:numel(files)
    img_pths=get_dirs(files{f},'.tif');
    img=0;
for i=1:1%numel(img_pths)%
    img_tmp=double(read(Tiff(img_pths{i},'r')));
% %     finfo=czifinfo(img_pths{i});
%     img_tmp=ReadImage6D(img_pths{i}, true, 1);
%     img_tmp(img_tmp>thresh_max)=thresh_max;
%     img_tmp=normalize(img_tmp, 'scale');
    img_tmp=img_normalize(img_tmp);
    imshow(img_tmp,[]);
%     mask_black=get_white(img_tmp,[0.5,0.5,0.5,0.01]);
%     img_tmp=histeq(img_tmp).*mask_black;
%     for s =1:3
%     img_tmp(:,:,s)=normalize(img_tmp(:,:,s), 'range');
%     end

    img_tmp=img_tmp.*(~get_red(img_tmp,thresh_red));
    mask_white=get_white(img_tmp,thresh_white);
    img_tmp(:,:,2:3)=img_tmp(:,:,2:3).*(~mask_white);
%     img_tmp=imgaussfilt(img_tmp,sigma);
img=img+img_tmp;
end
% img=normalize(img, 'range');
% imshow(img./numel(img_pths),[]);%
imshow(img,[]);%





% figure(2*i-1)
% imshow(img,[]);
img_max=256;%max(max(max(img)));
img=double(img)/double(img_max);
% img=white_bal(img,ch_w);
% imshow(img,[]);
% figure(1);
R=img(:,:,1);
G=img(:,:,2);
B=img(:,:,3);
% mask=double((R>G)&(R>0.8*B)&(B>1.*G)&(R>thresh0)&(G<thresh1)&(B<thresh2)).*circ_mask(size(R),circ_mask_size);
filtered_img=cat(3,mask,img(:,:,2),img(:,:,3));
% figure(2*i)
% imshow(filtered_img,[]);

imwrite(filtered_img,[img_pths{i}],'png');
areas(i)=sum(sum(mask));
ratio(i)=sum(sum(mask))/numel(mask);
end
table=[img_pths,num2cell(areas),num2cell(ratio)];
writecell(table,[folder_path,'_stat.xlsx'],'Sheet',1)%,'Range','D1')
 

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
% function [mask]=get_white(img,thresh)
% mask=(img(:,:,1)>thresh)&(img(:,:,2)>thresh)&(img(:,:,3)>thresh);
% end

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

%%
function r = bfGetReader(varargin)
% BFGETREADER return a reader for a microscopy image using Bio-Formats
% 
% SYNOPSIS  r = bfGetReader()
%           r = bfGetReader(path)
%
% Input 
%
%    id - (Optional - string) A valid path to the microscopy image
%
%    stichFiles (Optional - scalar). Toggle the grouping of similarly
%    named files into a single dataset based on file numbering.
%    Default: false;
%
% Output
%
%    r - A reader object of class extending loci.formats.ReaderWrapper
%
% Adapted from bfopen.m

% OME Bio-Formats package for reading and converting biological file formats.
%
% Copyright (C) 2012 - 2014 Open Microscopy Environment:
%   - Board of Regents of the University of Wisconsin-Madison
%   - Glencoe Software, Inc.
%   - University of Dundee
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

% Input check
ip = inputParser;
ip.addOptional('id', '', @ischar);
ip.addOptional('stitchFiles', false, @isscalar);
ip.parse(varargin{:});
id = ip.Results.id;

% verify that enough memory is allocated
bfCheckJavaMemory();

% load the Bio-Formats library into the MATLAB environment
status = bfCheckJavaPath();
assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
    'to the static Java path or add it to the Matlab path.']);

% Prompt for a file if not input
isFile = exist(id, 'file') == 2;
isFake = ischar(id) && strcmp(id(end-3:end), 'fake');
if nargin == 0 || (~isFile && ~isFake)
    [file, path] = uigetfile(bfGetFileExtensions, 'Choose a file to open');
    id = [path file];
    if isequal(path, 0) || isequal(file, 0), return; end
elseif isFile && ~verLessThan('matlab', '7.9')
    [~, f] = fileattrib(id);
    id = f.Name;
end

% set LuraWave license code, if available
if exist('lurawaveLicense')
    path = fullfile(fileparts(mfilename('fullpath')), 'lwf_jsdk2.6.jar');
    javaaddpath(path);
    java.lang.System.setProperty('lurawave.license', lurawaveLicense);
end

r = loci.formats.ChannelFiller();
r = loci.formats.ChannelSeparator(r);
if ip.Results.stitchFiles
    r = loci.formats.FileStitcher(r);
end

OMEXMLService = loci.formats.services.OMEXMLServiceImpl();
r.setMetadataStore(OMEXMLService.createOMEXMLMetadata());
r.setId(id);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File       : GetOMEData.m
% Version    : 1.3
% Author     : czsrh
% Date       : 13.09.2019
% Institution : Carl Zeiss Microscopy GmbH
%
% Simple script to get "some" not, but not all metainformation using
% the MATLAB wrapper for the BioFormats library.
%
% Use at your own risk.
%
% Copyright(c) 2019 Carl Zeiss AG, Germany. All Rights Reserved.
%
% Permission is granted to use, modify and distribute this code,
% as long as this copyright notice remains part of the code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function OMEData = GetOMEData(filename)
% Get OME Meta Information using BioFormats Library

    % To access the file reader without loading all the data, use the low-level bfGetReader.m function:
    reader = bfGetReader(filename);

    % You can then access the OME metadata using the getMetadataStore() method:
    omeMeta = reader.getMetadataStore();

    % get ImageCount --> currently only reading one image is supported
    imagecount = omeMeta.getImageCount();
    % create empty cell array to store the image IDs
    imageIDs_str = cell(1, imagecount);
    imageIDs = cell(1, imagecount);

    % try to get all the imageIDs as strings and numbers (zero-based)
    try
        for id = 1:imagecount
            imageIDs_str{id} = omeMeta.getImageID(id-1);
            imageIDs{id} = id-1; 
        end
        % store ine OMEData
        OMEData.ImageIDs = imageIDs;
        OMEData.ImageIDstrings = imageIDs_str;
    catch
        OMEData.ImageIDs = 'na';
        OMEData.ImageIDstrings = 'na';
        msg = 'No suitable ImageIDs found.';
        warning(msg);
    end

    % use default imageID to read other metadata
    imageID = imageIDs{1};

    % get the actual metadata and store them in a structured array
    [pathstr,name,ext] = fileparts(filename);
    OMEData.FilePath = pathstr;
    OMEData.Filename = strcat(name, ext);

    % Get dimension order
    OMEData.DimOrder = char(omeMeta.getPixelsDimensionOrder(imageID).getValue());

    % Number of series inside the complete data set
    OMEData.SeriesCount = reader.getSeriesCount();

    % Dimension Sizes C - T - Z - X - Y
    OMEData.SizeC = omeMeta.getPixelsSizeC(imageID).getValue();
    OMEData.SizeT = omeMeta.getPixelsSizeT(imageID).getValue();
    OMEData.SizeZ = omeMeta.getPixelsSizeZ(imageID).getValue();
    OMEData.SizeX = omeMeta.getPixelsSizeX(imageID).getValue();
    OMEData.SizeY = omeMeta.getPixelsSizeY(imageID).getValue();

    % Scaling XYZ
    try
        OMEData.ScaleX = round(double(omeMeta.getPixelsPhysicalSizeX(imageID).value()),3); % in micron
    catch
        msg = 'Problem getting X-Scaling. Use Default = 1';
        warning(msg);
        OMEData.ScaleX = 1;
    end

    try
        OMEData.ScaleY = round(double(omeMeta.getPixelsPhysicalSizeY(imageID).value()),3); % in micron
    catch
        msg = 'Problem getting Y-Scaling. Use Default = 1';
        warning(msg);
        OMEData.ScaleY = 1;
    end


    try
        OMEData.ScaleZ = round(double(omeMeta.getPixelsPhysicalSizeZ(imageID).value()),3); % in micron
    catch
        % in case of only a single z-plane set to 1 micron ...
            msg = 'Problem getting Z-Scaling. Use Default = 1';
        warning(msg);
        OMEData.ScaleZ = 1;
    end

    % read relevant objective information from metadata
    try
        % get the correct objective ID (the objective that was used to acquire the image)
        tmp = char(omeMeta.getInstrumentID(imageID));
        OMEData.InstrumentID = str2double(tmp(end));
        tmp = char(omeMeta.getObjectiveSettingsID(OMEData.InstrumentID));
        objID = str2double(tmp(end));
        % error handling --> sometime only one objective is there with ID > 0
        numobj = omeMeta.getObjectiveCount(OMEData.InstrumentID);
        if numobj == 1
            objID = 0;
        end

        OMEData.ObjID = objID; 
    catch
            msg = 'No suitable instrument and objective ID found.';
            warning(msg);
    end

    try
        % get objective immersion
        OMEData.ObjImm = char(omeMeta.getObjectiveImmersion(OMEData.InstrumentID, OMEData.ObjID).getValue());
    catch
        msg = 'Problem getting immersion type.';
        warning(msg);
        OMEData.ObjImm = 'na';
    end

    try
        % get objective lens NA
        OMEData.ObjNA = round(omeMeta.getObjectiveLensNA(OMEData.InstrumentID, OMEData.ObjID).doubleValue(),2);
    catch
        msg = 'Problem getting objective NA.';
        warning(msg);
        OMEData.ObjNA = 'na';
    end

    try
        % get objective magnification
        OMEData.ObjMag = round(omeMeta.getObjectiveNominalMagnification(OMEData.InstrumentID, OMEData.ObjID).doubleValue(),2); 
    catch
        msg = 'Problem getting objective magnification.';
        warning(msg);
        OMEData.ObjMag = 'na';
    end

    try
        % get objective model
        OMEData.ObjModel = char(omeMeta.getObjectiveModel(OMEData.InstrumentID, OMEData.ObjID));
    catch
        msg = 'Problem getting objective model.';
        warning(msg);
        OMEData.ObjModel = 'na';
    end

    % get excitation and emission wavelengths for all channels
    for c = 1:OMEData.SizeC
        try
            OMEData.WLEx{c} = round(omeMeta.getChannelExcitationWavelength(imageID, c-1).value().doubleValue());
            OMEData.WLEm{c} = round(omeMeta.getChannelEmissionWavelength(imageID, c-1).value().doubleValue());
        catch
            %msg = 'Problem getting excitation and emission wavelengths. Set to zero.';
            %warning(msg);
            OMEData.WLEx{c} = 'na';
            OMEData.WLEm{c} = 'na';
        end

        try
            OMEData.Channels{c} = char(omeMeta.getChannelName(imageID, c-1));
            OMEData.Dyes{c} = char(omeMeta.getChannelFluor(imageID, c-1));
        catch
            msg = 'No Metadata for current channel available.';
            warning(msg);
            OMEData.Channels{c} = 'na';
            OMEData.Dyes{c} = 'na';
        end

        try
            OMEData.LaserIntensity{c} = omeMeta.getChannelLightSourceSettingsAttenuation(imageID, c-1);
        catch
            msg = 'Readout of Intensity Values failed.';
            warning(msg);
            OMEData.LaserIntensity{c} = 'na';
        end

    end

    % get the number of instruments
    OMEData.NumberOfInstruments = omeMeta.getInstrumentCount();
    %OMEData.NumberOfInstruments = 'na';
    OMEData.InstrumentID = {};
    OMEData.NumberofLightsources = {};
    OMEData.LaserID = {};
    OMEData.LaserPower = {};

    % get all instrument IDs
    for num = 1: OMEData.NumberOfInstruments
        OMEData.InstrumentID{num} = omeMeta.getInstrumentID(num-1);
        OMEData.NumberofLightsources{num} = omeMeta.getLightSourceCount(num-1);
        if omeMeta.getLightSourceCount(num-1) == 0
            msg = join(['No LightSource found in MetaData for Indstrument: ', num2str(num)]);
            warning(msg);
        end
    end

    for num = 1:OMEData.NumberOfInstruments
        try
            OMEData.LaserID{num} = omeMeta.getLaserID(num, OMEData.NumberofLightsources{num});
        catch
            OMEData.LaserID{num} = 'na';
        end

        try
            OMEData.LaserPower{num} = omeMeta.getLaserPower(num, OMEData.NumberofLightsources{num});
        catch
            OMEData.LaserPower{num} = 'na';
        end

    end

    % close BioFormats Reader
    reader.close()
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File       : ReadImage6D.m
% Version    : 1.2
% Author     : czsrh
% Date       : 13.09.2019
% Institution : Carl Zeiss Microscopy GmbH
%
% Simple script to read CZI image data incl metaiformation using
% the MATLAB wrapper for the BioFormtas library.
% The actual image pixel data with be stored inside a 6D array.
% Limitations: - For tile images one has to define the actual image series
%              - Dimension Order : XYCZT
%
% Use at your own risk.
%
% Copyright(c) 2019 Carl Zeiss AG, Germany. All Rights Reserved.
%
% Permission is granted to use, modify and distribute this code,
% as long as this copyright notice remains part of the code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Simple usage example:
out = ReadImage6D(filename, true, 1);
metadata = out{2};
image6d = out{1};
img = image6d(1,1,1,1,:,:);
figure()
img = squeeze(img);
imagesc(img);
axis equal
axis tight
%}


function out = ReadImage6D(filename, useSeriesID, seriesID)
    switch nargin
        case 1
            useSeriesID = true;
            seriesID = 1;
    end

    % Get OME Meta-Information
    MetaData = GetOMEData(filename);

    % The main inconvenience of the bfopen.m function is that it loads all
    % the content of an image regardless of its size.
    
    % Initialize BioFormtas Reader
    reader = bfGetReader(filename);

    % add progress bar
    h = waitbar(0,'Processing Data ...');
    totalframes = MetaData.SeriesCount * MetaData.SizeC * MetaData.SizeZ * MetaData.SizeT;
    framecounter = 0;

    % Pre-allocate array with size (Series, SizeT, SizeZ, SizeC, SizeX, SizeY) 
    image6d = zeros(MetaData.SeriesCount, MetaData.SizeT, MetaData.SizeZ,  MetaData.SizeC, MetaData.SizeY, MetaData.SizeX);

    if useSeriesID == true
        series = seriesID;
    end
    
    if useSeriesID == false
        series = MetaData.SeriesCount;
    end

    % for series = 1 : MetaData.SeriesCount
    for series = 1 : series

        % set reader to current series
        reader.setSeries(series-1);
        for timepoint = 1: MetaData.SizeT
            for zplane = 1: MetaData.SizeZ
                for channel = 1: MetaData.SizeC

                    framecounter = framecounter + 1;
                    % update waitbar
                    wstr = {'Reading Images: ', num2str(framecounter), ' of ', num2str(totalframes),'Frames' };
                    waitbar(framecounter / totalframes, h, strjoin(wstr))

                    % get linear index of the plane (1-based)
                    iplane = loci.formats.FormatTools.getIndex(reader, zplane - 1, channel - 1, timepoint -1) +1;
                    % get frame for current series
                    image6d(series, timepoint, zplane, channel, :, :) = bfGetPlane(reader, iplane);

                end
            end
        end
    end

    % close waitbar
    close(h)

    % close BioFormats Reader
    reader.close();

    % store image data and meta information in cell array
    out = {};
    % store the actual image data as 6d array
    out{1} = image6d;
    % store the image metainformation
    out{2} = MetaData;
    
end
    
%%

