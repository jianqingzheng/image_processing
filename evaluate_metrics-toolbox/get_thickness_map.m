function D = get_thickness_map(matrix)
% matrix is a 3D matrix representing the volume
% with the thin plate defined by a boundary value of 1

% Calculate the distance transform of the boundary points
bvol=matrix>0.5;
bvol_lap=lap_conv3d(bvol)>0.5;

D = (2 *bwdist(bvol_lap)).*bvol;
map_msk=lap_conv3d(D)>2.0;
D = D.*map_msk.*(D>1.0) ;
% D = bwdist(bvol).*bvol;
% Find the maximum distance along each axis
% p1=find_max_coordinate(D);
% p2=find_closest_point_to_volume(p1, bvol_lap);
% 
% % D=find_max_projection_and_set_to_zero(D,p1-p2);
% d=max_intensity_projection(D, p1-p2);

max_dist = max(max(max(D)));

% Calculate the thickness of the plate
thickness = max_dist;
end

%%

function [proj, proj_coords] = max_intensity_projection(volume, normal)
    % Find the two axes that lie in the plane
    [~, axis1] = max(abs(normal));
    axis2 = mod(axis1, 3) + 1;

    % Get the dimensions of the projection
    proj_dims = size(volume);
    proj_dims(axis1) = [];
    proj_dims(axis2) = [];

    % Project the volume onto the plane
    proj = zeros(proj_dims);
    [proj_coords1, proj_coords2] = ndgrid(1:proj_dims(1), 1:proj_dims(2));
    proj_coords = zeros([numel(proj_coords1), 3]);
    proj_coords(:, axis1) = proj_coords1(:);
    proj_coords(:, axis2) = proj_coords2(:);
    proj_coords(:, 3) = dot(normal, bsxfun(@minus, proj_coords, [1, 1, 1]), 2) ./ dot(normal, [0,0,1]');

    % Find the maximum intensity projection
    proj_coords = round(proj_coords);
    for i = 1:proj_dims(1)
        for j = 1:proj_dims(2)
            proj(i,j) = max(volume(proj_coords(i,j,1),:,proj_coords(i,j,3)));
        end
    end
end




% function [MIP, MIP_coords] = max_intensity_projection(volume, normal)
% % volume: 3D volume
% % normal: normal vector of the 2D plane
% normal=normal./norm(normal);
% % Find the two axes that lie in the plane
% [~, axis1] = max(abs(normal));
% axis2 = mod(axis1, 3) + 1;
% 
% % Calculate the coordinates of the projection
% [x, y] = meshgrid(1:size(volume, axis2), 1:size(volume, axis1));
% proj_coords = zeros([numel(x), 3]);
% proj_coords(:, axis1) = reshape(y, [], 1);
% proj_coords(:, axis2) = reshape(x, [], 1);
% proj_coords(:, 3) = round((normal* proj_coords' - dot(normal, [1,1,1]')) ./ dot(normal, [0,0,1]'));
% 
% % Clip coordinates outside the volume
% valid_indices = all(proj_coords > 0, 2) & all(proj_coords <= repmat(size(volume), [numel(x), 1]), 2);
% proj_coords = proj_coords(valid_indices, :);
% 
% % Compute the maximum intensity projection
% MIP = zeros(size(volume, axis1), size(volume, axis2));
% MIP_coords = zeros(size(MIP, 1), size(MIP, 2), 3);
% for i = 1:size(proj_coords, 1)
%     x = proj_coords(i, axis2);
%     y = proj_coords(i, axis1);
%     z = proj_coords(i, 3);
%     intensity = volume(y, x, z);
%     if intensity > MIP(y, x)
%         MIP(y, x) = intensity;
%         MIP_coords(y, x, :) = [y, x, z];
%     end
% end
% end


function coords = get_voxel_coords(size)
% Helper function to compute the coordinates of all voxels in a volume
[x, y, z] = ndgrid(1:size(1), 1:size(2), 1:size(3));
coords = [x(:), y(:), z(:)];
end

function R = get_rot(u, v)
    % Calculate the cross product of u and v
    w = cross(u, v);
    
    % Calculate the angle between u and v
    theta = acos(dot(u, v));
    
    % Check if the vectors are parallel and handle the special case
    if theta < eps
        R = eye(3);
    else
        % Calculate the rotation matrix using the axis-angle representation
        w_skew = [0, -w(3), w(2); w(3), 0, -w(1); -w(2), w(1), 0];
        R = eye(3) + w_skew*sin(theta) + w_skew^2*(1-cos(theta));
    end
end

function max_proj_volume = find_max_projection_and_set_to_zero(volume, direction)
    
    direction = direction / norm(direction);
    
    % Compute the maximum intensity projection of the volume in the direction vector
    projection = sum(volume .* repmat(direction', size(volume)), 4);
    [max_projection, idx] = max(projection, [], 'all', 'linear');
    
    % Set all voxels in the volume that are not in the projection plane to 0
    [xx, yy, zz] = ind2sub(size(volume), idx);
    [~, ~, V] = svd(normal);
    x_range = 1:size(volume, 1);
    y_range = 1:size(volume, 2);
    z_range = 1:size(volume, 3);
    x_range(xx) = [];
    y_range(yy) = [];
    z_range(zz) = [];
    [X, Y, Z] = ndgrid(x_range, y_range, z_range);
    coords = [X(:), Y(:), Z(:)];
    projection_plane = [xx*ones(size(coords, 1), 1), yy*ones(size(coords, 1), 1), zz*ones(size(coords, 1), 1)];
    coords_on_plane = coords - projection_plane;
    coords_on_plane_rot = (V' * coords_on_plane')';
    max_dist = max(abs(coords_on_plane_rot), [], 2);
    mask = max_dist < 1e-3;  % tolerance for numerical errors
    volume(~mask) = 0;
    
    % Display the maximum intensity projection
    figure;
    imshow(max_projection, []);
    
    
    
    
    
    % volume is the 3D matrix representing the volume
    % direction is a 3-element vector specifying the direction
    direction = direction./sqrt(sum((direction.^2)));
    std_vec=[1,0,0];
    R=get_rot(std_vec,direction);
    
    % % Determine the direction of the parallel view
    % view_direction = [1, 0, 0]; % for example, this is the x-axis direction
    
    % Calculate the maximum 2D projection in the specified direction
    max_proj = max(volume, [], find(direction));
    
    % Create a binary mask of the max projection
    max_proj_mask = max_proj > 0;
    
    % Set all voxels in the volume that are not in the max projection to 0
    volume(~max_proj_mask) = 0;
    
    % Display the max projection
    imshow(max_proj);
    
    
    
    % Compute the dot product between each coordinate and the direction vector
    [X, Y, Z] = ndgrid(1:size(volume,1), 1:size(volume,2), 1:size(volume,3));
    coords = [X(:) Y(:) Z(:)];
    proj_coords = coords * direction(:);
    
    % Reshape the projection coordinates into a 3D matrix
    proj_coords = reshape(proj_coords, size(volume));
    
    % Compute the maximum value projection along the appropriate dimension
    [max_proj, idx] = max(proj_coords, [], 'all', 'linear');
    [max_proj_i, max_proj_j, max_proj_k] = ind2sub(size(volume), idx);
    
    % Create a binary mask of voxels that are part of the maximum projection
    max_proj_mask = zeros(size(volume));
    max_proj_mask(max_proj_i, max_proj_j, max_proj_k) = 1;
    % max_proj_mask(max_proj_j, max_proj_i, max_proj_k) = 1;
    
    % Set all voxels that are not part of the maximum projection to zero
    max_proj_volume = volume .* max_proj_mask;

end

% % Create a random 3D matrix of size 10x10x10
% volume = rand(10, 10, 10);
% 
% % Specify the direction as a unit vector
% direction = [1, 1, 1] / sqrt(3);
% 
% % Find the maximum value projection along the direction vector and set all other voxels to zero
% max_proj_volume = find_max_projection_and_set_to_zero(volume, direction);
% 
% % Display the resulting volume using the maximum value as the colormap limit
% max_val = max(max_proj_volume(:));
% figure;
% sliceViewer(max_proj_volume, [], [], [], [], [], [], max_val);
% This code creates a random 3D matrix of size 10x10x10, specifies the direction as a unit vector, and computes the maximum value projection in the volume along that direction using the find_max_projection_and_set_to_zero function. The resulting volume is then displayed using the sliceViewer function from the Image Processing Toolbox, with the maximum value as the colormap limit.






function max_coord = find_max_coordinate(volume)
% volume is the 3D matrix representing the volume

% Find the maximum value and its index in the volume
[max_value, max_idx] = max(volume(:));

% Convert the 1D index to 3D coordinates
[x, y, z] = ind2sub(size(volume), max_idx);
max_coord = [x, y, z];
end

function [closest_point, distance] = find_closest_point_to_volume(point, volume, voxel_size)
if nargin<3
    voxel_size=[1,1,1];
end
% point is a 3-element vector representing the query point
% volume is a 3D matrix representing the volume
% voxel_size is a 3-element vector representing the size of a voxel in each dimension

% Compute the coordinates of all voxels in the volume
[x, y, z] = ndgrid(1:size(volume, 1), 1:size(volume, 2), 1:size(volume, 3));
voxel_coords = [x(:), y(:), z(:)];

% Compute the physical coordinates of each voxel
phys_coords = bsxfun(@times, voxel_coords - 1, voxel_size(:)') + voxel_size(:)'/2;

% Compute the distances between the query point and each voxel
distances = sqrt(sum((phys_coords - point(:)').^2, 2));

% Find the index of the closest voxel
[~, idx] = min(distances);

% Convert the index to voxel coordinates
[x, y, z] = ind2sub(size(volume), idx);

% Compute the physical coordinates of the closest voxel
closest_voxel = bsxfun(@times, [x, y, z] - 1, voxel_size(:)') + voxel_size(:)'/2;

% Compute the distance between the query point and the closest voxel
distance = sqrt(sum((closest_voxel - point(:)').^2));

% Return the physical coordinates of the closest point
closest_point = closest_voxel;

end


% function closest_point = find_closest_point_to_volume(point, volume)
% % point is the 3D coordinate of the specific point
% % volume is a 3D matrix representing the volume
% 
% % Compute the distance from the point to each point in the volume
% distances = sqrt(sum((volume - repmat(point, size(volume(:), 1), 1)).^2, 2));
% 
% % Find the point with the smallest distance to the point
% [~, idx] = min(distances);
% 
% % Convert the 1D index to 3D coordinates
% [x, y, z] = ind2sub(size(volume), idx);
% closest_point = [x, y, z];
% 
% % Check if the closest point is inside the volume
% if volume(closest_point(1), closest_point(2), closest_point(3)) == 1
%     % The point is inside the volume, return the point
%     return
% end
% 
% % Find the boundary of the volume
% boundary = zeros(size(volume));
% boundary(1,:,:) = 1;
% boundary(end,:,:) = 1;
% boundary(:,1,:) = 1;
% boundary(:,end,:) = 1;
% boundary(:,:,1) = 1;
% boundary(:,:,end) = 1;
% 
% % Project the closest point onto the boundary
% distances = sqrt(sum((boundary - repmat(closest_point, size(boundary(:), 1), 1)).^2, 2));
% [~, idx] = min(distances);
% [x, y, z] = ind2sub(size(volume), idx);
% closest_point = [x, y, z];
% end



% % Create a 3D matrix with a thin plate of arbitrary orientation
% matrix = zeros(10, 10, 10);
% [x,y,z] = meshgrid(1:10);
% matrix(sqrt((x-5).^2 + (y-5).^2 + (z-5).^2) <= 3) = 1;
% 
% % Calculate the thickness of the plate
% thickness = plate_thickness(matrix);
% disp(thickness);




% 
% 
% 
% function asd=avg_surf_dist(vol1,vol2,side_idx,sample_rate)
% if nargin<4
%     sample_rate=1;
% if nargin<3
%     side_idx=1;
% end
% end
% % % lap_k1=[];
% % lap_k0=[0,0,0;
% %         0,-1,0;
% %         0,0,0];
% % lap_k2=[0,-1,0;
% %         -1,6,-1;
% %         0,-1,0];
% % lap_kernel=cat(3,lap_k0,lap_k2,lap_k0);
% % v1=convn(double(vol1),lap_kernel,'same')>0.5;
% % v2=convn(double(vol2),lap_kernel,'same')>0.5;
% v1=lap_conv3d(vol1)>0.5;
% v2=lap_conv3d(vol2)>0.5;
% % volshow(v2);
% [~,coord1]=img2ind(v1);
% [~,coord2]=img2ind(v2);
% dist_mat=pdist2(coord1(1:sample_rate:end,:),coord2(1:sample_rate:end,:),"euclidean");
% asd=mean(min(dist_mat,[],side_idx));
% 
% end
% 
% function [ind,pnd]=img2ind(img,thresh)
% if nargin<2
%     thresh=0.000000001;
% end
% img_sz=size(img);
% ind=find(img>thresh);
% pnd=[];
% id=ind;
% for s=img_sz
%     coord=mod(id,s);
%     pnd=[pnd,coord];
%     id=floor(id/s);
% end
% end