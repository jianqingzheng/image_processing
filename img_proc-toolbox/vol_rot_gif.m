function []=vol_rot_gif(fig_id,gif_file,cam_angles)
if nargin<3
    rot_num=5;
%     cam_angles=[[0:rot_num:360]+90;zeros(1,360/rot_num+1)+5];
    cam_angles=[zeros(1,360/rot_num+1)+90;[0:rot_num:360]+90];
    if nargin<2
        gif_file=['vol_rot_gif','.gif'];
    end
end
if nargin>=1
    fig=figure(fig_id);
end
loopnumb=Inf;
delay=0.1;

for ang_id=1:size(cam_angles,2)
    view(double(cam_angles(:,ang_id))');
    camlight;%
    ax = gca;
%     ax.Units = 'pixels';
%     ti = ax.TightInset;
%     pos=ax.Position;
%     rect = [0, 0, floor(pos(3)), floor(pos(4))];
    axis vis3d
    drawnow
    f=getframe(ax.Parent);
    img=frame2im(f);
    delete(findall(gcf,'Type','light'));%
    [X,map]=rgb2ind(img,256);
    if ang_id==1
        imwrite(X,map,gif_file,'gif','LoopCount',loopnumb,'DelayTime',delay);
    else
        imwrite(X,map,gif_file,'WriteMode','Append','DelayTime',delay);
    end
end
end

