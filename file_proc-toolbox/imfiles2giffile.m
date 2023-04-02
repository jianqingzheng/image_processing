function []=imfiles2giffile(imdirs,giffile,delay,loopnumb)
if nargin<4
    loopnumb=inf;
    if nargin<3
        delay=0.001;
    end
end
% im_dirs_temp=get_dirs(Label_path,'.obj',[Label_path,'/',sprintf('P%02d',p)]);
for i=1:numel(imdirs)
    img=imread(imdirs{i});
    drawnow
    f=getframe;
    f=frame2im(f);
    [X,map]=rgb2ind(img,256);
    if i==1
        imwrite(X,map,giffile,'gif','LoopCount',loopnumb,'DelayTime',delay);
    else
        imwrite(X,map,giffile,'WriteMode','Append','DelayTime',delay);
    end
end
end