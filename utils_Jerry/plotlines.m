 function []=plotlines(xyz, Color)
    a = xyz(:,1)';
    b = xyz(:,2)';
    c = xyz(:,3)';
    Color = [Color; Color(end, :)];
    Color = reshape(Color, 1, size(Color, 1), size(Color, 2));

    %空间曲线，色彩值以zdata为准
    patch('XData',[a, nan],'YData', [b,nan], 'ZData', [c,nan], ... 
                    'CData', Color, ... 
                    'facecolor','none','edgecolor','interp');   %色彩使用插值策略平滑
    material dull; light; grid on; xlabel('x'); ylabel('y'); zlabel('z');
    view([60,30]); axis equal; axis manual;
 end
