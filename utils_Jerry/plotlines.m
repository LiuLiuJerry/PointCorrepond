 function []=plotlines(xyz, Color)
    a = xyz(:,1)';
    b = xyz(:,2)';
    c = xyz(:,3)';
    Color = [Color; Color(end, :)];
    Color = reshape(Color, 1, size(Color, 1), size(Color, 2));

    %�ռ����ߣ�ɫ��ֵ��zdataΪ׼
    patch('XData',[a, nan],'YData', [b,nan], 'ZData', [c,nan], ... 
                    'CData', Color, ... 
                    'facecolor','none','edgecolor','interp');   %ɫ��ʹ�ò�ֵ����ƽ��
    material dull; light; grid on; xlabel('x'); ylabel('y'); zlabel('z');
    view([60,30]); axis equal; axis manual;
 end
