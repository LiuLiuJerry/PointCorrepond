for i=0:50
    suf = strcat(num2str(i,'%d'), '.obj');
    name = strcat('D:\GitHub\data\Octopus\octopus\frame2\frame_', suf);
%     [V,F]=readOBJfast(name);
    [V, F] = read_obj(name);
    F = [F(1,:);F(3,:);F(5,:)];
    suf = strcat(num2str(i,'%04d'), '.off');
    name1=strcat('D:\GitHub\data\Octopus\octopus\frame2_off\mesh_',suf);
    writeOFF(name1,V',F');
end