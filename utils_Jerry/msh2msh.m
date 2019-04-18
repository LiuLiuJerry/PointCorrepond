for i=0:249
    suf = strcat(num2str(i,'%04d'), '.obj');
    name = strcat('G:/Result/walk/trans/mesh_', suf);
    [V,F]=readOBJfast(name);
    suf = strcat(num2str(i,'%04d'), '.off');
    name1=strcat('G:/Result/walk/trans_off/mesh_',suf);
    writeOFF(name1,V,F);
end