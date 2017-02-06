function [] = grid_filter( infile , dx)
% filter data according to self-correlation length scale

xyz=zeros(1,3);
fid=fopen(infile);
tline=fgets(fid);
c=1;
while ischar(tline)
    if length(tline)>4
        if strcmp(tline(1:4),'ATOM')
            tmp=textscan(tline,'ATOM      %s  CA  UNK X   %d     %f %f %f');
            xyz(c,1)=tmp{3}; xyz(c,2)=tmp{4}; xyz(c,3)=tmp{5};
            c=c+1;
        end
    end
    tline=fgets(fid);
end
% sort the data
tmp=sqrt(xyz(:,1).^2+xyz(:,2).^2+xyz(:,3).^2);
[~,I]=sort(tmp);
xyz=xyz(I,:);
clear I

N=size(xyz,1);
gridbox=struct;
x=ceil(xyz(1,:)/dx);
gridbox.index=x;
gridbox.xyz=xyz(1,:);
for i=2:N
    x=ceil(xyz(i,:)/dx);
    L=size(gridbox,1);
    locate=0;
    for j=L:-1:1
        r=gridbox(j).index-x;
        if r(1)==0 && r(2)==0 && r(3)==0
            gridbox(j).xyz=[gridbox(j).xyz;xyz(i,:)];
            locate=1;
            break
        end
    end
    
    if locate==0
        tmp=struct;
        tmp.index=x;
        tmp.xyz=xyz(i,:);
        gridbox=[gridbox;tmp];
        clear tmp
    end
end

% size(gridbox,1)
% s=0;
% for i=1:size(gridbox,1)
%    s=s+size(gridbox(i).xyz,1); 
% end
% s

xyz_cluster=[];
for i=1:size(gridbox,1)
%     if size(gridbox(i).xyz,1)>cutoff
        xyz_cluster=[xyz_cluster;mean(gridbox(i).xyz,1),size(gridbox(i).xyz,1)];
%     end
end

% writepdb for small_cluster, store number of points in a cluster in beta
data=xyz_cluster;
fid=fopen([infile(1:end-4),'_small.pdb'],'w');
fprintf(fid,'CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n');
for i=1:size(data,1)
    if i>99999
        fprintf(fid,'ATOM%7s  CA  UNK X   0    %8.3f%8.3f%8.3f  1.00%6d\n',lower(dec2hex(i)),data(i,:));
    else
        fprintf(fid,'ATOM%7d  CA  UNK X   0    %8.3f%8.3f%8.3f  1.00%6d\n',i,data(i,:));
    end
end
fprintf(fid,'END\n');
fclose(fid);
clear data

end