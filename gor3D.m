%%
% This function calculate the 3D radial distribution of particle
% usage: gor3D(sel1, sel2, rmax, dr, V)
% -sel1 sel2 are N by 3 vectors containing the x,y,z postions of N particles, g of r is calculated between
% the 2 sets of particles, leave sel2=[] for self-radial distribution; if sel1 sel2 are N by 2 vectors,
% 2D gor will be calculated
% -rmax is the maximun radius for calculation
% -dr is the bin size for histogram
% -V is the total volume of the box for normalization, leave V=[] will use the span of particles in space to
% determine the volume


function gr = gor3D( sel1, sel2, rmax, dr ,V)
if size(sel1,2)==3
    
    if isempty(V)
        V=(max(sel1(:,1))-min(sel1(:,1)))*(max(sel1(:,2))-min(sel1(:,2)))*(max(sel1(:,3))-min(sel1(:,3)));
        %V=4/3*3.14159265*rmax^3;
    end
    x=(0:dr:rmax)';
    rup=rmax+dr/2;
    if isempty(sel2)
        N=size(sel1,1); count=1; r=0;
        Np=N*(N-1)/2;
        for i=1:N-1
            for j=i+1:N
                rtmp=sqrt((sel1(i,1)-sel1(j,1))^2+(sel1(i,2)-sel1(j,2))^2+(sel1(i,3)-sel1(j,3))^2);
                if rtmp<rup
                    r(count)=rtmp;
                    count=count+1;
                end
            end
        end
    else
        count=1; r=0;
        Np=size(sel1,1)*size(sel2,1);
        for i=1:size(sel1,1)
            for j=1:size(sel2,1)
                rtmp=sqrt((sel1(i,1)-sel2(j,1))^2+(sel1(i,2)-sel2(j,2))^2+(sel1(i,3)-sel2(j,3))^2);
                if rtmp<rup
                    r(count)=rtmp;
                    count=count+1;
                end
            end
        end
    end
    gr=hist(r,x);
    rho=Np/V;
    %rho=length(r)/V;
    gr(1)=gr(1)/(x(1)+dr/2)^3;
    for i=2:length(gr)
        gr(i)=gr(i)/((x(i)+dr/2)^3-(x(i)-dr/2)^3);
    end
    gr=[x,gr'/rho/4/3.14159265*3];
    
elseif size(sel1,2)==2
    
    if isempty(V)
        V=(max(sel1(:,1))-min(sel1(:,1)))*(max(sel1(:,2))-min(sel1(:,2)));
        %V=3.14159265*rmax^2;
    end
    x=(0:dr:rmax)';
    rup=rmax+dr/2;
    if isempty(sel2)
        N=size(sel1,1);  count=1; r=0;
        Np=N*(N-1)/2;
        for i=1:N-1
            for j=i+1:N
                rtmp=sqrt((sel1(i,1)-sel1(j,1))^2+(sel1(i,2)-sel1(j,2))^2);
                if rtmp<rup
                    r(count)=rtmp;
                    count=count+1;
                end
            end
        end
    else
        Np=size(sel1,1)*size(sel2,1);
        count=1; r=0;
        for i=1:size(sel1,1)
            for j=1:size(sel2,1)
                rtmp=sqrt((sel1(i,1)-sel2(j,1))^2+(sel1(i,2)-sel2(j,2))^2);
                if rtmp<rup
                    r(count)=rtmp;
                    count=count+1;
                end
            end
        end
    end
    gr=hist(r,x);
    rho=Np/V;
    %rho=length(r)/V;
    gr(1)=gr(1)/(x(1)+dr/2)^2;
    for i=2:length(gr)
        gr(i)=gr(i)/((x(i)+dr/2)^2-(x(i)-dr/2)^2);
    end
    gr=[x,gr'/rho/3.14159265];
end

return

