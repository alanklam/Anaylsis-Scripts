Nsample=7;
%% 1. convert data in txt to pdb format (note: 100nm = 1A in pdb unit)
for i=1:Nsample
   infile=[num2str(i),'/corrected-nmdar-',num2str(i),'.txt'];
   outfile=[num2str(i),'/nmdar.pdb'];
   txt2pdb(infile,outfile);
   infile=[num2str(i),'/corrected-ampar-',num2str(i),'.txt'];
   outfile=[num2str(i),'/ampar.pdb'];
   txt2pdb(infile,outfile);
   infile=[num2str(i),'/corrected-homer-',num2str(i),'.txt'];
   outfile=[num2str(i),'/homer.pdb'];
   txt2pdb(infile,outfile);
end
%% 2. exponential fit and filter data using self-correlation length dx
namelist=['nmdar';'ampar';'homer'];
% for i=1:Nsample
%     for j=1:3;
for i=7:7
    for j=1:1;
        infile=[num2str(i),'/gor_',namelist(j,:),'.dat'];
        gor=importdata(infile);
        A0=max(gor(:,2));
        r0=max(gor(:,1));
        fo = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0],...
            'Upper',[A0*2,r0],...
            'StartPoint',[A0 r0/10]);
        ft = fittype('A*exp(-x/r)','options',fo);
        range=find(gor(:,2)>0);
        fitted=fit(gor(range,1),(gor(range,2)),ft);
        plot(fitted,gor(:,1),(gor(:,2)),'x'); xlim([0,r0]); ylim([0,A0]);
        set(gca,'FontSize',25,'LineWidth',3)
        xlabel('r (100 nm)'); ylabel('g(r)');
        title('Exponential fit')
        cv=coeffvalues(fitted);
        disp(namelist(j,:))
        dx=cv(2);
%         tic
%        grid_filter([num2str(i),'/',namelist(j,:),'.pdb'],dx);
%         toc
    end
    disp(['finished cell',num2str(i)])
end
%% 3. Plot Inter-cluster distance statistics
bin=10; %set bin size, 10=100nm
dn=importdata('output/dist_NH_all_samples.dat'); da=importdata('output/dist_AH_all_samples.dat'); dan=importdata('output/dist_AN_all_samples.dat'); 
x=0:bin:80;
figure(1); y1=hist(dn(:,1),x); y2=hist(da(:,1),x); y3=hist(dan(:,1),x);
xi=0:1:80; 
y1i=interp1(x,y1/sum(y1),xi,'spline'); y1i(y1i<0)=0;
y2i=interp1(x,y2/sum(y2),xi,'spline'); y2i(y2i<0)=0;
y3i=interp1(x,y3/sum(y3),xi,'spline'); y3i(y3i<0)=0;
plot(xi/10,y1i,'r',xi/10,y2i,'k',xi/10,y3i,'b',x/10,y1/sum(y1),'rx',x/10,y2/sum(y2),'kx',x/10,y3/sum(y3),'bx','LineWidth',2,'MarkerSize',15)
grid on;
% plot(x/10,y1/sum(y1),'rx-',x/10,y2/sum(y2),'kx-',x/10,y3/sum(y3),'bx-','LineWidth',2,'MarkerSize',15)
set(gca,'FontSize',25,'LineWidth',3)
xlabel('r (100 nm)'); ylabel('Normalized count');
legend('NMDAR-HOMER','AMPAR-HOMER','NMDAR-AMPAR')
title(['Inter-cluster distance,',num2str(Nsample),' samples'])

%    Rg seems to have no correlation with the inter-cluster distance, shown
%    in the following plots
% figure; plot(dn(:,1)/10,dn(:,2),'xr','MarkerSize',15); ylim([0,3.5])
% set(gca,'FontSize',25,'LineWidth',3)
% title(['NMDAR Rg against NMDAR-HOMER distance']); xlabel('r (100 nm)'); ylabel('Rg');
% figure; plot(da(:,1)/10,da(:,2),'xr','MarkerSize',15); ylim([0,3.5])
% set(gca,'FontSize',25,'LineWidth',3)
% title(['AMPAR Rg against AMPAR-HOMER distance']); xlabel('r (100 nm)'); ylabel('Rg');
% figure; plot(dan(:,1)/10,dan(:,2),'xr','MarkerSize',15); ylim([0,3.5])
% set(gca,'FontSize',25,'LineWidth',3)
% title(['AMPAR Rg against NMDAR-AMPAR distance']); xlabel('r (100 nm)'); ylabel('Rg');
%% 3.1 Plot Inter-cluster distance statistics for different sample
bin=8; %set bin size, 10=100nm
for i=1:Nsample
    dn=importdata([num2str(i),'/dist_NH.dat']); da=importdata([num2str(i),'/dist_AH.dat']); dan=importdata([num2str(i),'/dist_AN.dat']);
    x=0:bin:80;
    figure(1); y1=hist(dn(:,1),x); y2=hist(da(:,1),x); y3=hist(dan(:,1),x);
    xi=0:1:80;
    y1i=interp1(x,y1/sum(y1),xi,'spline'); y1i(y1i<0)=0;
    y2i=interp1(x,y2/sum(y2),xi,'spline'); y2i(y2i<0)=0;
    y3i=interp1(x,y3/sum(y3),xi,'spline'); y3i(y3i<0)=0;
    figure(i); plot(xi/10,y1i,'r',xi/10,y2i,'k',xi/10,y3i,'b',x/10,y1/sum(y1),'rx',x/10,y2/sum(y2),'kx',x/10,y3/sum(y3),'bx','LineWidth',2,'MarkerSize',15)
    grid on;
    % plot(x/10,y1/sum(y1),'rx-',x/10,y2/sum(y2),'kx-',x/10,y3/sum(y3),'bx-','LineWidth',2,'MarkerSize',15)
    set(gca,'FontSize',25,'LineWidth',3)
    xlabel('r (100 nm)'); ylabel('Normalized count');
    legend('NMDAR-HOMER','AMPAR-HOMER','NMDAR-AMPAR')
    title(['Inter-cluster distance, cell ',num2str(i)])
end
%% 4. Plot Rg statistics
bin=0.3; %set bin size, 1=100nm
NmdarR=[]; AmparR=[]; HomerR=[];
for i=1:Nsample
    NmdarR=[NmdarR,importdata([num2str(i),'/Nmdar_R.dat'])]; AmparR=[AmparR,importdata([num2str(i),'/Ampar_R.dat'])]; HomerR=[HomerR,importdata([num2str(i),'/Homer_R.dat'])];
end
x=0:bin:3;
figure(1); b_n=hist(NmdarR,x);  std(NmdarR)
b_a=hist(AmparR,x);  std(AmparR)
b_h=hist(HomerR,x);  std(HomerR)

xi=0:0.01:3;
yn=interp1(x,b_n/sum(b_n),xi,'spline'); yn(yn<0)=0;
ya=interp1(x,b_a/sum(b_a),xi,'spline'); ya(ya<0)=0;
yh=interp1(x,b_h/sum(b_h),xi,'spline'); yh(yh<0)=0;
plot(xi,yn,'r',xi,ya,'k',xi,yh,'b--',x,b_n/sum(b_n),'rx',x,b_a/sum(b_a),'kx',x,b_h/sum(b_h),'bx','LineWidth',2,'MarkerSize',15)
xlim([0,max(x)]); grid on;
legend('NMDAR','AMPAR','HOMER')
set(gca,'FontSize',25,'LineWidth',3)
xlabel('Rg (100 nm)'); ylabel('Normalized count');
title(['Radius of gyration,',num2str(Nsample),' samples'])
%% 5. Plot statistics of angle between N-H-A pair
bin=20; %set bin size, in degree
angle=[];
for s=1:Nsample
 angle=[angle;importdata([num2str(s),'/angle.dat'])]; 
end
angle=angle/pi*180;
x=0:bin:180;
angle_h=hist(angle,x);  std(angle)
xi=0:1:180;
yi=interp1(x,angle_h/sum(angle_h),xi,'spline'); yi(yi<0)=0;
plot(xi,yi,x,angle_h/sum(angle_h),'bx','LineWidth',2,'MarkerSize',15)
set(gca,'FontSize',25,'LineWidth',3)
xlabel('angle (degree)'); ylabel('Normalized count'); xlim([0,180]);
title(['Angle between N-H-A,',num2str(Nsample),' samples'])
%% 6. Plot cross-correlation between protein position, with exponential fit
namelist=['NA';'NH';'AH']; 
for i=1:Nsample
    figure;
    for j=1:3;
        infile=['output/gor_',namelist(j,:),'_sample',num2str(i),'.dat'];
        %infile=['output/gor_',namelist(j,:),'_sample',num2str(i),'_raw.dat'];
        gor=importdata(infile);
        A0=max(gor(:,2));
        r0=max(gor(:,1));
        %single exp fit
        fo = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0],...
            'Upper',[A0*2,r0],...
            'StartPoint',[A0 r0/10]);
        ft = fittype('A*exp(-x/r)','options',fo);
        
        range=find(gor(:,1)>0.7); %fit only data for r>70 nm seems to give a better fit
        fitted=fit(gor(range,1),(gor(range,2)),ft);
        subplot(1,3,j); plot(fitted,gor(:,1),(gor(:,2)),'x'); xlim([0,r0]); ylim([0,A0]);
        set(gca,'FontSize',25,'LineWidth',3)
        xlabel('r (100 nm)'); ylabel('g(r)ij');
        title(['g(r) ',namelist(j,:),' sample ',num2str(i)])
        cv=coeffvalues(fitted);
        
        A=cv(1); dx=cv(2);
        txt = ['A=',num2str(A),' dx=',num2str(dx)];
        text(2,A0/2,txt); 
        set(gca,'FontSize',25,'LineWidth',3)
    end
end
%% 7 poisson point process simulation
X = zeros(10^5,2); %initialise the points
kappa = 20; alpha = 50; r = 0.1; %parameters 
meanpts=kappa*(1 + 2*r)^2;
N = poissrnd(meanpts); %number of cluster centers 
C = rand(N,2)*(1+2*r) - r; %draw cluster centers 
total_so_far = 0;
for c=1:N
    NC = poissrnd(alpha); %number of points in cluster
    k = 0;
    while k < NC %draw uniformly in the n-ball via accept-reject
        Y = 2*r*rand(1,2) - r; %candidate point
        if norm(Y) < r
            X(total_so_far+k+1,:) = C(c,:) + Y;
            k = k+1;
        end
    end
    total_so_far = total_so_far + NC;
end
X = X(1:total_so_far,:); %cut off unused rows
plot(X(:,1),X(:,2),'x') 
axis([0, 1,0, 1])
%% correlated point process 2D, unconstrianted
xmax=800; ymax=800; %zmax=10;
Ngen=10000; Nnow=1000; xi=4;

pt=zeros(Ngen,2);
pt(1:Nnow,:)=[rand(Nnow,1)*xmax,rand(Nnow,1)*ymax];
while Nnow<Ngen
    id=randi([1,Nnow]);
    angle=rand*2*3.14159265;
    d=-xi*log(1-rand);
    xnew=pt(id,1)+d*sin(angle); ynew=pt(id,2)+d*cos(angle);
    if xnew>=0 && xnew<=xmax && ynew>=0 && ynew<=ymax
        Nnow=Nnow+1;
        pt(Nnow,:)=[xnew,ynew];
    end
end
scatter(pt(:,1),pt(:,2),'x')
%% correlated point process 3D, one after one
xmax=800; ymax=800; zmax=200;
Ngen=1000; Nnow=10; xi=3;
pt=zeros(Ngen,3); 
pt(1:Nnow,:)=[rand(Nnow,1)*xmax,rand(Nnow,1)*ymax,rand(Nnow,1)*zmax];

while Nnow<Ngen
    id=randi([1,Nnow]);
    phi=rand*2*3.14159265;
    theta=rand*3.14159265;
    d=0;
    while d<0.5
     d=-xi*log(1-rand);
    end
    
    xnew=pt(id,1)+d*sin(theta)*cos(phi); ynew=pt(id,2)+d*sin(theta)*sin(phi);
    znew=pt(id,3)+d*cos(theta);
    if xnew>=0 && xnew<=xmax && ynew>=0 && ynew<=ymax && znew>=0 && znew<=zmax
        Nnow=Nnow+1;
        pt(Nnow,:)=[xnew,ynew,znew];
    end
    
end

scatter3(pt(:,1),pt(:,2),pt(:,3),'x')
%% correlated point process 3D
xmax=800; ymax=800; zmax=50;
Ngen=10000; Ncen=5000; Nnow=Ncen; xi=1.5;
center=[rand(Ncen,1)*xmax,rand(Ncen,1)*ymax,rand(Ncen,1)*zmax];

pt=zeros(Ngen,3); 

pt(1:Ncen,:)=center;

for id=1:Ncen;
    Ncount=0;
    while Ncount<(Ngen/Ncen-1)
        phi=rand*2*3.14159265;
        theta=rand*3.14159265;
        d=-xi*log(1-rand)+0.2;
        
        xnew=pt(id,1)+d*sin(theta)*cos(phi); ynew=pt(id,2)+d*sin(theta)*sin(phi); 
        znew=pt(id,3)+d*cos(theta);
        if xnew>=0 && xnew<=xmax && ynew>=0 && ynew<=ymax && znew>=0 && znew<=zmax
            Nnow=Nnow+1;
            Ncount=Ncount+1;
            pt(Nnow,:)=[xnew,ynew,znew];
        end
        
    end
end
scatter3(pt(:,1),pt(:,2),pt(:,3),'x')
%% correlated point process 3D, constrianted to sphere surface
xmax=800; ymax=800; 
Ngen=50; Ncen=1; Nnow=Ncen; xi=1.5;
center=[rand(Ncen,1)*xmax,rand(Ncen,1)*ymax,repmat(5,Ncen,1)];
radius=5;

pt=zeros(Ngen,3); ptcen=zeros(Ngen,1);
cosz=rand(Ncen,1)*2-1; sinz=sqrt(1-cosz.^2);
phi=rand(Ncen,1)*2*3.14159265;
pt(1:Ncen,:)=center+radius*[sinz.*cos(phi),sinz.*sin(phi),cosz];
ptcen(1:Ncen)=(1:1:Ncen)';
t=1;
while Nnow<Ngen
    id=randi([1,Nnow]);
    angle=rand*2*3.14159265;
    d=0;
    while d==0 || d>=2*radius
       d=-xi*log(1-rand);
    end
    x0=pt(id,:)-center(ptcen(id),:); 
    if x0(2)>=0
        anglez=acos(x0(1)/norm(x0)); 
    else
        anglez=-acos(x0(1)/norm(x0));
    end
    if x0(1)>=0
        angley=-asin(x0(3)/norm(x0));
    elseif x0(3)>=0
        angley=-(3.14159265-asin(x0(3)/norm(x0)));
    elseif x0(3)<0 
        angley=3.14159265+asin(x0(3)/norm(x0));
    end
    Rz=[cos(anglez),-sin(anglez),0;sin(anglez),cos(anglez),0;0,0,1];
    Ry=[cos(angley),0,sin(angley);0,1,0;-sin(angley),0,cos(angley)];
    xnew=[radius;0;0]+[-d^2/2/radius;d*sqrt(1-d^2/4/radius^2)*sin(angle);d*sqrt(1-d^2/4/radius^2)*cos(angle)];
    xnew=(Rz*Ry*xnew)'+center(ptcen(id),:); %tranform back to original orientation
    
    Nnow=Nnow+1;
    pt(Nnow,:)=xnew;
    ptcen(Nnow)=ptcen(id);
end
scatter3(pt(:,1),pt(:,2),pt(:,3),'x')
%% homogeneous point process 3D, constrianted to sphere surface
xmax=800; ymax=800; 
Ngen=50; Ncen=1; Nnow=0; 
radius=5;
center=[rand(Ncen,1)*xmax,rand(Ncen,1)*ymax,repmat(radius,Ncen,1)];

pt=zeros(Ngen,3);
while Nnow<Ngen
    id=randi([1,Ncen]);
    %theta=rand*3.14159265; 
    cosz=rand*2-1; sinz=sqrt(1-cosz^2);
    phi=rand*2*3.14159265;
    %radius2=rand*radius;
    %xnew=center(id,:)+radius*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
    xnew=center(id,:)+radius*[sinz*cos(phi),sinz*sin(phi),cosz];
    Nnow=Nnow+1;
    pt(Nnow,:)=xnew;
end
scatter3(pt(:,1),pt(:,2),pt(:,3),'x')
%% homogeneous point process 3D
xmax=800; ymax=800; zmax=10;
%Ngen=round((1+sqrt(1+8*xmax*ymax*zmax/1.5^3))/2); 
Ngen=20000;
pt=[rand(Ngen,1)*xmax,rand(Ngen,1)*ymax,rand(Ngen,1)*zmax];
figure; scatter3(pt(:,1),pt(:,2),pt(:,3),'x')
%% homogeneous point process 2D
xmax=800; ymax=800; 
Ngen=10000; 
pt=[rand(Ngen,1)*xmax,rand(Ngen,1)*ymax];
scatter(pt(:,1),pt(:,2),'x')
%% compute g(r) for 2D or 3D simulations
% gor3D( sel1, sel2, rmax, dr ,V): sel2=[] for self-correlation, V=[] to
% use data dimension
%V=3.14*radius^3*4/3*100;

gor3=gor3D(pt,[],3,0.01,[],[xmax,0;ymax,0;zmax0]);
figure; scatter(gor3(:,1),gor3(:,2)); 

%% Expoential fit
gor=gor3;
A0=max(gor(:,2));
r0=max(gor(:,1));
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0],...
    'Upper',[A0*2,r0],...
    'StartPoint',[A0 r0/10]);
ft = fittype('A*exp(-x/r)+1','options',fo);
range=find(gor(:,2)>0); 
fitted=fit(gor(10:end,1),(gor(10:end,2)),ft);
figure; plot(fitted,gor(:,1),(gor(:,2)),'x'); xlim([0,r0]); ylim([0,20]);
hold on; plot([0,r0],[1,1],'k');
set(gca,'FontSize',25,'LineWidth',5)
xlabel('r (100 nm)'); ylabel('g(r)');
title('Exponential fit')
cv=coeffvalues(fitted);
A=cv(1); dx=cv(2);
txt = ['A=',num2str(A),' dx=',num2str(dx)];
text(2,20/2,txt);
set(gca,'FontSize',25,'LineWidth',3)