clc
close all
clear


nx=250;
ny=750;

% Base directory (current working directory)
base = pwd;

% Go three levels up: ../../../
data_dir = fullfile(base, '..', '..', '..','..', ...
                    'database', ...
                    'dataRe40Bo40_EUV/data');

data_dir_G = fullfile(base, '..', '..', '..','..', ...
                    'database', ...
                    'dataRe40Bo40_EG/data');


% Read X.csv, Y.csv
x = readmatrix(fullfile(data_dir_G,'x.csv'));
x = x(:);

y = readmatrix(fullfile(data_dir_G,'y.csv'));
y = y(:);

G = readmatrix(fullfile(data_dir_G,'Stacked_1.csv'));
G = G(:);

G=reshape(G,[size(y,1),size(x,1)]);


Velocity = readmatrix(fullfile(data_dir,'Stacked_1.csv'));
Velocity = Velocity(:);

U=reshape(Velocity(1:size(y,1)*size(x,1)),[size(y,1),size(x,1)]);

V=reshape(Velocity(size(y,1)*size(x,1)+1:end),[size(y,1),size(x,1)]);


[XGG,YGG]=meshgrid(x,y);

xmax=max(x(:));
xmin=min(x(:));

ymax=max(y(:));
ymin=min(y(:));

x=linspace(xmin,xmax,nx);
y=linspace(ymin,ymax,ny);

[XG,YG]=meshgrid(x,y);

T=zeros(size(x,2)*size(y,2),5);

T(:,1)=XG(:);
T(:,2)=YG(:);

Gnew=interp2(XGG,YGG,G,XG,YG,"spline");
Unew=interp2(XGG,YGG,U,XG,YG,"spline");
Vnew=interp2(XGG,YGG,V,XG,YG,"spline");

T(:,3)=Unew(:);
T(:,4)=Vnew(:);
T(:,5)=Gnew(:);

G_shaped=reshape(Gnew,[size(y,2),size(x,2)]);
U_shaped=reshape(Unew,[size(y,2),size(x,2)]);
V_shaped=reshape(Vnew,[size(y,2),size(x,2)]);

figure
contourf(x,y,G_shaped,'LineStyle','none')
clim([0 1])
daspect([1,1,1])

figure
contourf(x,y,U_shaped,'LineStyle','none')
daspect([1,1,1])


figure
contourf(x,y,V_shaped,'LineStyle','none')
daspect([1,1,1])