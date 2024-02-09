clear all;
close all;
clc;

%%
%naca 2412 nomenclature
ymc_c=0.02;
xmc_c=0.4;
tm=0.12;
alpha=10*pi/180;%in radian

%%
%cosine clustring
%choose even number of points on circle 
n=200;
theta=2*pi/(n-1);

%%
%taking projection on x axis
i=1:n/2;
x_c=0.5*(1-cos((i-0.5)*theta));

%%
%camber line formation
for j=1:length(x_c)
    if(x_c(j)>=0) && (x_c(j)<=xmc_c)
        ycx_c(j)=ymc_c*(2*(x_c(j)/xmc_c)-(x_c(j)/xmc_c)^2);
        dyc_dx(j)=ymc_c*(2*(1/xmc_c)-(2*(x_c(j))/(xmc_c)^2));
    else
        ycx_c(j)=ymc_c*(2*((1-x_c(j))/(1-xmc_c))-((1-x_c(j))/(1-xmc_c))^2);
        dyc_dx(j)=ymc_c*(2*(-1/(1-xmc_c))+(2*(1-x_c(j))/(1-xmc_c)^2));
    end
end

%%
%thickness
for j=1:length(x_c)
    tx(j)=tm*(2.969*sqrt(x_c(j))-1.260*x_c(j)-3.516*x_c(j)*x_c(j)+2.843*x_c(j)*x_c(j)*x_c(j)-1.015*x_c(j)*x_c(j)*x_c(j)*x_c(j));
end

%%
%upper and lower surface
for j=1:length(x_c)
    xu(j)=x_c(j)-0.5*tx(j)*dyc_dx(j)/sqrt(1+dyc_dx(j)*dyc_dx(j));
    yu(j)=ycx_c(j)+0.5*tx(j)/sqrt(1+dyc_dx(j)*dyc_dx(j));
    xl(j)=x_c(j)+0.5*tx(j)*dyc_dx(j)/sqrt(1+dyc_dx(j)*dyc_dx(j));
    yl(j)=ycx_c(j)-0.5*tx(j)/sqrt(1+dyc_dx(j)*dyc_dx(j));
end

%%
%normal and axial force
x=x_c(1:end-1);
for j=1:length(x)
    dN(j)=-20000*(x(j)-1)*(x(j)-1)+119000+(288*(yu(j+1)-yu(j))/(x_c(j+1)-x_c(j))+731*(yl(j+1)-yl(j))/(x_c(j+1)-x_c(j)))*x(j)^(-0.2);
    dA(j)=1019*x(j)^(-0.2)+(4*(x(j)-1)*(x(j)-1)+5.4)*10000*(yu(j+1)-yu(j))/(x_c(j+1)-x_c(j))-(2*(x(j)-1)*(x(j)-1)+17.3)*10000*(yl(j+1)-yl(j))/(x_c(j+1)-x_c(j));
end
N = trapz(x,dN,2);
A = trapz(x,dA,2);
b=[N;A];
a=[cos(alpha),sin(alpha);-sin(alpha),cos(alpha)];
F = a\b;
L=F(1);
D=F(2);

%%
%moments
for j=1:length(x)
    dM_le(j)=(20000*(x(j)-1)*(x(j)-1)-119000-(288*(yu(j+1)-yu(j))/(x_c(j+1)-x_c(j))+731*(yl(j+1)-yl(j))/(x_c(j+1)-x_c(j)))*x(j)^(-0.2))*x(j)+((40000*(x(j)-1)^2 + 54000)*(yu(j+1)-yu(j))/(x_c(j+1)-x_c(j))+288*x(j)^(-0.2))*yu(j)+(-20000*(x(j)-1)*(x(j)-1)*(yl(j+1)-yl(j))/(x_c(j+1)-x_c(j))+731*x(j)^(-0.2))*yl(j);
end
M_le=trapz(x,dM_le,2);
M_qc=M_le+N*0.25;

%%
%center of pressure
x_cp=-M_le/N;
