clear all;
close all;
clc;

%%
%naca 2412 nomenclature
ymc_c=0.02;
xmc_c=0.4;
tm=0.12;

%%
%equal spacing
n=200;
h=1/(n-1);
x_c=0:h:1;

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
x=[fliplr(xl),xu];
y=[fliplr(yl),yu];
plot(x,y);
xlabel('--- x/c ---->');
ylabel('---- y/c ---->');
title('naca 2412');

