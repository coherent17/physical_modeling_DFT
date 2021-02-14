clc;
clear;

eigenvector1=[0.72292,0];
eigenvector2=[0.68996,0.03650];
eigenvector3=[-0.68996,0.03650];
eigenvector4=[0.72292,0];

data1=importdata("k=( 1.698, 0.000, 0.000)_1BS_n1.dat");
data2=importdata("k=( 1.698, 0.000, 0.000)_2BS_n1.dat");
data3=importdata("k=( 1.698, 0.000, 0.000)_1BS_n2.dat");
data4=importdata("k=( 1.698, 0.000, 0.000)_2BS_n2.dat");

isovalue_AbsSqr=0.45*0.0047;
tolerance_AbsSqr=0.0017;
isovalue_AbsReIm=0.4*0.0685;
tolerance_AbsReIm=0.0060;

%n=1 part
flag1=0;
for i=1:+1:length(data1.data(:,1))
    if data1.data(i,1)<0.36 && data1.data(i,1)>-0.45 && data1.data(i,2)>1.06 && data1.data(i,2)<3.1 && data1.data(i,3)<12.5 && data1.data(i,3)>9
        flag1=flag1+1;
        psi1(1,flag1)=(eigenvector1(1)+eigenvector1(2)*1i)*(data1.data(i,4)+data1.data(i,5)*1i)+(eigenvector2(1)+eigenvector2(2)*1i)*(data2.data(i,4)+data2.data(i,5)*1i);
        x1(1,flag1)=data1.data(i,1);
        y1(1,flag1)=data1.data(i,2);
        z1(1,flag1)=data1.data(i,3);
    else 
        flag1=flag1+0;
    end
end

figure(1)
for i=1:+1:flag1
   if  abs(real(psi1(1,i)))<(isovalue_AbsReIm+tolerance_AbsReIm) && abs(real(psi1(1,i)))>(isovalue_AbsReIm-tolerance_AbsReIm)
       if real(psi1(1,i))>0
           plot3(x1(1,i),y1(1,i),z1(1,i),'r.','MarkerSize',15);hold on;
       elseif real(psi1(1,i))<0
           plot3(x1(1,i),y1(1,i),z1(1,i),'b.','MarkerSize',15);hold on;
       end
   end
end
plot_carbon_atoms
xlabel('x(Å)');
ylabel('y(Å)');
zlabel('z(Å)');
title('(n=1)Re(Ψ) at K point');
axis([-2 2 0 4 8 12])
grid on

figure(2)
for i=1:+1:flag1
   if  abs(imag(psi1(1,i)))<(isovalue_AbsReIm+tolerance_AbsReIm) && abs(imag(psi1(1,i)))>(isovalue_AbsReIm-tolerance_AbsReIm)
       if imag(psi1(1,i))>0
           plot3(x1(1,i),y1(1,i),z1(1,i),'r.','MarkerSize',15);hold on;
       elseif imag(psi1(1,i))<0
           plot3(x1(1,i),y1(1,i),z1(1,i),'b.','MarkerSize',15);hold on;
       end
   end
end
plot_carbon_atoms
xlabel('x(Å)');
ylabel('y(Å)');
zlabel('z(Å)');
title('(n=1)Im(Ψ) at K point');
axis([-2 2 0 4 8 12])
grid on

figure(3)
for i=1:+1:flag1
   if  (sqrt(real(psi1(1,i))^2+imag(psi1(1,i))^2))^2<(isovalue_AbsSqr+tolerance_AbsSqr) && (sqrt(real(psi1(1,i))^2+imag(psi1(1,i))^2))^2>(isovalue_AbsSqr-tolerance_AbsSqr)
       if sqrt(real(psi1(1,i))^2+imag(psi1(1,i))^2)>0
           plot3(x1(1,i),y1(1,i),z1(1,i),'r.','MarkerSize',15);hold on;
       elseif sqrt(real(psi1(1,i))^2+imag(psi1(1,i))^2)<0
           plot3(x1(1,i),y1(1,i),z1(1,i),'b.','MarkerSize',15);hold on;
       end
   end
end
plot_carbon_atoms
xlabel('x(Å)');
ylabel('y(Å)');
zlabel('z(Å)');
title('(n=1)|Ψ|^2 at K point');
axis([-2 2 0 4 8 12])
grid on

%n=2 part
flag2=0;
for i=1:+1:length(data3.data(:,1))
    if data3.data(i,1)<0.62 && data3.data(i,1)>-0.62 && data3.data(i,2)>0.8 && data3.data(i,2)<3.4 && data3.data(i,3)<11.5 && data3.data(i,3)>9
        flag2=flag2+1;
        psi2(1,flag2)=(eigenvector3(1)+eigenvector3(2)*1i)*(data3.data(i,4)+data3.data(i,5)*1i)+(eigenvector4(1)+eigenvector4(2)*1i)*(data4.data(i,4)+data4.data(i,5)*1i);
        x2(1,flag2)=data3.data(i,1);
        y2(1,flag2)=data3.data(i,2);
        z2(1,flag2)=data3.data(i,3);
    else 
        flag2=flag2+0;
    end
end

figure(4)
for i=1:+1:flag2
   if  abs(real(psi2(1,i)))<(isovalue_AbsReIm+tolerance_AbsReIm) && abs(real(psi2(1,i)))>(isovalue_AbsReIm-tolerance_AbsReIm)
       if real(psi2(1,i))>0
           plot3(x2(1,i),y2(1,i),z2(1,i),'r.','MarkerSize',15);hold on;
       elseif real(psi2(1,i))<0
           plot3(x2(1,i),y2(1,i),z2(1,i),'b.','MarkerSize',15);hold on;
       end
   end
end
plot_carbon_atoms
xlabel('x(Å)');
ylabel('y(Å)');
zlabel('z(Å)');
title('(n=2)Re(Ψ) at K point');
axis([-2 2 0 4 8 12])
grid on

figure(5)
for i=1:+1:flag2
   if  abs(imag(psi2(1,i)))<(isovalue_AbsReIm+tolerance_AbsReIm) && abs(imag(psi2(1,i)))>(isovalue_AbsReIm-tolerance_AbsReIm)
       if imag(psi2(1,i))>0
           plot3(x2(1,i),y2(1,i),z2(1,i),'r.','MarkerSize',15);hold on;
       elseif imag(psi2(1,i))<0
           plot3(x2(1,i),y2(1,i),z2(1,i),'b.','MarkerSize',15);hold on;
       end
   end
end
plot_carbon_atoms
xlabel('x(Å)');
ylabel('y(Å)');
zlabel('z(Å)');
title('(n=2)Im(Ψ) at K point');
axis([-2 2 0 4 8 12])
grid on

figure(6)
for i=1:+1:flag2
   if  (sqrt(real(psi2(1,i))^2+imag(psi2(1,i))^2))^2<(isovalue_AbsSqr+tolerance_AbsSqr) && (sqrt(real(psi2(1,i))^2+imag(psi2(1,i))^2))^2>(isovalue_AbsSqr-tolerance_AbsSqr)
       if sqrt(real(psi2(1,i))^2+imag(psi2(1,i))^2)>0
           plot3(x2(1,i),y2(1,i),z2(1,i),'r.','MarkerSize',15);hold on;
       elseif sqrt(real(psi2(1,i))^2+imag(psi2(1,i))^2)<0
           plot3(x2(1,i),y2(1,i),z2(1,i),'b.','MarkerSize',15);hold on;
       end
   end
end
plot_carbon_atoms
xlabel('x(Å)');
ylabel('y(Å)');
zlabel('z(Å)');
title('(n=2)|Ψ|^2 at K point');
axis([-2 2 0 4 8 12])
grid on

function plot_carbon_atoms
    plot3(0,1.4242654,10.25,'k.','MarkerSize',80);hold on;
    plot3(0,2.8485308,10.25,'k.','MarkerSize',80);hold on;
end