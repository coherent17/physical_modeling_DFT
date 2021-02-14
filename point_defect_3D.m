clc;
clear;

%n=1
% eigenvector1=[-0.18861,-0.00000];
% eigenvector2=[0.38850, 0.00001];
% eigenvector3=[0.44208, 0.00000];
% eigenvector4=[0.49666,-0.00000];
% eigenvector5=[-0.18569, 0.32161];
% eigenvector6=[-0.22599, 0.00000];
% eigenvector7=[0.42710, 0.00000];

%n=2
% eigenvector1=[-0.20577,-0.35516];
% eigenvector2=[0.42646,-0.00059];
% eigenvector3=[0.18567,-0.32160];
% eigenvector4=[0.01484, 0.02573];
% eigenvector5=[-0.22102,-0.38286];
% eigenvector6=[0.20210, 0.35101];
% eigenvector7=[-0.38909,-0.00050];

%n=3
% eigenvector1=[0.09336,-0.16281];
% eigenvector2=[0.38896,-0.00028];
% eigenvector3=[-0.22104,-0.38284];
% eigenvector4=[-0.24828, 0.43010];
% eigenvector5=[0.37135, 0.00001];
% eigenvector6=[0.11392,-0.19625];
% eigenvector7=[0.42665, 0.00027];

%n=4
eigenvector1=[-0.20576, 0.35515];
eigenvector2=[0.42646, 0.00060];
eigenvector3=[0.18567, 0.32160];
eigenvector4=[0.01484,-0.02573];
eigenvector5=[0.44208, 0.00000];
eigenvector6=[0.20210,-0.35101];
eigenvector7=[-0.38909, 0.00050];

%n=5
% eigenvector1=[0.09336, 0.16280];
% eigenvector2=[0.38897, 0.00028];
% eigenvector3=[-0.22104, 0.38284];
% eigenvector4=[-0.24828,-0.43010];
% eigenvector5=[-0.18567,-0.32160];
% eigenvector6=[0.11392, 0.19625];
% eigenvector7=[0.42665,-0.00027];

%n=6
% eigenvector1=[0.40936,-0.00000];
% eigenvector2=[0.42748, 0.00000];
% eigenvector3=[-0.37137,-0.00000];
% eigenvector4=[-0.02973,-0.00000];
% eigenvector5=[-0.22104, 0.38284];
% eigenvector6=[-0.40586,-0.00000];
% eigenvector7=[-0.38825,-0.00000];

%n=7
% eigenvector1=[0.62406, 0.00000];
% eigenvector2=[0.00102,-0.00001];
% eigenvector3=[-0.00001,-0.00000];
% eigenvector4=[0.50737, 0.00000];
% eigenvector5=[-0.00000, 0.00000];
% eigenvector6=[0.59424, 0.00000];
% eigenvector7=[-0.00091, 0.00000];

data1=importdata("k=( 0.849, 0.000, 0.000)_1BS.dat");
data2=importdata("k=( 0.849, 0.000, 0.000)_2BS.dat");
data3=importdata("k=( 0.849, 0.000, 0.000)_3BS.dat");
data4=importdata("k=( 0.849, 0.000, 0.000)_4BS.dat");
data5=importdata("k=( 0.849, 0.000, 0.000)_5BS.dat");
data6=importdata("k=( 0.849, 0.000, 0.000)_6BS.dat");
data7=importdata("k=( 0.849, 0.000, 0.000)_7BS.dat");

isovalue_AbsSqr=0.45*0.0047;
tolerance_AbsSqr=0.00175;
isovalue_AbsReIm=0.4*0.0685;
tolerance_AbsReIm=0.013;

%n=4 part
flag1=0;
for i=1:+1:length(data1.data(:,1))
    if data1.data(i,1)<1.9 && data1.data(i,1)>-1.85 && data1.data(i,2)>1.68 && data1.data(i,2)<7.3 && data1.data(i,3)<12.5 && data1.data(i,3)>9
        flag1=flag1+1;
        psi1(1,flag1)=(eigenvector1(1)+eigenvector1(2)*1i)*(data1.data(i,4)+data1.data(i,5)*1i)+ ...
                      (eigenvector2(1)+eigenvector2(2)*1i)*(data2.data(i,4)+data2.data(i,5)*1i)+ ...
                      (eigenvector3(1)+eigenvector3(2)*1i)*(data3.data(i,4)+data3.data(i,5)*1i)+ ...
                      (eigenvector4(1)+eigenvector4(2)*1i)*(data4.data(i,4)+data4.data(i,5)*1i)+ ...
                      (eigenvector5(1)+eigenvector5(2)*1i)*(data5.data(i,4)+data5.data(i,5)*1i)+ ...
                      (eigenvector6(1)+eigenvector6(2)*1i)*(data6.data(i,4)+data6.data(i,5)*1i)+ ...
                      (eigenvector7(1)+eigenvector7(2)*1i)*(data7.data(i,4)+data7.data(i,5)*1i);
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
title('(n=4)Re(Ψ) at K point');
axis([-4 4 0 8])
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
title('(n=4)Im(Ψ) at K point');
axis([-4 4 0 8])
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
title('(n=4)|Ψ|^2 at K point');
axis([-4 4 0 8])
grid on

function plot_carbon_atoms
    a0=2.4669;
    a1=a0*[0.5,sqrt(3)/2,0];
    a2=a0*[-0.5,sqrt(3)/2,0];

%the carbon atoms in the unit cell

    %point defect:
    X1(1)=(a1(1)+a2(1))/3;
    Y1(1)=(a1(2)+a2(2))/3;
    X1(2)=(a1(1)+a2(1))/3*2;
    Y1(2)=(a1(2)+a2(2))/3*2;
    Z1(1)=10.25;
    Z1(2)=10.25;
    
    X2(1)=a1(1)/3+4*a2(1)/3;
    Y2(1)=a1(2)/3+4*a2(2)/3;
    X2(2)=2*a1(1)/3+5*a2(1)/3;
    Y2(2)=2*a1(2)/3+5*a2(2)/3;
    Z2(1)=10.25;
    Z2(2)=10.25;
    
    X3(1)=4*a1(1)/3+a2(1)/3;
    Y3(1)=4*a1(2)/3+a2(2)/3;
    X3(2)=5*a1(1)/3+2*a2(1)/3;
    Y3(2)=5*a1(2)/3+2*a2(2)/3;
    Z3(1)=10.25;
    Z3(2)=10.25;
    
    X4(1)=4*a1(1)/3+4*a2(1)/3;
    Y4(1)=4*a1(2)/3+4*a2(2)/3;
    X4(2)=5*a1(1)/3+5*a2(1)/3;
    Y4(2)=5*a1(2)/3+5*a2(2)/3;
    Z4(1)=10.25;
    Z4(2)=10.25;
    
    plot3(X1(2),Y1(2),Z1(2),'k.','Markersize',80);hold on;
    plot3(X2,Y2,Z2,'k.','Markersize',80);hold on;
    plot3(X3,Y3,Z3,'k.','Markersize',80);hold on;
    plot3(X4,Y4,Z4,'k.','Markersize',80);hold on;
end