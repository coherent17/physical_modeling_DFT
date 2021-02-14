%cut the surface of the z=10.591667
clc;
clear;

z=10.591667;

%n=1
% eigenvector1=[ 0.40675,-0.01097];
% eigenvector2=[0.27679, 0.06658];
% eigenvector3=[-0.35016,-0.07603];
% eigenvector4=[0.24530,-0.05126];
% eigenvector5=[-0.21638, 0.37472];
% eigenvector6=[-0.38793,-0.12720];
% eigenvector7=[0.13743, 0.05263];
% eigenvector8=[-0.43143,-0.02241];

%n=2
% eigenvector1=[ 0.61184, 0.00000];
% eigenvector2=[0.01321,-0.01577];
% eigenvector3=[0.00330, 0.00496];
% eigenvector4=[0.43288, 0.00000];
% eigenvector5=[0.16706,-0.18675];
% eigenvector6=[0.23439, 0.08745];
% eigenvector7=[0.31740, 0.07764];
% eigenvector8=[0.45351, 0.00000];

%n=3
% eigenvector1=[ 0.39865, 0.01291];
% eigenvector2=[-0.10592,-0.43757];
% eigenvector3=[0.10967, 0.03476];
% eigenvector4=[-0.07810, 0.23798];
% eigenvector5=[-0.21630,-0.37484];
% eigenvector6=[-0.07715, 0.33159];
% eigenvector7=[-0.34332,-0.19770];
% eigenvector8=[-0.30120,-0.10671];

%n=4
eigenvector1=[ 0.08748,-0.16745];
eigenvector2=[0.04215,-0.39644];
eigenvector3=[0.38092, 0.18797];
eigenvector4=[-0.21632, 0.37471];
eigenvector5=[0.07816, 0.23797];
eigenvector6=[0.13073,-0.38546];
eigenvector7=[0.36389, 0.21275];
eigenvector8=[-0.14652,-0.10090];

%n=5
% eigenvector1=[0.41863,-0.00185];
% eigenvector2=[-0.14442, 0.33944];
% eigenvector3=[0.24709, 0.05119];
% eigenvector4=[-0.16689,-0.18672];
% eigenvector5=[0.43279, 0.00000];
% eigenvector6=[-0.00366,-0.37932];
% eigenvector7=[-0.42887,-0.01025];
% eigenvector8=[-0.17432, 0.12905];

%n=6
% eigenvector1=[ 0.10700, 0.17765];
% eigenvector2=[-0.03470, 0.23175];
% eigenvector3=[0.51545,-0.10705];
% eigenvector4=[-0.21632,-0.37472];
% eigenvector5=[-0.24516,-0.05125];
% eigenvector6=[-0.01612, 0.30600];
% eigenvector7=[0.45222, 0.00000];
% eigenvector8=[-0.22936, 0.15488];

%n=7
% eigenvector1=[ -0.01091,-0.00081];
% eigenvector2=[0.27097,-0.08232];
% eigenvector3=[0.54283, 0.00000];
% eigenvector4=[0.24514,-0.05125];
% eigenvector5=[-0.21636, 0.37467];
% eigenvector6=[-0.26822, 0.03972];
% eigenvector7=[-0.36132,-0.08254];
% eigenvector8=[0.39799,-0.07634];

%n=8
% eigenvector1=[ -0.21629,-0.01181];
% eigenvector2=[0.53452, 0.00000];
% eigenvector3=[0.18937,-0.08101];
% eigenvector4=[0.43270,-0.00000];
% eigenvector5=[0.16697,-0.18665];
% eigenvector6=[0.42174, 0.00000];
% eigenvector7=[-0.09360,-0.04777];
% eigenvector8=[-0.41993, 0.09873];

data1=importdata("k=( 0.849, 0.000, 0.000)_1BS.dat");
data2=importdata("k=( 0.849, 0.000, 0.000)_2BS.dat");
data3=importdata("k=( 0.849, 0.000, 0.000)_3BS.dat");
data4=importdata("k=( 0.849, 0.000, 0.000)_4BS.dat");
data5=importdata("k=( 0.849, 0.000, 0.000)_5BS.dat");
data6=importdata("k=( 0.849, 0.000, 0.000)_6BS.dat");
data7=importdata("k=( 0.849, 0.000, 0.000)_7BS.dat");
data8=importdata("k=( 0.849, 0.000, 0.000)_8BS.dat");

%n=4 part
flag1=0;
for i=1:+1:length(data1.data(:,1))
    if data1.data(i,3)==z
        flag1=flag1+1;
        psi1(1,flag1)=(eigenvector1(1)+eigenvector1(2)*1i)*(data1.data(i,4)+data1.data(i,5)*1i)+ ...
                      (eigenvector2(1)+eigenvector2(2)*1i)*(data2.data(i,4)+data2.data(i,5)*1i)+ ...
                      (eigenvector3(1)+eigenvector3(2)*1i)*(data3.data(i,4)+data3.data(i,5)*1i)+ ...
                      (eigenvector4(1)+eigenvector4(2)*1i)*(data4.data(i,4)+data4.data(i,5)*1i)+ ...
                      (eigenvector5(1)+eigenvector5(2)*1i)*(data5.data(i,4)+data5.data(i,5)*1i)+ ...
                      (eigenvector6(1)+eigenvector6(2)*1i)*(data6.data(i,4)+data6.data(i,5)*1i)+ ...
                      (eigenvector7(1)+eigenvector7(2)*1i)*(data7.data(i,4)+data7.data(i,5)*1i)+ ...
                      (eigenvector8(1)+eigenvector8(2)*1i)*(data8.data(i,4)+data8.data(i,5)*1i);
        x1(1,flag1)=data1.data(i,1);
        y1(1,flag1)=data1.data(i,2);
    else 
        flag1=flag1+0;
    end
end

figure(1)
scatter(x1(1,1:flag1),y1(1,1:flag1),5,real(psi1(1,1:flag1)));hold on;
colorbar
title(colorbar,'Re(Ψ)');
plot_carbon_atoms;
title(['(n=4)Re(Ψ) at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -4 12])

figure(2)
scatter(x1(1,1:flag1),y1(1,1:flag1),5,imag(psi1(1,1:flag1)));hold on;
colorbar
title(colorbar,'Im(Ψ)');
plot_carbon_atoms;
title(['(n=4)Im(Ψ) at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -4 12])


figure(3)
scatter(x1(1,1:flag1),y1(1,1:flag1),5,(real(psi1(1,1:flag1)).^2+imag(psi1(1,1:flag1)).^2));hold on;
colorbar
title(colorbar,'|Ψ|^2');
plot_carbon_atoms;
title(['(n=4)|Ψ|^2 at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -4 12])

function plot_carbon_atoms
    a0=2.4669;
    a1=a0*[0.5,sqrt(3)/2,0];
    a2=a0*[-0.5,sqrt(3)/2,0];

%the carbon atoms in the unit cell
    X1(1)=(a1(1)+a2(1))/3;
    Y1(1)=(a1(2)+a2(2))/3;
    X1(2)=(a1(1)+a2(1))/3*2;
    Y1(2)=(a1(2)+a2(2))/3*2;

    X2(1)=a1(1)/3+4*a2(1)/3;
    Y2(1)=a1(2)/3+4*a2(2)/3;
    X2(2)=2*a1(1)/3+5*a2(1)/3;
    Y2(2)=2*a1(2)/3+5*a2(2)/3;

    X3(1)=4*a1(1)/3+a2(1)/3;
    Y3(1)=4*a1(2)/3+a2(2)/3;
    X3(2)=5*a1(1)/3+2*a2(1)/3;
    Y3(2)=5*a1(2)/3+2*a2(2)/3;

    X4(1)=4*a1(1)/3+4*a2(1)/3;
    Y4(1)=4*a1(2)/3+4*a2(2)/3;
    X4(2)=5*a1(1)/3+5*a2(1)/3;
    Y4(2)=5*a1(2)/3+5*a2(2)/3;

    plot(X1,Y1,'k.','Markersize',10);hold on;
    plot(X2,Y2,'k.','Markersize',10);hold on;
    plot(X3,Y3,'k.','Markersize',10);hold on;
    plot(X4,Y4,'k.','Markersize',10);hold on;
    
end