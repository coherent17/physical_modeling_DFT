%cut the surface of the z=10.591667
clc;
clear;
tic;

z=10.591667;

%n=1
% eigenvector1=[0.21459,-0.00001];
% eigenvector2=[-0.15951, 0.27631];
% eigenvector3=[0.30512,-0.52854];
% eigenvector4=[-0.24113,-0.41758];
% eigenvector5=[0.20353,-0.35253];
% eigenvector6=[-0.28545,-0.00000];

%n=2
% eigenvector1=[0.51172, 0.00000];
% eigenvector2=[-0.24973, 0.43260];
% eigenvector3=[0.02362,-0.04092];
% eigenvector4=[-0.04289,-0.07428];
% eigenvector5=[-0.24656, 0.42701];
% eigenvector6=[0.48571, 0.00000];

%n=3
% eigenvector1=[0.37956, 0.21915];
% eigenvector2=[0.19283, 0.33396];
% eigenvector3=[-0.30654, 0.17697];
% eigenvector4=[0.51005, 0.00000];
% eigenvector5=[-0.15095,-0.26147];
% eigenvector6=[-0.37007,-0.21366];

%n=4
eigenvector1=[0.10728, 0.18585];
eigenvector2=[0.31904,-0.00001];
eigenvector3=[0.61031, 0.00000];
eigenvector4=[-0.24109, 0.41757];
eigenvector5=[-0.40708, 0.00000];
eigenvector6=[-0.14274,-0.24720];

%n=5
% eigenvector1=[0.25583, 0.44317];
% eigenvector2=[0.49951, 0.00000];
% eigenvector3=[0.04726,-0.00001];
% eigenvector4=[-0.04289, 0.07427];
% eigenvector5=[0.49309, 0.00000];
% eigenvector6=[0.24286, 0.42063];

%n=6
% eigenvector1=[-0.00001, 0.43830];
% eigenvector2=[0.19279,-0.33395];
% eigenvector3=[-0.30657,-0.17698];
% eigenvector4=[-0.25503,-0.44169];
% eigenvector5=[-0.15097, 0.26148];
% eigenvector6=[-0.00000,-0.42732];

data1=importdata("k=( 0.849, 0.000, 0.000)_1BS.dat");
data2=importdata("k=( 0.849, 0.000, 0.000)_2BS.dat");
data3=importdata("k=( 0.849, 0.000, 0.000)_3BS.dat");
data4=importdata("k=( 0.849, 0.000, 0.000)_4BS.dat");
data5=importdata("k=( 0.849, 0.000, 0.000)_5BS.dat");
data6=importdata("k=( 0.849, 0.000, 0.000)_6BS.dat");

% %--- find the z
% list_z = unique(data1.data(:,3));

toc;

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
                      (eigenvector6(1)+eigenvector6(2)*1i)*(data6.data(i,4)+data6.data(i,5)*1i);
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
axis([-8 8 -2 14]);

figure(2)
scatter(x1(1,1:flag1),y1(1,1:flag1),5,imag(psi1(1,1:flag1)));hold on;
colorbar
title(colorbar,'Im(Ψ)');
plot_carbon_atoms;
title(['(n=4)Im(Ψ) at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -2 14]);

figure(3)
scatter(x1(1,1:flag1),y1(1,1:flag1),5,(real(psi1(1,1:flag1)).^2+imag(psi1(1,1:flag1)).^2));hold on;
colorbar
title(colorbar,'|Ψ|^2');
plot_carbon_atoms;
title(['(n=4)|Ψ|^2 at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -2 14]);

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

    plot(X1(2),Y1(2),'k.','Markersize',10);hold on;
    plot(X2,Y2,'k.','Markersize',10);hold on;
    plot(X3(2),Y3(2),'k.','Markersize',10);hold on;
    plot(X4,Y4,'k.','Markersize',10);hold on;
end