%cut the surface of the z=10.591667
clc;
clear;
tic;

z=10.591667;

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
% eigenvector1=[-0.20576, 0.35515];
% eigenvector2=[0.42646, 0.00060];
% eigenvector3=[0.18567, 0.32160];
% eigenvector4=[0.01484,-0.02573];
% eigenvector5=[0.44208, 0.00000];
% eigenvector6=[0.20210,-0.35101];
% eigenvector7=[-0.38909, 0.00050];

%n=5
eigenvector1=[0.09336, 0.16280];
eigenvector2=[0.38897, 0.00028];
eigenvector3=[-0.22104, 0.38284];
eigenvector4=[-0.24828,-0.43010];
eigenvector5=[-0.18567,-0.32160];
eigenvector6=[0.11392, 0.19625];
eigenvector7=[0.42665,-0.00027];

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
                      (eigenvector6(1)+eigenvector6(2)*1i)*(data6.data(i,4)+data6.data(i,5)*1i)+ ...
                      (eigenvector7(1)+eigenvector7(2)*1i)*(data7.data(i,4)+data7.data(i,5)*1i);
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
title(['(n=5)Re(Ψ) at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -2 14])

figure(2)
scatter(x1(1,1:flag1),y1(1,1:flag1),5,imag(psi1(1,1:flag1)));hold on;
colorbar
title(colorbar,'Im(Ψ)');
plot_carbon_atoms;
title(['(n=5)Im(Ψ) at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -2 14])


figure(3)
scatter(x1(1,1:flag1),y1(1,1:flag1),5,(real(psi1(1,1:flag1)).^2+imag(psi1(1,1:flag1)).^2));hold on;
colorbar
title(colorbar,'|Ψ|^2');
plot_carbon_atoms;
title(['(n=5)|Ψ|^2 at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -2 14])

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
    plot(X3,Y3,'k.','Markersize',10);hold on;
    plot(X4,Y4,'k.','Markersize',10);hold on;
    
end