%cut the surface of the z=___
clc;
clear;

z=11.104167;
eigenvector1=[0.72292,0];
eigenvector2=[0.68996,0.03650];
eigenvector3=[-0.68996,0.03650];
eigenvector4=[0.72292,0];

data1=importdata("k=( 1.698, 0.000, 0.000)_1BS_n1.dat");
data2=importdata("k=( 1.698, 0.000, 0.000)_2BS_n1.dat");
data3=importdata("k=( 1.698, 0.000, 0.000)_1BS_n2.dat");
data4=importdata("k=( 1.698, 0.000, 0.000)_2BS_n2.dat");

%n=1 part
flag1=0;
for i=1:+1:length(data1.data(:,1))
    if data1.data(i,3)==z
        flag1=flag1+1;
        psi1(1,flag1)=(eigenvector1(1)+eigenvector1(2)*1i)*(data1.data(i,4)+data1.data(i,5)*1i)+(eigenvector2(1)+eigenvector2(2)*1i)*(data2.data(i,4)+data2.data(i,5)*1i);
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
title(['(n=1)Re(Ψ) at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -8 8]);

figure(2)
scatter(x1(1,1:flag1),y1(1,1:flag1),5,imag(psi1(1,1:flag1)));hold on;
colorbar
title(colorbar,'Im(Ψ)');
plot_carbon_atoms;
title(['(n=1)Im(Ψ) at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -8 8]);

figure(3)
scatter(x1(1,1:flag1),y1(1,1:flag1),5,sqrt(real(psi1(1,1:flag1)).^2+imag(psi1(1,1:flag1)).^2).^2);hold on;
% scatter(x1(1,1:flag1),y1(1,1:flag1),5,(real(psi1(1,1:flag1)).^2+imag(psi1(1,1:flag1)).^2));hold on;
colorbar
title(colorbar,'|Ψ|^2');
plot_carbon_atoms;
title(['(n=1)|Ψ|^2 at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -8 8]);

%n=2 part
flag2=0;
for i=1:+1:length(data3.data(:,1))
    if data3.data(i,3)==z
        flag2=flag2+1;
        psi2(1,flag2)=(eigenvector3(1)+eigenvector3(2)*1i)*(data3.data(i,4)+data3.data(i,5)*1i)+(eigenvector4(1)+eigenvector4(2)*1i)*(data4.data(i,4)+data4.data(i,5)*1i);
        x2(1,flag2)=data1.data(i,1);
        y2(1,flag2)=data1.data(i,2);
    else 
        flag2=flag2+0;
    end
end

figure(4)
scatter(x2(1,1:flag2),y2(1,1:flag2),5,real(psi2(1,1:flag2)));hold on;
colorbar
title(colorbar,'Re(Ψ)');
plot_carbon_atoms;
title(['(n=2)Re(Ψ) at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -8 8]);

figure(5)
scatter(x2(1,1:flag2),y2(1,1:flag2),5,imag(psi2(1,1:flag2)));hold on;
colorbar
title(colorbar,'Im(Ψ)');
plot_carbon_atoms;
title(['(n=2)Im(Ψ) at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -8 8]);

figure(6)
scatter(x2(1,1:flag2),y2(1,1:flag2),5,sqrt(real(psi2(1,1:flag2)).^2+imag(psi2(1,1:flag2)).^2).^2);hold on;
colorbar
title(colorbar,'|Ψ|^2');
plot_carbon_atoms;
title(['(n=2)|Ψ|^2 at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -8 8]);

function plot_carbon_atoms
    a0=2.4669;
    a1=a0*[0.5,sqrt(3)/2,0];
    a2=a0*[-0.5,sqrt(3)/2,0];
    X1(1)=(1/3)*(a1(1)+a2(1));
    Y1(1)=(1/3)*(a1(2)+a2(2));
    X2(1)=(2/3)*(a1(1)+a2(1));
    Y2(1)=(2/3)*(a1(2)+a2(2));
    X3(1)=a1(1);
    Y3(1)=a1(2);
    X4(1)=X3(1)+a2(1);
    Y4(1)=Y3(1)+a1(2);
    for i=1:+1:5
        theta=(60*i);
        X1(i+1)=X1(1)*cosd(theta)-Y1(1)*sind(theta);
        Y1(i+1)=Y1(1)*cosd(theta)+X1(1)*sind(theta);
        X2(i+1)=X2(1)*cosd(theta)-Y2(1)*sind(theta);
        Y2(i+1)=Y2(1)*cosd(theta)+X2(1)*sind(theta);
        X3(i+1)=X3(1)*cosd(theta)-Y3(1)*sind(theta);
        Y3(i+1)=Y3(1)*cosd(theta)+X3(1)*sind(theta);
        X4(i+1)=X4(1)*cosd(theta)-Y4(1)*sind(theta);
        Y4(i+1)=Y4(1)*cosd(theta)+X4(1)*sind(theta);
    end
    plot(X1,Y1,'k.','MarkerSize',20);hold on;
    plot(X2,Y2,'k.','MarkerSize',20);hold on;
    plot(0,0,'r.','MarkerSize',10);hold on;
    plot(X3,Y3,'r.','MarkerSize',10);hold on;
    plot(X4,Y4,'r.','MarkerSize',10);hold on;
end