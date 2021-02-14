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

%n=1 part
flag1=0;
for i=1:+1:length(data1.data(:,1))
    if data1.data(i,1)==0
        flag1=flag1+1;
        psi1(1,flag1)=(eigenvector1(1)+eigenvector1(2)*1i)*(data1.data(i,4)+data1.data(i,5)*1i)+(eigenvector2(1)+eigenvector2(2)*1i)*(data2.data(i,4)+data2.data(i,5)*1i);
        y1(1,flag1)=data1.data(i,2);
        z1(1,flag1)=data1.data(i,3);
    else 
        flag1=flag1+0;
    end
end

figure(1)
scatter(y1(1,1:flag1),z1(1,1:flag1),5,real(psi1(1,1:flag1)));hold on;
colorbar
title(colorbar,'Re(Ψ)');
title('(n=1)Re(Ψ) at K point on the surface x=0(Å)');
xlabel('y(Å)');
ylabel('z(Å)');


figure(2)
scatter(y1(1,1:flag1),z1(1,1:flag1),5,imag(psi1(1,1:flag1)));hold on;
colorbar
title(colorbar,'Im(Ψ)');
title('(n=1)Im(Ψ) at K point on the surface x=0(Å)');
xlabel('y(Å)');
ylabel('z(Å)');

figure(3)
scatter(y1(1,1:flag1),z1(1,1:flag1),5,sqrt(real(psi1(1,1:flag1)).^2+imag(psi1(1,1:flag1)).^2).^2);hold on;
colorbar
title(colorbar,'|Ψ|^2');
title('(n=1)|Ψ|^2 at K point on the surface x=0(Å)');
xlabel('y(Å)');
ylabel('z(Å)');

%n=2 part
flag2=0;
for i=1:+1:length(data3.data(:,1))
    if data3.data(i,1)==0
        flag2=flag2+1;
        psi2(1,flag2)=(eigenvector3(1)+eigenvector3(2)*1i)*(data3.data(i,4)+data3.data(i,5)*1i)+(eigenvector4(1)+eigenvector4(2)*1i)*(data4.data(i,4)+data4.data(i,5)*1i);
        y2(1,flag2)=data3.data(i,2);
        z2(1,flag2)=data3.data(i,3);
    else 
        flag2=flag2+0;
    end
end

figure(4)
scatter(y2(1,1:flag2),z2(1,1:flag2),5,real(psi2(1,1:flag2)));hold on;
colorbar
title(colorbar,'Re(Ψ)');
title('(n=2)Re(Ψ) at K point on the surface x=0(Å)');
xlabel('y(Å)');
ylabel('z(Å)');


figure(5)
scatter(y2(1,1:flag2),z2(1,1:flag2),5,imag(psi2(1,1:flag2)));hold on;
colorbar
title(colorbar,'Im(Ψ)');
title('(n=2)Im(Ψ) at K point on the surface x=0(Å)');
xlabel('y(Å)');
ylabel('z(Å)');


figure(6)
scatter(y2(1,1:flag2),z2(1,1:flag2),5,sqrt(real(psi2(1,1:flag2)).^2+imag(psi2(1,1:flag2)).^2).^2);hold on;
colorbar
title(colorbar,'|Ψ|^2');
title('(n=2)|Ψ|^2 at K point on the surface x=0(Å)');
xlabel('y(Å)');
ylabel('z(Å)');


