clc;
clear;

%real space
figure (1)
%supercell(2*2)
S_a0=2.4669*2;
S_a1=S_a0*[0.5,sqrt(3)/2,0];
S_a2=S_a0*[-0.5,sqrt(3)/2,0];
S_a3=S_a0*[0,0,20.5];
S_X1(1)=S_a1(1);
S_Y1(1)=S_a1(2);
S_X1(2)=S_a2(1);
S_Y1(2)=S_a2(2);
DrawLine(0,0,S_X1(1),S_Y1(1),'g',1);
DrawLine(0,0,S_X1(2),S_Y1(2),'g',1);
text(10*S_a1(1)/9-0.25,10*S_a1(2)/9,'S a_1');
text(10*S_a2(1)/9,10*S_a2(2)/9,'S a_2');
%normal cell
a0=2.4669;
a1=a0*[0.5,sqrt(3)/2,0];
a2=a0*[-0.5,sqrt(3)/2,0];
a3=a0*[0,0,20.5];
X1(1)=a1(1);
Y1(1)=a1(2);
X1(2)=a2(1);
Y1(2)=a2(2);
DrawLine(0,0,X1(1),Y1(1),'r',5);
DrawLine(0,0,X1(2),Y1(2),'r',5);
text(10*a1(1)/9-0.25,10*a1(2)/9,'a1');
text(10*a2(1)/9,10*a2(2)/9,'a2');

title("Real space lattice constant (2*2 supercell)");
xlabel("x(Å)");
ylabel("y(Å)");
axis([-4 4 0 8]);


%reciprocal space
figure(2)
%normal cell
omega=dot(cross(a1,a2),a3);
b1=2*pi*cross(a2,a3)/omega;
b2=2*pi*cross(a3,a1)/omega;
b3=2*pi*cross(a1,a2)/omega;
X2(1)=b1(1);
Y2(1)=b1(2);
X2(2)=b2(1);
Y2(2)=b2(2);
DrawLine(0,0,X2(1),Y2(1),'r',1);
DrawLine(0,0,X2(2),Y2(2),'r',1);
text(10*b1(1)/9-0.25,10*b1(2)/9,'b1');
text(10*b2(1)/9,10*b2(2)/9,'b2');
%supercell(2*2)
S_omega=dot(cross(S_a1,S_a2),S_a3);
S_b1=2*pi*cross(S_a2,S_a3)/S_omega;
S_b2=2*pi*cross(S_a3,S_a1)/S_omega;
S_b3=2*pi*cross(S_a1,S_a2)/S_omega;
S_X2(1)=S_b1(1);
S_Y2(1)=S_b1(2);
S_X2(2)=S_b2(1);
S_Y2(2)=S_b2(2);
DrawLine(0,0,S_X2(1),S_Y2(1),'g',5);
DrawLine(0,0,S_X2(2),S_Y2(2),'g',5);
text(5*S_b1(1)/4-0.45,5*S_b1(2)/4+0.25,'S b_1');
text(5*S_b2(1)/4,5*S_b2(2)/4+0.25,'S b_2');

title("Reciprocal space lattice constant (2*2 supercell)");
xlabel("x(Å^-^1)");
ylabel("y(Å^-^1)");
axis([-4 4 0 8])

function DrawLine(x1,y1,x2,y2,color,w)
    plot(x1,y1);hold on;
    plot(x2,y2);hold on;
    line([x1,x2],[y1,y2],'Color',color,'LineWidth',w);hold on;
end