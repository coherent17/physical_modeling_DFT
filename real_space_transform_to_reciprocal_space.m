clear;
clc;
%---Real space
a1=[0.5,sqrt(3)/2,0];
a2=[-0.5,sqrt(3)/2,0];
a3=[0,0,20];
X0=0;
Y0=0;
X1=[0,0,0,0,0,0];
Y1=[0,0,0,0,0,0];
X2=[0,0,0,0,0,0];
Y2=[0,0,0,0,0,0];
X3=[0,0,0,0,0,0];
Y3=[0,0,0,0,0,0];
%X4=[0,0,0,0,0,0];
%Y4=[0,0,0,0,0,0];
X5=[0,0,0,0,0,0];
Y5=[0,0,0,0,0,0];
s=[0,0,0,0,0,0];
%--(X1,Y1)is the first circle of the atom C
%--(X2,Y2)is the second circle of the atom C
X1(1)=(1/3)*(a1(1)+a2(1));
Y1(1)=(1/3)*(a1(2)+a2(2));
X2(1)=(2/3)*(a1(1)+a2(1));
Y2(1)=(2/3)*(a1(2)+a2(2));
for i=1:+1:5
    theta=(60*i);
    X1(i+1)=X1(1)*cosd(theta)-Y1(1)*sind(theta);
    Y1(i+1)=Y1(1)*cosd(theta)+X1(1)*sind(theta);
    X2(i+1)=X2(1)*cosd(theta)-Y2(1)*sind(theta);
    Y2(i+1)=Y2(1)*cosd(theta)+X2(1)*sind(theta);
end
%--(X3,Y3)is the lattice point 
X3(1)=a1(1);
Y3(1)=a1(2);
theta=0;
for k=1:+1:5
    theta=(60*k);
    X3(k+1)=X3(1)*cosd(theta)-Y3(1)*sind(theta);
    Y3(k+1)=Y3(1)*cosd(theta)+X3(1)*sind(theta);
end
%--(X4,Y4)is the point on the mid of the vertices of the WS cell
[X4,Y4]=GetMidPoint(X3,Y3,X0,Y0);
for i=1:+1:6
   s(i)=GetSlope(X4(i),Y4(i),0,0); 
end
%--(X5,Y5)use the slope and the point on the curve to solve the intersection 
[X5(1),Y5(1)]=GetPoint(X4(1),Y4(1),X4(2),Y4(2),s(1),s(2));
[X5(2),Y5(2)]=GetPoint(X4(2),Y4(2),X4(3),Y4(3),s(2),s(3));
[X5(3),Y5(3)]=GetPoint(X4(3),Y4(3),X4(4),Y4(4),s(3),s(4));
[X5(4),Y5(4)]=GetPoint(X4(4),Y4(4),X4(5),Y4(5),s(4),s(5));
[X5(5),Y5(5)]=GetPoint(X4(5),Y4(5),X4(6),Y4(6),s(5),s(6));
[X5(6),Y5(6)]=GetPoint(X4(6),Y4(6),X4(1),Y4(1),s(6),s(1));
%---Reciprocal space
omega=dot(cross(a1,a2),a3);
b1=2*pi*cross(a2,a3)/omega;
b2=2*pi*cross(a3,a1)/omega;
b3=2*pi*cross(a1,a2)/omega;
R_X0=0;
R_Y0=0;
R_X3=[0,0,0,0,0,0];
R_Y3=[0,0,0,0,0,0];
%R_X4=[0,0,0,0,0,0];
%R_Y4=[0,0,0,0,0,0];
R_X5=[0,0,0,0,0,0];
R_Y5=[0,0,0,0,0,0];
R_s=[0,0,0,0,0,0];
%---(R_X3,R_Y3)is the lattice point in the reciprocal space
R_X3(1)=b1(1);
R_Y3(1)=b1(2);
for k=1:+1:5
    R_theta=(60*k);
    R_X3(k+1)=R_X3(1)*cosd(R_theta)-R_Y3(1)*sind(R_theta);
    R_Y3(k+1)=R_Y3(1)*cosd(R_theta)+R_X3(1)*sind(R_theta);
end
[R_X4,R_Y4]=GetMidPoint(R_X3,R_Y3,R_X0,R_Y0);
for i=1:+1:6
   R_s(i)=GetSlope(R_X4(i),R_Y4(i),0,0); 
end
[R_X5(1),R_Y5(1)]=GetPoint(R_X4(1),R_Y4(1),R_X4(2),R_Y4(2),R_s(1),R_s(2));
[R_X5(2),R_Y5(2)]=GetPoint(R_X4(2),R_Y4(2),R_X4(3),R_Y4(3),R_s(2),R_s(3));
[R_X5(3),R_Y5(3)]=GetPoint(R_X4(3),R_Y4(3),R_X4(4),R_Y4(4),R_s(3),R_s(4));
[R_X5(4),R_Y5(4)]=GetPoint(R_X4(4),R_Y4(4),R_X4(5),R_Y4(5),R_s(4),R_s(5));
[R_X5(5),R_Y5(5)]=GetPoint(R_X4(5),R_Y4(5),R_X4(6),R_Y4(6),R_s(5),R_s(6));
[R_X5(6),R_Y5(6)]=GetPoint(R_X4(6),R_Y4(6),R_X4(1),R_Y4(1),R_s(6),R_s(1));
%---Real space
figure (1)
axis([-1.5 1.5 -1.5 1.5]);
line([a1(1),X0],[a1(2),Y0],'Color','b','LineWidth',5);hold on;
line([a2(1),X0],[a2(2),Y0],'Color','b','LineWidth',5);hold on;
text(5*a1(1)/4,5*a1(2)/4,'a1');
text(5*a2(1)/4,5*a2(2)/4,'a2');
plot(X0,Y0,'g.','MarkerSize',50);hold on;
plot(X1,Y1,'k.','MarkerSize',80);hold on;
plot(X2,Y2,'k.','MarkerSize',80);hold on;
DrawLine(X1,Y1,'k',1);
for j=1:+1:6
    line([X1(j),X2(j)],[Y1(j),Y2(j)],'Color','k','LineWidth',4);hold on;
end
plot(X3,Y3,'g.','MarkerSize',50);hold on;
plot(X4,Y4,'bx','MarkerSize',10);hold on;
plot(X5,Y5,'bo','MarkerSize',20);hold on;
DrawLine(X5,Y5,'r',10);

%---Reciprocal space
figure(2)
axis([-8 8 -8 8])
plot(R_X0,R_Y0,'go','MarkerSize',8);hold on;
line([b1(1),R_X0],[b1(2),R_Y0],'Color','b','LineWidth',5);hold on;
line([b2(1),R_X0],[b2(2),R_Y0],'Color','b','LineWidth',5);hold on;
text(7*b1(1)/6,7*b1(2)/6,'b1');
text(7*b2(1)/6,7*b2(2)/6,'b2');
plot(R_X3,R_Y3,'g.','MarkerSize',50);hold on;
plot(R_X4,R_Y4,'bx','MarkerSize',10);hold on;
plot(R_X5,R_Y5,'k.','MarkerSize',80);hold on;
DrawLine(R_X5,R_Y5,'r',10);

%-----function part
function DrawLine(x,y,color,w)
    line(x,y,'Color',color,'LineWidth',w);hold on;
    line([x(1),x(6)],[y(1),y(6)],'Color',color,'LineWidth',w);hold on;
end
function [x,y]=GetMidPoint(a,b,c,d)
    x=(a+c)/2;
    y=(b+d)/2;
end
function s=GetSlope(a,b,c,d)
    s=-(a-c)/(b-d);
end
function [x,y]=GetPoint(x1,y1,x2,y2,s1,s2)
    if s1==Inf
        x=x1;
        y=y2-s2*(x2-x1);
    elseif s2==Inf
        x=x2;
        y=y1-s1*(x1-x2);
    elseif s2==-Inf
        x=x2;
        y=y1+s1*(x2-x1);
    elseif s1==-Inf
        x=x1;
        y=y2+s2*(x1-x2);
    else
        delta=(s2-s1)/(s1*s2);
        delta_x=(x2-y2/s2)/s1-(x1-y1/s1)/s2;
        delta_y=(x2-y2/s2)-(x1-y1/s1);
        x=delta_x/delta;
        y=delta_y/delta;
    end
end