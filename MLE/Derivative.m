pkg load optim
clear
%derivative1
f=@(x) 0.3*x.^3-0.5*x.^2-3*x+2
df = @(x) 0.9*x.^2-x-3
%Step
xmin=0
xmax=5
x=xmin:(xmax-xmin)/20:xmax;
ylim([-4,12])

%a
x0=4.5
grid
subplot(111)
L=2
xlabel('variable x')
ylabel('function f(x)')
xlim([xmin,xmax])
hold
plot(x,f(x),'b','linewidth',2)
plot(x0,f(x0),'r','markersize',15)
plot([x0-L,x0+L],[f(x0)-df(x0)*L,f(x0)+df(x0)*L],'r','linewidth',1.5)
plot(x0-[f(x0)/df(x0)],0,'g','markersize',15)
hold off

%b
x0=1.
grid
subplot(111)
L=1 
xlabel('variable x')
ylabel('function f(x)')
xlim([xmin,xmax])
hold
plot(x,f(x),'b','linewidth',2)
plot(x0,f(x0),'r','markersize',15)
plot([x0-L,x0+L],[f(x0)-df(x0)*L,f(x0)+df(x0)*L],'r','linewidth',1.5)
plot(x0-[f(x0)/df(x0)],0,'g','markersize',15)
hold off


clear
%derivative2
g=9.81
L=5
H=3
beta=-20

%solve a & b
%10a +b = tan(beta) , 30/25 - 50b/25 + b = tan(beta)
%25a + 5b = 3 , a= (3-5b)/25
b= ((tand(beta)-2*L*(H/25)))/-1
a= (H-L*b)/L^2

%solve alpha and v0
alpha = atand(b)
v0 = sqrt(-g/(2*a*(cosd(alpha))^2))

x = 0:L/100:L;
y = @(x) a*x.^2 + b*x
dy = @(x) 2*a*x+b
grid



plot(x, y(x),'b','linewidth',1.5)
hold
plot([0,L,L],[0,0,y(L)],'k','linewidth',1)
plot([L,L+1],[y(L),y(L+1)],'r','linewidth',1.5)
hold off





clear
%derivative3


%Matrix values
x1 = 1
y1 = 2
k1 = -0.25
m1 = 1
x2 = 5
y2 = 4
k2 = 0.5
m2 = -2

%A leftside matrix
A = [x1^5,x1^4,x1^3,x1^2,x1,1
      5*x1^4,4*x1^3,3*x1^2,2*x1,1,0
      20*x1^3,12*x1^2,6*x1,2,0,0
      x2^5,x2^4,x2^3,x2^2,x2,1
      5*x2^4,4*x2^3,3*x2^2,2*x2,1,0
      20*x2^3,12*x2^2,6*x2,2,0,0]
      
%B right handside
B =[
    y1
    k1
    m1
    y2
    k2
    m2]
      
%solution i.e the coefficients a,b,c,d,e,f
sol=A\B
%%
a=sol(1)
b=sol(2)
c=sol(3)
d=sol(4)
e=sol(5)
f=sol(6)

% y and x
x=1:5/100:5
y = @(x) a*x.^5 + b*x.^4 +c*x.^3 + d*x.^2 + e*x + f
dy = @(x) 5*a*x.^4 +4*b*x.^3 +3*c*x.^2 +2*d*x + e
d2y = @(x) 20*a*x.^3 + 12*b*x.^2 + 6*c*x + 2*d
k = @(x) d2y(x)./(sqrt(1+dy(x).^2).^3)


figure(1)
grid
plot(x,dy(x),'b','linewidth',1.5)
hold on
plot(x2,dy(x2),'g','markersize',15)
plot(x1,dy(x1),'r','markersize',15)
hold off
title('derivative')


figure(2)
grid
plot(x,d2y(x),'b','linewidth',1.5)
hold on
plot(x2,d2y(x2),'g','markersize',15)
plot(x1,d2y(x1),'r','markersize',15)
hold off
title('second derivative')

figure(3)
grid
plot(x,k(x),'b','linewidth',1.5)
hold on
plot(x2,k(x2),'g','markersize',15)
plot(x1,k(x1),'r','markersize',15)
hold off
title('kaarevuus')

figure(4)
%k = @(x) d2y(x)./((sqrt(1+dy(x).^2)).^3)

%radius of curvator circle
R = @(x) abs(1./k(x))
R1 = R(x1)
R2 = R(x2)
%polar angle of the tangent i.e vector [1,f'(x0)]
th= @(x)atan2d(dy(x),1) 
th1 = th(x1)
th2 = th(x2)
%center of the circle of curvature
k(x1)
if k(x1)>=0
    xc = @(x) x +R(x).*cosd(th(x)+90)
    yc = @(x) y(x) +R(x).*sind(th(x)+90)
else
    xc = @(x) x +R(x).*cosd(th(x)-90)
    yc = @(x) y(x) +R(x).*sind(th(x)-90)
end

%points on the circle of curvature
t=0:360; 
kx= @(x) xc(x)+R(x)*cosd(t);
ky= @(x) yc(x)+R(x)*sind(t);

plot(x,y(x),'b','linewidth',1.5)
grid
hold on
plot(kx(x1),ky(x1),'r','linewidth',1.5)
plot(xc(x1),yc(x1),'g','linewidth',1.5)
L=0.5 %horizontal width =2L
plot([x1,x1+L],[y(x1),y(x1)+dy(x1)*L],'r','linewidth',1.5)
plot([x1,xc(x1)],[y(x1),yc(x1)],'c','linewidth',1.5)
plot(x1,y(x1),'b.','markersize',20)
plot(xc(x1),yc(x1),'g.','markersize',20)
hold off

k(x2)
if k(x2)>=0
    xc = @(x) x +R(x).*cosd(th(x)+90)
    yc = @(x) y(x) +R(x).*sind(th(x)+90)
else
    xc = @(x) x +R(x).*cosd(th(x)-90)
    yc = @(x) y(x) +R(x).*sind(th(x)-90)
end
%points on the circle of curvature
t=0:360; 
kx= @(x) xc(x)+R(x)*cosd(t);
ky= @(x) yc(x)+R(x)*sind(t);

hold on
plot(kx(x2),ky(x2),'g','linewidth',1.5)
plot(xc(x2),yc(x2),'g','linewidth',1.5)
L=0.5 %horizontal width =2L
plot([x2,x2+L],[y(x2),y(x2)+dy(x2)*L],'g','linewidth',1.5)
plot([x2,xc(x2)],[y(x2),yc(x2)],'c','linewidth',1.5)
plot(x2,y(x2),'b.','markersize',20)
plot(xc(x2),yc(x2),'g.','markersize',20)
hold off

clear
%derivative4
A=1
a=1.5
gamma = 1.57
B=1.5
b=1
T=12.57
t=0:T/100:T
t0 = 4.5

%normal function
x= @(t) A*sin(a.*t + gamma)
y = @(t) B*sin(b.*t)

%derivative function
dx = @(t) A*a*cos(a.*t + gamma)
dy = @(t) B*b*cos(b.*t)

%second derivative
d2x = @(t) -A*a*a*sin(a.*t + gamma)
d2y = @(t) -B*b*b*sin(b.*t)

%Curvature function
k = @(t) (dx(t).*d2y(t)-d2x(t).*dy(t))./(sqrt((dx(t)).^2+(dy(t)).^2)).^3

%Radius
R = @(t) 1./abs(k(t))

%theta
th = @(t) atan2d(dy(t),dx(t))

%tagent vector T 
Tl = @(t) sqrt((dy(t)).^2+(dx(t)).^2)

figure(1)
grid
plot(t/ pi,k(t),'b','linewidth',2)

figure(2)
grid
plot(t,Tl(t),'b','linewidth',2)

figure(3)
k(t0)
R(t0)
if k(t0)>=0
    x0 =  x(t0) +R(t0).*cosd(th(t0)+90)
    y0 =  y(t0) +R(t0).*sind(th(t0)+90)
else
    x0 =  x(t0) +R(t0).*cosd(th(t0)-90)
    y0 =  y(t0) +R(t0).*sind(th(t0)-90)
end
%points on the circle of curvature
d=0:360; 
kx0=  x0+R(t0)*cosd(d);
ky0=  y0+R(t0)*sind(d);

grid
plot(x(t),y(t),'b','linewidth',2)
hold on
plot(kx0,ky0,'g','linewidth',1.5)
plot(x0,y0,'g','linewidth',1.5)
L=0.5 %horizontal width =2L
plot([x(t0),x(t0+L)],[y(t0),y(t0)+dy(t0)],'r','linewidth',1.5)
plot([x(t0),x0],[y(t0),y0],'c','linewidth',1.5)
plot(x(t0),y(t0),'b.','markersize',20)
plot(x0,y0,'g.','markersize',20)
hold off





