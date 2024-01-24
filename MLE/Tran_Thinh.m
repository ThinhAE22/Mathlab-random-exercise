%List of excercise
%Acceleration velocity position: 1,2
%Derivatative: 1,2,3,4
%Extrume: 1,2,3,4,5
%Integral: 1,2,3,4,5,6
%position velocity acceleration: 1,2,3,5,6,4 
%(I choose the easy ex to do first)
%Taylor Polynomial: 1

%-----------------------------------------------------%
%AVP ex
%ex 1
%initial
clear
close all
T1 = 1
T2 = 2
T3 = 3
amax = 2
T = 4*T1 + 2*T2 + T3

%basic calculation

vmax = (T1 + T2)*amax
smax = (2*T1 +T2 + T3)*vmax

%Formular

%accel formular
%T11
t11 = 0:T1/100:T1;
aT11 = @(t11) 2*t11

%T21
t21 = T1:T2/100:T1+T2;
aT21 = @(t21) amax

%T12
t12 = T1+T2:T1/100:T2+2*T1;
aT12 = @(t12) -2*(t12-4)

%T3
t3 = T2+2*T1:T3/100:T2+2*T1+T3;
aT3 = @(t3) 0


%aT13
t13 = T2+2*T1+T3:T1/100:T2+3*T1+T3;
aT13 = @(t13) -2*(t13-7)

%T22
t22 = T2+3*T1+T3:T2/100:2*T2+3*T1+T3;
aT22 = @(t22) -amax

%aT14
t14 = 2*T2+3*T1+T3:T1/100:2*T2+4*T1+T3;
aT14 = @(t14) 2*(t14-11)





%Speed equation
%T11
t11 = 0:T1/100:T1;
vT11 = @(t11) t11.^2

%T21
t21 = T1:T2/100:T1+T2;
vT21 = @(t21) amax*t21 
vT21f = @(t21) amax*t21 + vT11(T1) -vT21(T1)
%T12
t12 = T1+T2:T1/100:T2+2*T1;
vT12 = @(t12) -t12.*(t12-8) 
vT12f = @(t12) -t12.*(t12-8) + vT21f(T1+T2) -vT12(T1+T2)
%T3
t3 = T2+2*T1:T3/100:T2+2*T1+T3;
vT3 = @(t3) vmax


%T13
t13 = T2+2*T1+T3:T1/100:T2+3*T1+T3;
vT13 = @(t13) -t13.*(t13-14) 
vT13f = @(t13) -t13.*(t13-14) + vmax - vT13(T2+2*T1+T3)

%T22
t22 = T2+3*T1+T3:T2/100:2*T2+3*T1+T3;
vT22 = @(t22) -amax*t22 
vT22f = @(t22) -amax*t22 + vT13f(T2+3*T1+T3)-vT22(T2+3*T1+T3)
%T14
t14 = 2*T2+3*T1+T3:T1/100:2*T2+4*T1+T3;
vT14 = @(t14) t14.*(t14-22) 
vT14f = @(t14) t14.*(t14-22) + vT22f(2*T2+3*T1+T3)-vT14(2*T2+3*T1+T3)


%Location
%T11
t11 = 0:T1/100:T1;
sT11 = @(t11) (t11.^3)./3
%T21
t21 = T1:T2/100:T1+T2;
sT21 = @(t21) (t21.*(amax*t21-2))./2 
sT21f = @(t21) (t21.*(amax*t21-2))./2  + sT11(T1)-sT21(T1)
%T12
t12 = T1+T2:T1/100:T2+2*T1;
sT12 = @(t12) -(t12.*(t12.^2-12*t12+30))./3  
sT12f = @(t12) -(t12.*(t12.^2-12*t12+30))./3 + sT21f(T1+T2) - sT12(T1+T2)

%T3
t3 = T2+2*T1:T3/100:T2+2*T1+T3;
sT3 = @(t3) vmax*t3
sT3f = @(t3) vmax*t3 + sT12f(T2+2*T1) - sT3(T2+2*T1)
%T13
t13 = T2+2*T1+T3:T1/100:T2+3*T1+T3;
sT13 = @(t13) -(t13.*(t13.^2-21*t13+129))./3
sT13f = @(t13) -(t13.*(t13.^2-21*t13+129))./3 + sT3f(T2+2*T1+T3)- sT13(T2+2*T1+T3)
%T22
t22 = T2+3*T1+T3:T2/100:2*T2+3*T1+T3;
sT22 = @(t22) -(t22.*(amax*t22-42))./2
sT22f  = @(t22) -(t22.*(amax*t22-42))./2 + sT13f(T2+3*T1+T3)-sT22(T2+3*T1+T3)
%T14
t14 = 2*T2+3*T1+T3:T1/100:2*T2+4*T1+T3;
sT14 = @(t14) (t14.*(t14.^2-33*t14+363))./3
sT14f = @(t14) (t14.*(t14.^2-33*t14+363))./3 + sT22f(2*T2+3*T1+T3)-sT14(2*T2+3*T1+T3)






%graph
subplot(3,1,1)
plot(t11,aT11(t11),'b','linewidth',1.5)
hold
plot(t21,aT21(t21),'b','linewidth',1.5)
plot(t12,aT12(t12),'b','linewidth',1.5)
plot(t3,aT3(t3),'b','linewidth',1.5)
plot(t13,aT13(t13),'b','linewidth',1.5)
plot(t22,aT22(t22),'b','linewidth',1.5)
plot(t14,aT14(t14),'b','linewidth',1.5)
title('acceleration')
ylim([-amax-1,amax+1])
xlim([0,2*T2+4*T1+T3])
grid
hold off

subplot(3,1,2)
plot(t11,vT11(t11),'b','linewidth',1.5)
hold
plot(t21,vT21f(t21),'b','linewidth',1.5)
plot(t12,vT12f(t12),'b','linewidth',1.5)
plot(t3,vT3(t3),'b','linewidth',1.5)
plot(t13,vT13f(t13),'b','linewidth',1.5)
plot(t22,vT22f(t22),'b','linewidth',1.5)
plot(t14,vT14f(t14),'b','linewidth',1.5)
ylim([0,vmax+1])
xlim([0,2*T2+4*T1+T3])
title('speed')
grid
hold off

subplot(3,1,3)
plot(t11,sT11(t11),'b','linewidth',1.5)
hold
plot(t21,sT21f(t21),'b','linewidth',1.5)
plot(t12,sT12f(t12),'b','linewidth',1.5)
plot(t3,sT3f(t3),'b','linewidth',1.5)
plot(t13,sT13f(t13),'b','linewidth',1.5)
plot(t22,sT22f(t22),'b','linewidth',1.5)
plot(t14,sT14f(t14),'b','linewidth',1.5)
ylim([0,(smax + 1)])
xlim([0,2*T2+4*T1+T3])
title('location')
grid
hold off


%ex 2
clear
close all
%Initial
m = 1
v0 = 5
alfa = 65
b = 1
g = 9.81

%Intial speed
v0x = v0*cosd(alfa)
v0y = v0*sind(alfa)

%i
%Formalar
dt =  0.01
%distant
s1x(1)=0 
s1y(1)=0
%distant
v1x(1) = v0x 
v1y(1) = v0y
k = 1
while  s1y(k) >= 0
      F1=0; %air resistance
      %accel
      a1x(k)=1/m*F1;
      a1y(k) = 1/m*(F1-m*g);
      %position
      s1x(k+1)=s1x(k)+v1x(k)*dt+1/2*a1x(k)*dt^2;
      s1y(k+1)=s1y(k)+v1y(k)*dt+1/2*a1y(k)*dt^2;
      %velocity
      v1x(k+1)=v1x(k)+a1x(k)*dt;
      v1y(k+1)=v1y(k)+a1y(k)*dt;
      k=k+1;
end 

t1=0:dt:(k-1)*dt;
T1=t1(end) 

%ii
%Formalar
dt =  0.01
%distant
s2x(1)=0 
s2y(1)=0
%distant
v2x(1) = v0x 
v2y(1) = v0y
k = 1
while  s2y(k) >= 0
      F2x = -b*v2x(k);
      F2y = -b*v2y(k);
      %accel
      a2x(k)=1/m*F2x;
      a2y(k) = 1/m*(F2y-m*g);
      %position
      s2x(k+1)=s2x(k)+v2x(k)*dt+1/2*a2x(k)*dt^2;
      s2y(k+1)=s2y(k)+v2y(k)*dt+1/2*a2y(k)*dt^2;
      %velocity
      v2x(k+1)=v2x(k)+a2x(k)*dt;
      v2y(k+1)=v2y(k)+a2y(k)*dt;
      k=k+1;
end 

t2=0:dt:(k-1)*dt;
T2=t2(end) 

%iii
%Formalar
dt =  0.01
%distant
s3x(1)=0 
s3y(1)=0
%distant
v3x(1) = v0x 
v3y(1) = v0y
v3(1) = sqrt((v0x)^2+(v0y)^2)
k = 1
while  s3y(k) >= 0
      F3x = -b*v3x(k)*v3(k);
      F3y = -b*v3y(k)*v3(k);
      %accel
      a3x(k)=1/m*F3x;
      a3y(k) = 1/m*(F3y-m*g);
      %position
      s3x(k+1)=s3x(k)+v3x(k)*dt+1/2*a3x(k)*dt^2;
      s3y(k+1)=s3y(k)+v3y(k)*dt+1/2*a3y(k)*dt^2;
      %velocity
      v3x(k+1)=v3x(k)+a3x(k)*dt;
      v3y(k+1)=v3y(k)+a3y(k)*dt;
      v3(k+1)= sqrt((v3x(k+1))^2+(v3y(k+1))^2);
      k=k+1;
end 

t3=0:dt:(k-1)*dt;
T3=t3(end) 


plot(s1x,s1y,'r','linewidth',1.5)
hold
grid
plot(s2x,s2y,'g','linewidth',1.5)
plot(s3x,s3y,'b','linewidth',1.5)
hold off


%-----------------------------------------------------
%Derivatative ex
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

%---------------------------------------------------
%Extrume ex
%Extreme num
%ex 1
clear
close all

%Initial

E = 1
Rt = 15
Rk = 3

R = 0:(Rt*4)/1000:Rt*4;

E0 = 10
Rt0 = 3.6
Rk0 = 0.7

R0 = 0:(Rt0*4)/1000:Rt0*4;


% Formular

%1

u = ((E*R)./(R + Rk)) - ((E*R)./(R + Rt))

funcform = @(R) ((E*R)./(R + Rk)) - ((E*R)./(R + Rt))

%Since the R > 0 
%du = 0 when:
rExtreme = sqrt(Rt*Rk)

%Consider the value of d2u at rExtreme is greater than 0
% rExtreme and u(rExtreme) is the peak of the graph


%2

u0 = ((E0*R0)./(R0 + Rk0)) - ((E0*R0)./(R0 + Rt0))

funcform0 = @(R0) ((E0*R0)./(R0 + Rk0)) - ((E0*R0)./(R0 + Rt0))

%Since the R0 > 0 
%du0 = 0 when:
r0Extreme = sqrt(Rt0*Rk0)

%Consider the value of d2u at rExtreme is less than 0
% r0Extreme and u0(rExtreme) is the peak of the graph


%graph
subplot(2,1,1)
plot(R,u,'linewidth',1.5)
hold
plot(rExtreme,funcform(rExtreme),'r','markersize',20)
hold off

subplot(2,1,2)
plot(R0,u0,'linewidth',1.5)
hold
plot(r0Extreme,funcform0(r0Extreme),'r','markersize',20)
hold off

%ex 2
clear
close all
%Initial
r = 10000
kv = 50
ke = 100
q = kv:1:800;

%Formular

c = @(q) kv*(q./2) + ke*(r./q)

%dc = 0 when
q0 = sqrt(2*ke*(r/kv))
%consider d2c is greater than 0 we have a min valley

%Graph
plot(q,c(q)./10^4,'linewidth',1.5)
hold
plot(q0,c(q0)./10^4,'r','markersize',20)
ylim([0.5,2])
hold off

%ex 3
clear
close all

%Initial
a = 5
b = 9
H = sqrt(a^2 - (b/2)^2)
h = 0:H/100:H
Ax = -b/2
Ay = 0
Bx = b/2
By = 0
Cx = 0
Cy = H

%Formular
AM = @(h) sqrt(h.^2 + (b/2)^2) 
BM = @(h) sqrt(h.^2 + (b/2)^2)
CM = @(h) H - h

s = @(h) AM(h) + BM(h) + CM(h)
%ds = 0 when:
h0 = sqrt(b^2/(4*3))

%Graph
if h0 < H
subplot(2,1,1)
plot(h,s(h),'linewidth',1.5)
hold
plot(h0,s(h0),'r','markersize',20)
hold off
subplot(2,1,2)
plot([Ax,Bx,Cx,Ax],[Ay,By,Cy,Ay],'linewidth',1.5)
hold
plot([Ax,0],[Ay,h0],'r','linewidth',1.5)
plot([Bx,0],[By,h0],'r','linewidth',1.5)
plot([Cx,0],[Cy,h0],'r','linewidth',1.5)
plot(0,h0,'r','markersize',20)
hold off
else
s = @(h) AM(h) + BM(h) + CM(h)
subplot(2,1,1)
plot(h,s(h),'linewidth',1.5)
hold
plot(H,s(H),'r','markersize',20)
hold off

subplot(2,1,2)
plot([Ax,Bx],[Ay,By],'linewidth',1.5)
hold
plot([Ax,0],[Ay,h0],'r','linewidth',1.5)
plot([Bx,0],[By,h0],'r','linewidth',1.5)
plot(0,h0,'r','markersize',20)
hold off
end


%ex 4
clear
close all

%Initial
h = 4
L = 10
v = 1.05
Ax = 0
Ay = 0
Px = 0
Py = h
Bx = L
By = 0
Qy = 0
Qx = 0:L/100:L;

%Formular
PQ = @(Qx) sqrt(Py^2 + Qx.^2)
QB = @(Qx) Bx - Qx
t = @(Qx) PQ(Qx)./1 + QB(Qx)./v

Qx0 = sqrt((Py^2)/(v^2-1))

vLimit = sqrt((Py^2)/(10^2)+1)

if v > vLimit
  %Graph 
  subplot(2,1,1)
  plot(Qx,t(Qx),'linewidth',1.5)
  hold
  plot(Qx0,t(Qx0),'r','markersize',20)
  hold off
  subplot(2,1,2)
  plot([Ax,Bx],[Ay,By],'k','linewidth',1.5)
  hold
  plot([Px,Qx0,Bx],[Py,Qy,By],'linewidth',1.5)
  plot(Px,Py,'r','markersize',20)
  plot(Qx0,Qy,'r','markersize',20)
  plot(Bx,By,'r','markersize',20)
  ylim([-2,Py])
  hold off

else
  %Graph 
  subplot(2,1,1)
  plot(Qx,t(Qx),'linewidth',1.5)
  hold
  plot(10,t(10),'r','markersize',20)
  hold off
  subplot(2,1,2)
  plot([Ax,Bx],[Ay,By],'k','linewidth',1.5)
  hold
  plot([Px,Bx],[Py,By],'linewidth',1.5)
  plot(Px,Py,'r','markersize',20)
  plot(Bx,By,'r','markersize',20)
  ylim([-2,Py])
  hold off
end

%ex 5
clear
close all

%Initial
m = 2
R = 1
k = 3
M = 10
b = 5
w = 0:M/100:M;

%Function
A = @(w) (m*R*w.^2)./sqrt((k-M*w.^2).^2 + (b*w).^2)
%solve w0 when A is max
w0 = (sqrt(2)*k)/sqrt(2*M*k-b^2) 

if -w0 > 0
  %graph
  plot(w,A(w),'linewidth',1.5)
  grid
else
  %graph
  plot(w,A(w),'linewidth',1.5)
  hold
  plot(w0,A(w0),'r','markersize',20)
  grid
  hold off
end

%---------------------------------------------------
%Integral
%Integral
%ex 1 
clear
close all

%Initial 
a = 0.7
h = 4
x = 0:h/100:h

%function
%1
f = 4
xf = 0

adf = @(x) 4*x


%2
g = @(x) a*x.^2

adg = @(x) (a/3)*x.^3

%cross point when  g(x) = f(x)
x0 = sqrt(4/a)

%Area 
A = adf(x0) -adg(x0) - (adf(0)-adg(0))

%Centroid calculate through internet
xp = 1/A * -(x0^2*(7*x0^2-80))/40
yp = 1/A * -(x0*(49*x0^4-8000))/1000


%graph
plot([0,x0], [4,4],'linewidth',1.5)
hold
plot(x, g(x),'r','linewidth',1.5)
plot([0,0], [0,4],'k','linewidth',1.5)
plot(xp,yp,'k','markersize',20)
xlim([-2,x0])
ylim([-1,5])
grid
axis square
hold off


%ex2
clear
close all

%Initial
A = 5
T = 4
t = 0:T/100:4

%Formular

%u(t) area
au = 2* 1/2*(T/2*A)

%(u(t))^2 area
%half of the shape is the integral of
%   (A^2/4)*x^2 = A^2 => intergral (6.25/3)*x^3

au2 = 2*(((A^2)/4)/3)*(T/2)^3)

%solution
uavg = 1/T * au

urms = sqrt(1/T * au2)


%ex 3
clear
close all

%Initial
V1 = 500
p1 = 30
V2 = 50
p3 = 600
V4 = 100
V5 = 200
k1 = 1.2
k2 = 1.8

%Hint

p2 = p1*(V1/V2)^k1
V3 = V2
p4 = p3
p5 = p4*(V4/V5)
V6 = V1
p6 = p5*(V5/V6)^k2

%Var range
V = V2:(V1-V2)/1000:V1;
%Formular

%Lower curve
pLower = @(V) p1*(V1./V).^k1


%upper curve

%V3 - V4, range 50 from 50 - 100
V34 = V3:(V4-V3)/100:V4;
pUpper3to4 = 600


%V4 - V5, range 200 from 100 to 200
V45 = V4:(V5-V4)/100:V5;
pUpper4to5 = @(V45) p4*(V4./V45)


%V5 - V6, range300 from 200 to 500
V56 = V5:(V6-V5)/100:V6;
pUpper5to6 = @(V56) p5*(V5./V56).^k2


%Area 
%From V2 to V4
A24 = -((V1^k1*V2*V4^k1-V1^k1*V2^k1*V4)*p1+V4^k1*(600*V2^(k1+1)-600*V2^k1*V4)*k1+V4^k1*(600*V2^k1*V4-600*V2^(k1+1)))/(V2^k1*V4^k1*(k1-1))

%From V4 to V5
A45 = ((V5^k1*(V4^(k1+1)*log(V5)-V4^(k1+1)*log(V4))*k1+V5^k1*(V4^(k1+1)*log(V4)-V4^(k1+1)*log(V5)))*p4+(V1^k1*V4^k1*V5-V1^k1*V4*V5^k1)*p1)/(V4^k1*V5^k1*(k1-1))

%From V5 to V6
A56 = (V6^(-k2-k1)*(((V5^(k1+1)*V6^(k2+k1)-V5^(k2+k1)*V6^(k1+1))*k1-V5^(k1+1)*V6^(k2+k1)+V5^(k2+k1)*V6^(k1+1))*p5+(V6^k2*(V1^k1*V5^k1*V6-V1^k1*V5*V6^k1)*k2+V6^k2*(V1^k1*V5*V6^k1-V1^k1*V5^k1*V6))*p1))/(V5^k1*(k1-1)*(k2-1))
Area = A24 + A45 + A56

%Graph

plot(V,pLower(V),'r','linewidth',1.5)
hold
plot([V3,V4],[pUpper3to4,pUpper3to4],'linewidth',1.5)
plot([V3,V2],[pUpper3to4,p2],'linewidth',1.5)
plot([V1,V6],[p1,p6],'linewidth',1.5)
plot(V45,pUpper4to5(V45),'linewidth',1.5)
plot(V56,pUpper5to6(V56),'linewidth',1.5)
grid
hold off


%ex 4
clear
close all

%Initial
h = 3
L = 4
xMax = sqrt(L)
%angle of the right tri on the zy plane
gamma = atan(h/L)

%Define range
x= 0:xMax/100:xMax;
y = 0:L/100:L;
z = 0:h/100:h;

%Define function
ygraph =  @(x) x.^2
n = length(x);
z1 = zeros(1,n);
%Formular

%Cross section height and width
w = @(y) sqrt(y)
h = @(y) (L-y).*tan(gamma)

%Cross section area
A = @(y) w(y).*h(y)

%Volume formular

V = (4*L^(5/2)*tan(gamma))/15

%graph
plot3(x,ygraph(x),z1,'b','linewidth',1.5)
hold
plot3(x,ygraph(x),h(ygraph(x)),'b','linewidth',1.5)
plot3([0,0,0,2],[0,0,4,4],[0,3,0,0],'b','linewidth',1.5)
plot3([0,0],[0,L],[0,0],'b','linewidth',1.5)
%cross section at x = 1
plot3([0,0],[ygraph(1),ygraph(1)],[0,h(ygraph(1))],'r','linewidth',1.5)
plot3([1,1],[ygraph(1),ygraph(1)],[0,h(ygraph(1))],'r','linewidth',1.5)
plot3([0,1],[ygraph(1),ygraph(1)],[0,0],'r','linewidth',1.5)
plot3([0,1],[ygraph(1),ygraph(1)],[h(ygraph(1)),h(ygraph(1))],'r','linewidth',1.5)
hold off
xlabel = 'x'
ylabel = 'y'
zlabel = 'z'
axis equal 

%ex 5
clear
close all

%Initial
v0 = 10
alfa = 30
h = 5
g = 9.81
a = (-g/(2*(v0*cosd(alfa))^2))
b = tand(alfa)

%Formular
xMax = ((-b-sqrt((b)^2-4*a*h))/(2*a))
dx = xMax/100
x = 0:dx:xMax;
y = @(x)   a* x.^2  + b*x + h
y1 =  a* x.^2  + b*x + h;
%length
s = (asinh(b+2*xMax*a)-asinh(b)+(b+2*xMax*a)*sqrt(b^2+4*xMax*a*b+4*xMax^2*a^2+1)-b*sqrt(b^2+1))/(4*a)

%2nd methods
n = length(x)
dy1=y1(2:n)-y1(1:n-1);%vertical distances of consecutive points
ds=sqrt(dx.^2+dy1.^2); %distances of consecutive points
sn=sum(ds) %length of the graph

%Graph

plot(x,y(x),'b','linewidth',1.5)
hold
plot([-1,xMax],[0,0],'k','linewidth',1.5)
plot([0,0],[-2,8],'k','linewidth',1.5)
hold off
ylim([-2,8])
xlim([-1,xMax])
grid
axis equal

%ex 6
clear
close all

%Initial
R = 5
L = 3
t = 0:2*pi/100:2*pi;

%equation
x = @(t)  R*cos(t) + L*cos(3*t)
y = @(t) R*sin(t) + L*sin(3*t)

x1 =   R*cos(t) + L*cos(3*t);
y1 =  R*sin(t) + L*sin(3*t);

n =length(t)

dx = x1(2:n) - x1(1:n-1);
dy = y1(2:n) - y1(1:n-1);
ds=sqrt(dx.^2+dy.^2); %distances of consecutive points
sn=sum(ds)

plot(x(t),y(t),'b','linewidth',1,5)
grid
axis equal

%-------------------------------------------------
%PVA
pkg install -forge optim
pkg load optim
%PVA ex 1
clear
close all
%when moving
s = @(t) sqrt((5 + 2*t).^2-4.^2)
v = @(t) (4*t + 10)./sqrt((2*t +5).^2-16)
a = @(t) -64./((2*t+5).^2-16).^(3/2)
Tstart = 0
Tend = 2
t = Tstart:(Tend-Tstart)/100:Tend;

figure(1)
subplot(311)
plot(t,s(t),'linewidth',1.5)
grid
xlim([Tstart,Tend])
title('position')
title('position')
subplot(312)
plot(t,v(t),'linewidth',1.5)
grid
xlim([Tstart,Tend])
title('velocity')
subplot(313)
plot(t,a(t),'linewidth',1.5)
grid
xlim([Tstart,Tend])
title('acceleration')
xlabel('time t')

%PVA ex 2
clear
close all
%initial
h = 5
th0 =  pi/6
w0 = 2*pi
alfa = -5*pi

% pos function
x = h*tan(th0)

dx = (x^2 + h^2)*w0/h

d2x = ((x^2+h^2)^2*alfa+h*(2*x*(dx)^2))/(h*(x^2+h^2))

%PVA ex 3
clear
close all
%Initial values
BC = 50
OA = sqrt(100^2+200^2)
OD = 600
DE = 400
OEy = 400
AB = sqrt(50^2+ 300^2)
AD = sqrt(400^2+100^2)
w = 2*pi %random value since the task was not given any
T=2*pi/w %period
dt=T/1000
t=0:dt:T;

% FIxed Coordinate
OCx = 200
OCy = 100

%alfa
alfa = w*t

%B coordinate
Bx = OCx + BC*sin(alfa)
By = OCy + BC*cos(alfa)
OB = sqrt(By.^2 + Bx.^2)

%Beta angle
Beta = asin(By./OB)

%calculate angle AOB and AOD from Cosine formular

AOB = acos((OA.^2+OB.^2-AB.^2)./(2*OA.*OB))

AOD = acos((OA.^2+OD.^2-AD.^2)./(2*OA.*OD))

%polar angle

gamma =  AOB- AOD + Beta

%Dy and Dx
Dy = OD*sin(gamma)
Dx = OD*cos(gamma)

%DEy 
DEy = Dy - OEy

%DEx
DEx = sqrt(400^2-(DEy).^2)

n = length(t)

s = DEx + OD*cos(gamma) %OEx

vk = (s(2:n)-s(1:n-1))/dt; % 2:n second charactem of the range

ak = (vk(2:n-1)-vk(1:n-2))/dt;



figure(1)
subplot(3,1,1)
plot(t,s/1000,'linewidth',1.5)
grid
title('position s(t)')

subplot(3,1,2)
plot(t(1:n-1),vk/1000,'linewidth',1.5)
grid
title('velocity v(t)')

subplot(3,1,3)
plot(t(1:n-2),ak/1000,'linewidth',1.5)
grid
title('acceleration a(t)')
xlabel('time t')

%Finding A
OAy = OA * sin(Beta + gamma)
OAx = OA * cos(Beta + gamma)

%O roots
Ox = 0
Oy = 0
%% animation
figure(2)
for k=1:5:n
plot([Ox,OAx(k),Bx(k),OCx],[Oy,OAy(k),By(k),OCy],'b.-','linewidth',1.5,'markersize',15)
hold
plot([OAx(k),Dx(k),s(k)],[OAy(k),Dy(k),OEy],'r.-','linewidth',1.5,'markersize',15)
hold off
grid
axis equal
%axis ([-4,4,-4,4])
pause(0.01)
end




%PVA ex5
clear
close all

%Initial data
Ox = 0
Oy = 0
OB = 250
OC = 600
BD = 900
w = 2*pi
T=2*pi/w %period
dt=T/1000
t=0:dt:T;
alfa = w*t


%Finding B coordinates
OBx = -OB*cos(alfa)
OBy = OB*sin(alfa)

%Finding C coordinates
OCx = -600
OCy = 0

%Finding CB
CB = sqrt((OBy).^2+(OCx-OBx).^2)

%Find CD
CD = BD - CB

%Find polar angle atan2 give the x y motion reference from the O
gamma = atan2((OCy - OBy),(OCx -OBx))


%Find D
ODx = BD*cos(gamma) + OBx
ODy = BD*sin(gamma) + OBy


%velocity and acceleration
n = length(t)
vkx = (ODx(2:n)-ODx(1:n-1))/dt;
vky = (ODy(2:n)-ODy(1:n-1))/dt;
vk = sqrt((vkx).^2 + (vky).^2);


akx = (vkx(2:n-1)-vkx(1:n-2))/dt;
aky = (vky(2:n-1)-vky(1:n-2))/dt;
ak = sqrt((akx).^2 + (aky).^2);


%aT and aN
t0 = 0.2
aT= @(t) (vkx(t).*akx(t)+vky(t).*aky(t))./vk(t)
aN= @(t) (vkx(t).*aky(t)-vky(t).*akx(t))./vk(t)


%graph
figure(1)
subplot(2,1,1)
plot(t(1:n-2),ak*10^4,'linewidth',1.5)
grid
title('acceleration a(t)')

subplot(2,1,2)
plot(t(1:n-1),vk,'linewidth',1.5)
grid
title('velocity v(t)')

figure(2)
vkmax = max(vk)
k = 200
p= 2/20
plot([Ox,OBx(k),ODx(k)],[Oy,OBy(k),ODy(k)],'linewidth',1.5,'markersize',15)
axis equal
hold
p2=plot([ODx(k),ODx(k)+p*akx(k)],[ODy(k),ODy(k)+p*aky(k)],'g','linewidth',2)
p3=plot([ODx(k),ODx(k)+p*aT(k)*vkx(k)/vk(k)],[ODy(k),ODy(k)+p*aT(k)*vky(k)/vk(k)],'c','linewidth',3)
p1=plot([ODx(k),ODx(k)+p*vkx(k)],[ODy(k),ODy(k)+p*vky(k)],'r','linewidth',2)
p4=plot([ODx(k),ODx(k)+p*aN(k)*(-vky(k))/vk(k)],[ODy(k),ODy(k)+p*aN(k)*vkx(k)/vk(k)],'m','linewidth',2)
p5=plot(ODx,ODy,'b.','linewidth',1.5)
plot(ODx(k),ODy(k),'b.','markersize',15)
hold off


%ex 6
clear
close all
% a
%Initial data 
R = 3
OM = R
ON = R
L = 5
MP = L
NP = L
Ax = 1
Ay = 6
Bx = 6
By = 4

T = 2
dt=T/1000
t=0:dt:T;
n=length(t)

%PVA
s=6*(t/T).^5-15*(t/T).^4+10*(t/T).^3; %s(t)
ds=(30*(t/T).^4-60*(t/T).^3+30*(t/T).^2)*1/T; %s'(t)
d2s=(120*(t/T).^3-180*(t/T).^2+60*(t/T))*(1/T)^2; %s''(t)

%position of P
x=Ax+s*(Bx-Ax);
y=Ay+s*(By-Ay);

vx=ds*(Bx-Ax);
vy=ds*(By-Ay);
ax=d2s*(Bx-Ax);
ay=d2s*(By-Ay);

v=sqrt(vx.^2+vy.^2); 
a=sqrt(ax.^2+ay.^2); 

%% inverse kinematics: x,y -> alfa ja beta

OP=sqrt(x.^2+y.^2);
phi=atan(y./x);
kO=acos((OM^2+OP.^2-MP^2)./(2*OM*OP));
kM=acos((OM^2+MP^2-OP.^2)./(2*OM*MP));
alfa=phi-kO;
beta=alfa + 2*kO;

%angular velocities and accelerations numerically

valfa=(alfa(2:n)-alfa(1:n-1))/dt; 
vbeta=(beta(2:n)-beta(1:n-1))/dt; 
aalfa=(valfa(2:n-1)-valfa(1:n-2))/dt; 
abeta=(vbeta(2:n-1)-vbeta(1:n-2))/dt; 

figure(2)
subplot(3,1,1)
plot(t,alfa,'r','linewidth',1.5)
hold on
plot(t,beta,'b','linewidth',1.5)
hold off
grid
title('angles \alpha and \beta (rad)')
legend({'\alpha','\beta'},'fontsize',11)


subplot(3,1,2)
plot(t(1:n-1),valfa,'r','linewidth',1.5)
hold on 
plot(t(1:n-1),vbeta,'b','linewidth',1.5)
hold off
grid
title('angular velocities (rad/sec)')
legend({'\alpha','\beta'},'fontsize',11)


subplot(3,1,3)
plot(t(1:n-2),aalfa,'r','linewidth',1.5)
hold on 
plot(t(1:n-2),abeta,'b','linewidth',1.5)
hold off
grid
title('angular accelerations (rad/sec^2)')
legend({'\alpha','\beta'},'fontsize',11)
xlabel('time t (sec)')

%% animation

%direct kinematics: alfa,beta -> Px,Py

Mx=OM.*cos(alfa);
My=OM.*sin(alfa);

Nx=ON.*cos(alfa+kO.*2);
Ny=ON.*sin(alfa+kO.*2); 

Px=OP.*cos(alfa+kO);
Py=OP.*sin(alfa+kO);

R=OM+MP+1 %limit for axis

figure(3)
for k=1:2:n
  plot(x(k),y(k),'gs','markersize',12) 
  hold on 
  plot([0,Mx(k)],[0,My(k)],'b-o','linewidth',2)
  plot([0,Nx(k)],[0,Ny(k)],'b-o','linewidth',2)
  plot([Mx(k),Px(k)],[My(k),Py(k)],'r-o','linewidth',2) 
  plot([Nx(k),Px(k)],[Ny(k),Py(k)],'r-o','linewidth',2) 
  plot(x,y,'k','linewidth',2)
  plot(Px(k),Py(k),'rs','linewidth',2) 
  hold off
  grid 
  axis([-1,R,-1,R])
  axis square 
  pause(0.001)
end 

clear
close all
% b
%Initial data 
R = 3
OM = R
ON = R
L = 5
MP = L
NP = L
Ax = 1
Ay = 6
Bx = 6
By = 4
T = 2

dt=T/1000
t=0:dt:T;
n=length(t)

s=6*(t/T).^5-15*(t/T).^4+10*(t/T).^3; %s(t)

% when t = 0 ==> OP = OA
x=Ax
y=Ay
OP=sqrt(x.^2+y.^2);
phi=atan(y./x);
kO=acos((OM^2+OP.^2-MP^2)./(2*OM*OP));
kM=acos((OM^2+MP^2-OP.^2)./(2*OM*MP));
alfaA=phi-kO;
betaA=alfaA + 2*kO;

% when t = T ==> OP = OB
x=Bx
y=By
OP=sqrt(x.^2+y.^2);
phi=atan(y./x);
kO=acos((OM^2+OP.^2-MP^2)./(2*OM*OP));
kM=acos((OM^2+MP^2-OP.^2)./(2*OM*MP));
alfaB=phi-kO;
betaB=alfaB + 2*kO;

alfa=alfaA+s*(alfaB-alfaA);
beta=betaA+s*(betaB-betaA);

figure(1)
plot(t,alfa,'r','linewidth',1.5)
hold 
plot(t,beta,'b','linewidth',1.5)
hold off 
grid
xlabel('time t')
legend({'\alpha','\beta'},'fontsize',12)

%direct kinematics: alfa,beta -> Px,Py

Mx=OM.*cos(alfa);
My=OM.*sin(alfa);

Nx=ON.*cos(alfa+kO.*2);
Ny=ON.*sin(alfa+kO.*2); 

Px=OP.*cos(alfa+kO);
Py=OP.*sin(alfa+kO);

R=OM+MP+1 %limit for axis

vx=(Px(2:n)-Px(1:n-1))/dt;
vy=(Py(2:n)-Py(1:n-1))/dt;
v=sqrt(vx.^2+vy.^2);
ax=(vx(2:n-1)-vx(1:n-2))/dt;
ay=(vy(2:n-1)-vy(1:n-2))/dt;
a=sqrt(ax.^2+ay.^2);

figure(2)
subplot(2,1,1)
plot(t(1:n-1),v,'linewidth',1.5)
grid
title('P:n vauhti v')
subplot(2,1,2)
plot(t(1:n-2),a,'linewidth',1.5)
grid
title('P:n kiihtyvyys a')
xlabel('aika t')

figure(3)
k = 500
plot(Ax,Ay,'r.',Bx,By,'g.','markersize',20)
hold on
plot([0,Mx(k)],[0,My(k)],'b-o','linewidth',2)
plot([0,Nx(k)],[0,Ny(k)],'b-o','linewidth',2)
plot([Mx(k),Px(k)],[My(k),Py(k)],'r-o','linewidth',2) 
plot([Nx(k),Px(k)],[Ny(k),Py(k)],'r-o','linewidth',2) 
plot(x,y,'k','linewidth',2)
plot(Px,Py,'k','linewidth',1.5)
plot([-1,R],[0,0],'k')
plot([0,0],[-1,R],'k')
hold off
grid
axis([-1,R,-1,R])
axis square
legend('A','B')
%%


%ex 4 PVA
clear
close all
%Initial values
r = 2
V = 3
T = 2*pi*r
dt = T/1000
t = 0:dt:T;
T0 = 4.1888
t0 = 0.3*T0


%Formular
%circle formular
yC = @(t) 2 + r*sin(V*t./r)
xC = @(t) V*t0 + r*cos(V*t./r)
%P
x = @(t) r*((V*t./r) -sin(V*t./r))
y = @(t) r*(1 -cos(V*t./r))

%V
vx = @(t) V-V*cos((V*t)./r)
vy = @(t)  V*sin((V*t)./r)
v = @(t) sqrt(vx(t).^2 + vy(t).^2)

%A
ax = @(t) (V^2*sin((V*t)./r))./r
ay = @(t) (V^2*cos((V*t)./r))./r
a = @(t) sqrt(ax(t).^2 + ay(t).^2)

%aT & aN
aT= @(t)(vx(t).*ax(t)+vy(t).*ay(t))./v(t)
aTx= @(t) aT(t).*vx(t)./v(t)
aTy=@(t) aT(t).*vy(t)./v(t)
aN= @(t)(vx(t).*ay(t)-vy(t).*ax(t))./v(t)
aNx= @(t) aN(t).*(-vy(t))./v(t)
aNy=@(t) aN(t).*vx(t)./v(t)

%graph
p = 5 % scale factor
figure(1)
plot(x(t),y(t),'b','linewidth',1.5)
hold
plot([0,T],[0,0],'k','linewidth',1.5)
plot([x(t0),x(t0)+p*vx(t0)/v(t0)],[y(t0),y(t0)+p*vy(t0)/v(t0)],'r','linewidth',3)
plot([x(t0),x(t0)+p*aTx(t0)/a(t0)],[y(t0),y(t0)+p*aTy(t0)/a(t0)],'c','linewidth',3)
plot([x(t0),x(t0)+p*aNx(t0)/a(t0)],[y(t0),y(t0)+p*aNy(t0)/a(t0)],'m','linewidth',3)
plot([x(t0),x(t0)+p*ax(t0)/a(t0)],[y(t0),y(t0)+p*ay(t0)/a(t0)],'g','linewidth',3)
plot(xC(t), yC(t), 'r','linewidth',1.5);
plot(V*t0, 2, 'r','markersize',15);
hold off
grid
axis([0,T,-1,7])

figure(2)
subplot(2,1,1)
plot(t,v(t),'linewidth',1.5)
axis([0,2*r*pi/V,0,6])
grid
title('P:n vauhti v')
subplot(2,1,2)
plot(t,a(t),'linewidth',1.5)
axis([0,2*r*pi/V,0,6])
grid
title('P:n kiihtyvyys a')
xlabel('aika t')

figure(3)
subplot(2,1,1)
plot(t,aT(t),'linewidth',1.5)
axis([0,2*r*pi/V,-5,5])
grid
title('P:n tangentiaalikiihtyvyys aT')
subplot(2,1,2)
plot(t,aN(t),'linewidth',1.5)
grid
axis([0,2*r*pi/V,-5,0])
title('P:n normaalikiihtyvyys aN')
xlabel('aika t')

%-------------------------------------------------------------------------
%Taylor polynominal
%ex1: log(1+x)
clear
close all
%degree n
n = 5
%graphs for x=a...b
a=-1
b=1
x=a:(b-a)/100:b;
f= @(x) log(1+x)
x0=0

Tn=0
%factorial(k)=k!
for k=1:n
    Tn=Tn+(-1)^(k-1)/k * x.^k; 
end

plot(x,f(x),'b','linewidth',1.5)
hold
plot(x,Tn,'r','linewidth',1.2)
hold off
grid
xlabel('x')
legend({'log(1+x)','T_n(x)'},'fontsize',12)
title(['degree n = ',num2str(n)])
ylim([-2,2])















