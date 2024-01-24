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
axis([0,(e/2)*pi,0,6])
grid
title('P:n vauhti v')
subplot(2,1,2)
plot(t,a(t),'linewidth',1.5)
axis([0,(e/2)*pi,0,6])
grid
title('P:n kiihtyvyys a')
xlabel('aika t')

figure(3)
subplot(2,1,1)
plot(t,aT(t),'linewidth',1.5)
axis([0,(e/2)*pi,-5,5])
grid
title('P:n tangentiaalikiihtyvyys aT')
subplot(2,1,2)
plot(t,aN(t),'linewidth',1.5)
grid
axis([0,(e/2)*pi,-5,0])
title('P:n normaalikiihtyvyys aN')
xlabel('aika t')











































