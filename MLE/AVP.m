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





