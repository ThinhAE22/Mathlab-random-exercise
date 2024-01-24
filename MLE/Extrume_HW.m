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







