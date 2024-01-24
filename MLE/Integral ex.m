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

%ex 7
clear 
close all



