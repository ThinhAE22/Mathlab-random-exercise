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
