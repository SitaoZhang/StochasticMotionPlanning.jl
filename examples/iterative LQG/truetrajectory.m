clc;
load u.mat
N = 180;
x1t=zeros(1,N+1);
x2t=zeros(1,N+1);
x1t(1,1)=-10;
x2t(1,1)=10;
dt=0.01;
for i=1:N
    r1=(randn(1,1));
    x1t(i+1)=x1t(i)+dt*((-x1t(i)+2*x2t(i))+0.5*(u(i)+2*sin(2*i*dt)-2*sin(2*(i-1)*dt)))+sqrt(dt)*0.5*(0.01*(x2t(i)^2)+0.01*(x1t(i)^2))*r1;
    x2t(i+1)=x2t(i)+dt*((-3*x1t(i)-x2t(i)-0.05*x2t(i)^3)-2*(u(i)+2*sin(2*i*dt)-2*sin(2*(i-1)*dt)))+sqrt(dt)*(-2)*(0.01*(x2t(i)^2)+0.01*(x1t(i)^2))*r1;
end
xt=[x1t;x2t];
save xt.mat xt

load x.mat
r = rectangle('Position',[-8 0 5 5.8]);
r.EdgeColor = 'k';
r.FaceColor = [0.99 0.99 0.99];
r.LineStyle = "--";

hold on
plot(x(1, :), x(2, :))
plot(x1t,x2t)

legend({'ideal sample path' 'true sample path'},'Location','northeast')
