V=10;
R=50;
L=R*pi;
N=160;
w=10;
dt0=L/((N-1)*V);
x0=zeros(3*N+1,1);
i=1;
while i<=N
    x0(i)=-w/2;
    i=i+1;
end
i=N+1;
while i<=2*N
    x0(i)=(i-N-1)*V*dt0;
    i=i+1;
end
x0(2*N+1)=dt0;
i=2*N+2;
while i<3*N+1
    x0(i+1)=x0(i)+V*dt0/R;
    i=i+1;
end
xlb=zeros(3*N+1,1);
xub=zeros(3*N+1,1);
i=1;
while i<N+1
    xlb(i)=-w/2;
    xub(i)=w/2;
    i=i+1;
end
i=N+1;
while i<2*N+1
    xlb(i)=0;
    xub(i)=L;
    i=i+1;
end
xlb(2*N+1)=0;
xub(2*N+1)=dt0+2*w/((N-1)*V);
i=2*N+2;
while i<3*N+2
    xlb(i)=-pi/2;
    xub(i)=pi/2;
    i=i+1;
end

options = optimset('Display','iter','TolX',1e-5,'TolFun', 1e-5, 'TolCon', 1e-5,...
        'MaxFunEval', 100000,'MaxIter', 1000);
[x_star,f]=fmincon(@(x)objfun_(x,N),x0,[],[],[],[],xlb,xub,@(x)confun_(x,N,V,R,L),options);

figure
plot(x_star(1:N))
title('s_n vs t')

figure
plot(x_star(N+1:2*N))
title('s_t vs t')

figure
stairs(x_star(2*N+2:3*N+1))
title('\theta vs t')

Tf=(N-1)*f
Tf0=(N-1)*dt0

x_cord=zeros(N);
x_out=zeros(N);
x_in=zeros(N);
y_cord=zeros(N);
y_out=zeros(N);
y_in=zeros(N);

i=1;
while i<N+1
   x_cord(i)=(-x_star(i)+R)*cos(x_star(i+N)/R);
   y_cord(i)=(-x_star(i)+R)*sin(x_star(i+N)/R);
   x_in(i)=(-w/2+R)*cos(x_star(i+N)/R);
   x_out(i)=(w/2+R)*cos(x_star(i+N)/R);
   y_in(i)=(-w/2+R)*sin(x_star(i+N)/R);
   y_out(i)=(w/2+R)*sin(x_star(i+N)/R);
    i=i+1;
end

figure
plot(x_in,y_in,'b')
hold on
plot(x_cord,y_cord,'r')
hold on
plot(x_out,y_out,'b')
hold off
axis([-R-w R+w 0 2*(R+w)])
title('Actual trajectory')