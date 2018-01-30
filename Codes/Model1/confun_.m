function [c,ceq] = confun_( x,N,V,R,L )
ceq=zeros(3*N+1,1);
c=zeros(2*N,1);
% c=[];
i=1;
while i<N
    ceq(i)=x(i+1)-x(i)-V*sin(x(2*N+1+i))*x(2*N+1);
    ceq(i+N-1)=x(i+N+1)-x(i+N)-V*cos(x(2*N+1+i))*x(2*N+1)/(1-x(i)/R);
    i=i+1;
end
ceq(3*N-1)=x(1);
ceq(3*N)=x(N+1);
ceq(3*N+1)=x(2*N)-L;
% ceq(2*N+2)=x(2*N+2);
%Additional constraints
i=2*N+2;
while i<=3*N+1
    c(i-2*N-1)=(x(i)-x(i-1))/x(2*N+1)-0.3;
    c(i-N-1)=-(x(i)-x(i-1))/x(2*N+1)-0.3;
    i=i+1;
end