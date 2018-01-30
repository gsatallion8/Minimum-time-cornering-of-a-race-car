function [c,ceq] = confun_( x,N,V0,R,L,g )
ceq=zeros(3*N+1,1);
c=zeros(4*N,1);
% c=[];
i=1;
while i<N
    ceq(i)=x(i+1)-x(i)-x(3*N+1+i)*sin(x(2*N+1+i))*x(2*N+1);
    ceq(i+N-1)=x(i+N+1)-x(i+N)-x(3*N+1+i)*cos(x(2*N+1+i))*x(2*N+1)/(1-x(i)/R);
    ceq(i+2*N-2)=x(3*N+2+i)-x(3*N+1+i)-x(4*N+1+i)*x(2*N+1);
    i=i+1;
end
ceq(3*N-2)=x(1);
ceq(3*N-1)=x(N+1);
ceq(3*N)=x(2*N)-L;
ceq(3*N+1)=x(3*N+2)-V0;
% ceq(2*N+2)=x(2*N+2);
%Additional constraints
i=2*N+3;
while i<=3*N+1
    c(i-2*N-1)=(x(i)-x(i-1))/x(2*N+1)-0.3;
    c(i-N-1)=-(x(i)-x(i-1))/x(2*N+1)-0.3;
    c(i-1)=(x(2*N+i)-x(2*N+i-1))/x(2*N+1)-5*g;
    c(i+N-1)=-(x(2*N+i)-x(2*N+i-1))/x(2*N+1)-5*g;
    i=i+1;
end