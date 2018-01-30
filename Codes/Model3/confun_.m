function [c,ceq] = confun_( x,N,V0,R,L,g,M,I,r,Cd,B,C,D,E,Pmax )
ceq=zeros(6*N+4,1);
c=zeros(5*N,1);
% c=[];
i=1;
while i<N
    ceq(i)=x(i+1)-x(i)-x(3*N+1+i)*sin(x(6*N+1+i)+x(5*N+1+i))*x(2*N+1);
    ceq(i+N-1)=x(i+N+1)-x(i+N)-x(3*N+1+i)*cos(x(6*N+1+i)+x(5*N+1+i))*x(2*N+1)/(1-x(i)/R);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     k=x(4*N+1+i)/(x(3*N+1+i)*cos(x(5*N+1+i)))-1;
%     lmda=x(5*N+1+i)-x(2*N+1+i);
%     sgmx=k/(1+k);
%     sgmy=tan(lmda)/(1+k);
%     sgm=sqrt(sgmx^2+sgmy^2);
%     S=sgmx*D*sin(C*atan(B*sgm-E*(B*sgm-atan(B*sgm))))/sgm;
%     f=S*sgmy/sgmx;
%     alpha=(x(7*N+1+i)-S*r)*M*g/I;

% either the above or the below block should be left uncommented

    S=x(7*N+1+i)/r;
    f=0;
    alpha=0;
    
    a=g*(S*cos(x(2*N+1+i)-x(5*N+1+i))-f*sin(x(2*N+1+i)-x(5*N+1+i)))-Cd*x(3*N+1+i)^2*cos(x(5*N+1+i))^3;
    b=g*(S*sin(x(2*N+1+i)-x(5*N+1+i))+f*cos(x(2*N+1+i)-x(5*N+1+i)))/x(3*N+1+i)-Cd*x(3*N+1+i)*sin(x(5*N+1+i))*cos(x(5*N+1+i))^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ceq(i+2*N-2)=x(3*N+2+i)-x(3*N+1+i)-a*x(2*N+1);
    ceq(i+3*N-3)=x(4*N+2+i)-x(4*N+1+i)-alpha*x(2*N+1);
    ceq(i+4*N-4)=x(5*N+2+i)-x(5*N+1+i)-b*x(2*N+1);
    ceq(i+5*N-5)=x(6*N+2+i)-x(6*N+1+i)+x(3*N+1+i)*cos(x(6*N+1+i)*x(5*N+1+i))/(R-x(i));
    i=i+1;
end
ceq(5*N-4)=x(1);
ceq(5*N-3)=x(N+1);
ceq(5*N-2)=x(2*N)-L;
ceq(5*N-1)=x(3*N+2)-V0;
ceq(5*N)=x(5*N+2);
ceq(5*N+1)=x(6*N+1);
ceq(5*N+2)=x(6*N+2);
ceq(5*N+3)=x(7*N+1);
ceq(6*N+4)=x(4*N+1)-V0/r;
%Additional constraints
i=2*N+3;
while i<=3*N+1
    c(i-2*N-1)=(x(i)-x(i-1))/x(2*N+1)-0.3;
    c(i-N-1)=-(x(i)-x(i-1))/x(2*N+1)-0.3;
    c(i-1)=x(5*N+i)*x(2*N+i)*M*g-Pmax;
    c(i+N-1)=(x(5*N+i)-x(5*N+i-1))/x(2*N+1)-Pmax/100;
    c(i+2*N-1)=-(x(5*N+i)-x(5*N+i-1))/x(2*N+1)-Pmax/100;
    i=i+1;
end