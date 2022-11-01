%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: 3Lemon2Orange (GitHub)
% About Explicit and Implicit scheme in 1D heat equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
format compact

starttime=cputime;

test=2; % 1-Explicit, 2-Inplicit

xmin=0;
xmax=1;
tstart=0;
tfinal=1;
N=20;
dx=(xmax-xmin)/N;
dt=0.01;

n=N+1;
m=(tfinal-tstart)/dt +1; 
c=dt/(dx^2); %gamma
x=[xmin:dx:xmax];
t=[tstart:dt:tfinal];
u=zeros(n,m);

%initial value
u(1,:)=zeros(1,m);
u(n,:)=zeros(1,m);

for i=2:N
    temp=(i-1)*dx;
    if temp <= 0.5
        u(i,1)=2*temp;
    else
        u(i,1)=2 - 2*temp;
    end
end

%Calculation
if test==1   %Explicit
    for k=1:m-1
        for j=2:N
            u(j,k+1)=(1-2*c)*u(j,k)+c*(u(j+1,k)+u(j-1,k));
        end
    end
end

if test==2   %Inplicit
    a=(1+2*c)*ones(N-1,1);
    b=-c*ones(N-1,1);
    A=spdiags([b a b],-1:1,N-1,N-1);
    A=full(A);
    for k=1:m-1
        v=A\u(2:n-1,k);
        u(2:n-1,k+1)=v;
    end
end

totalrunningtime=floor((cputime-starttime));
sprintf('Total running time is %d second(s)',totalrunningtime)

error=zeros(n,m);
if test==1
    for k=1:m-1
        for j=2:N
            error(j,k)=(u(j,k+1)-u(j,k))/dt - (u(j+1,k)-2*u(j,k)+u(j-1,k))/(dx)^2;
        end
    end
end

if test==2
    for k=1:m-1
        for j=2:N
            error(j,k+1)=(u(j,k+1)-u(j,k))/dt - (u(j+1,k+1)-2*u(j,k+1)+u(j-1,k+1))/(dx)^2;
        end
    end
end

if test==1
    figure(1)
    plot(x,u(:,1),x,u(:,2),x,u(:,3),x,u(:,4))
end

%plot mesh
figure(2)
mesh(x,t,u')
xlabel('x')
ylabel('t')
zlabel('u')

%plot different layers
figure(3)
for l=[1:5:m]
    hold on
    plot(x,u(:,l))
end
hold off
xlabel('x')
ylabel('t')

