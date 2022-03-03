clear all; clc; close(gcf); 
L=pi; % x in [0, L] 
H=3; % y in [0, H] 
N=101; % 
M=51; % 
T=pi; % 
hx=L/(N-2); % 
hy=H/(M-2); % 
%J=2+round(2*T*sqrt(1/hx^2+1/hy^2)); 
J=101; 
tau=T/(J-1); % 
x=zeros(N,1); % 
y=zeros(M,1); % 
t=zeros(J,1); % 
% --------------------------------------------------- 
for n=1:N 
    x(n)=(n-1)*hx-hx/2; 
end 
for m=1:M 
    y(m)=(m-1)*hy-hy/2; 
end 
for j=1:J 
    t(j)=(j-1)*tau; 
end 
%------------------------------------------------------------------ 
% ------------------------------------------- 
u=zeros(M,N,J); 
for j=1:J 
    for n=1:N 
        for m=1:M 
            u(m,n,j)=(exp(-4*t(j))+(sin(t(j))-cos(t(j))+exp(-4*t(j)))/2)*cos(x(n)); 
        end 
    end 
end 
%------------------------------------------------------------------ 
% ------------------------------------------------ 
v=zeros(M,N,J); % 
err=zeros(M,N,J); % 
% 
for n=1:N 
    for m=1:M 
        v(m,n,1)=cos(x(n)); 
    end 
end
w=zeros(M,N); %
alpha_x=zeros(N-1,1); 
beta_x=zeros(N-1,1); 
alpha_y=zeros(M-1,1); 
beta_y=zeros(M-1,1); 
Ax=1; 
Cx=2*(1+4*hx^2/tau); 
Ay=1; 
Cy=2*(1+4*hy^2/tau); 
for j=1:J-1 
    % -------------- 
    for m=2:M-1 
        alpha_x(1)=1;
        % ---------------------------- 
        for n=2:N-1 
            F=(v(m-1,n,j)-2*v(m,n,j)+v(m+1,n,j))*hx^2/hy^2+2*4*hx^2/tau*v(m,n,j)+4*hx^2*cos(x(n))*sin(t(j)+tau/2); 
            alpha_x(n)=Ax/(Cx-Ax*alpha_x(n-1)); 
            beta_x(n)=(F+Ax*beta_x(n-1))/(Cx-Ax*alpha_x(n-1)); 
        end 
        %------------------------------ 
         w(m,N)=beta_x(N-1)/(1-alpha_x(N-1)); 
        for n=N-1:-1:1 
            w(m,n)=alpha_x(n)*w(m,n+1)+beta_x(n); 
        end 
%         w(m,N)= w(m,N-1);
%         w(m,1) = w(m,2);
    end 
        % --------- 
        alpha_y(1)=1;
    for n=2:N-1 
        % -------------------------- 
        for m=2:M-1 
            F=(w(m,n-1)-2*w(m,n)+w(m,n+1))*hy^2/hx^2+2*4*hy^2/tau*w(m,n)+4*hy^2*cos(x(n))*sin(t(j)+tau/2); 
            alpha_y(m)=Ay/(Cy-Ay*alpha_y(m-1)); 
            beta_y(m)=(F+Ay*beta_y(m-1))/(Cy-Ay*alpha_y(m-1)); 
        end 
        % ------------------------------ 
        v(M,n,j+1)=beta_y(M-1)/(1-alpha_y(M-1)); 
        err(M,n,j+1)=v(M,n,j+1)-u(M,n,j+1); 
        for m=M-1:-1:1 
            v(m,n,j+1)=alpha_y(m)*v(m+1,n,j+1)+beta_y(m); 
            err(m,n,j+1)=v(m,n,j+1)-u(m,n,j+1); 
        end 
    end
    %
    for m=1:M
        v(m,1,j+1) = v(m,2,j+1);
        err(m,1,j+1)=v(m,1,j+1)-u(m,1,j+1); 
        v(m,N,j+1) = v(m,N-1,j+1); 
        err(m,N,j+1)=v(m,N,j+1)-u(m,N,j+1); 
    end 
    for n=1:N
        v(1,n,j+1) = v(2,n,j+1);
        err(1,n,j+1)=v(1,n,j+1)-u(1,n,j+1);
        v(M,n,j+1) = v(M-1,n,j+1); 
        err(M,n,j+1)=v(M,n,j+1)-u(M,n,j+1); 
    end
end
j=20; 
figure 
surf(x,y,u(:,:,j)) 
xlabel('x')
ylabel('y') 
zlabel('u') 
title('analitical solution') 
figure 
surf(x,y,v(:,:,j)) 
xlabel('x') 
ylabel('y') 
zlabel('u') 
title('calculated solution') 
figure 
surf(x,y,err(:,:,j)) 
xlabel('x') 
ylabel('y') 
zlabel('u') 
title('error')
