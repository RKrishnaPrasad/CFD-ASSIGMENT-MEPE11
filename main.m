close all;
clear all; clc;
%dU/dt= k*d/dn(|dU/dn|^(n-1)*dU/dn); dimensionless momentum equation
%plate velocity
Uo=10;
k=2;
omega=pi/2;

%power law model -> Mu(gamma_dot)= K*(gamma_dot)^(n-1), k-consistency
%inital reference shear rate = Uo*(w/vo)^
%coefficient and n= power law index(non-dimensional)
K=14.67e-3;
n=0.76;
epsilon_o = 0.012 ; % reference viscosity of the newtonian fluid
mu_o= 1.2*10^-3 ;       % zero shear viscosity

%no of meshes;
N=100;
%discretize our space 
y_vec= linspace(0,10,N);
dy=y_vec(2)-y_vec(1); %dx
%discretize our time
dt=0.5*(dy^2)/(2*k);
t_vec= 0:dt:10;

U_mat=zeros(length(y_vec),length(t_vec));
%%Boundary Condition
U_mat(:,1)=0;
U_mat(end,:)=0;
for i=1:length(t_vec)
    U_mat(1,i)=cos(t_vec(i)); % sin(tau) or cos(tau)
end

r=dt/(dy^2);
a=[];
f=[];
for j=1:length(t_vec)-1
    A=zeros(length(y_vec)-2,length(y_vec)-2);  %%%%depends on y_vec
    C=zeros(length(y_vec)-2,1);
    
    for n=2:length(y_vec)-1
        if n==2
          a1= K*(Uo*(omega/epsilon_o)^0.5)^(n-1);%K*((U_mat(n+1,j)-U_mat(n-1,j))/2)/(2*dy*mu_o);
          a2= K*(Uo*(omega/epsilon_o)^0.5)^(n-1);
        elseif n~=2
          a1= K*((U_mat(n+1,j)-U_mat(n-1,j))/2)/(2*dy*mu_o);
          a2= K*((U_mat(n,j)-U_mat(n-2,j))/2)/(2*dy*mu_o);
%         else
%           a2= K*(Uo*(omega/epsilon_o)^0.5)^(n-1);%K*((U_mat(n,j)-U_mat(n-2,j))/2)/(2*dy*mu_o);
        end
        %%%%%%%%
        if a1==0 && a2==0
            f1=0;
        else
            f1=(2*a1*a2/(a1+a2));
        end
        %%%%%%%%
          
        if n==2
          a3= K*(Uo*(omega/epsilon_o)^0.5)^(n-1);%K*((U_mat(n+2,j)-U_mat(n,j))/2)/(2*dy*mu_o);
        elseif n>2 && n<=length(y_vec)-2
          a3= K*((U_mat(n+2,j)-U_mat(n,j))/2)/(2*dy*mu_o);
        end
        %%%%%%%%%%%%%%%        
        if a1==0 && a3==0
            f2=0;
        else
            f2=(2*a1*a3/(a1+a3));
        end
        
        %FDM discretized equation for governing momentum equation;
        %power law fluid;
        %p*U_mat(n-1,j+1)+q*U_mat(n,j+1)+r*U_mat(n+1,j+1)=s;
       
        if n==2
          A(n-1,1)=   (1+ (r/2)*(f1+f2));
          A(n-1,2)=   -r/2*f2;
          C(n-1,1)=    U_mat(n,j)*(1-r)+r/2*(U_mat(n-1,j)+U_mat(n+1,j))+(r/2*f1*U_mat(n-1,j+1));
        elseif n==length(y_vec)-1
          A(n-1,n-2)= (r/2*f1*U_mat(n-1,j+1));
          A(n-1,n-1)= (1+ (r/2)*(f1+f2));
          C(n-1,1)=    U_mat(n,j)*(1-r)+r/2*(U_mat(n-1,j)+U_mat(n+1,j))+(r/2*f2*U_mat(n+1,j+1));
        elseif n<=length(y_vec)-2
          A(n-1,n-2)= (r/2*f1*U_mat(n-1,j+1));
          A(n-1,n-1)= (1+ (r/2)*(f1+f2));
          A(n-1,n)=   -r/2*f2;
          C(n-1,1)=   U_mat(n,j)*(1-r)+r/2*(U_mat(n-1,j)+U_mat(n+1,j));
        end
          
    end
    U=TDMA(A,C);
    for i=1:length(U)
        U_mat(i+1,j+1)=U(i,1);
    end
    
end

for j=1:length(t_vec)-1
    A1=zeros(length(y_vec)-2,length(y_vec)-2);  %%%%depends on y_vec
    C1=zeros(length(y_vec)-2,1);
    
    for n=2:length(y_vec)-1
        if n==2
          a1= K*(Uo*(omega/epsilon_o)^0.5)^(n-1);%K*((U_mat(n+1,j)-U_mat(n-1,j))/2)/(2*dy*mu_o);
          a2= K*(Uo*(omega/epsilon_o)^0.5)^(n-1);
        elseif n~=2
          a1= K*((U_mat(n+1,j)-U_mat(n-1,j))/2)/(2*dy*mu_o);
          a2= K*((U_mat(n,j)-U_mat(n-2,j))/2)/(2*dy*mu_o);
%         else
%           a2= K*(Uo*(omega/epsilon_o)^0.5)^(n-1);%K*((U_mat(n,j)-U_mat(n-2,j))/2)/(2*dy*mu_o);
        end
        %%%%%%%%
        if a1==0 && a2==0
            f1=0;
        else
            f1=(2*a1*a2/(a1+a2));
        end
        %%%%%%%%%%%%%%
          
        if n==2
          a3= K*(Uo*(omega/epsilon_o)^0.5)^(n-1);%K*((U_mat(n+2,j)-U_mat(n,j))/2)/(2*dy*mu_o);
        elseif n>2 && n<=length(y_vec)-2
          a3= K*((U_mat(n+2,j)-U_mat(n,j))/2)/(2*dy*mu_o);
        end
        %%%%%%%%%%%%%%%        
        if a1==0 && a3==0
            f2=0;
        else
            f2=(2*a1*a3/(a1+a3));
        end
        
        %FDM discretized equation for governing momentum equation;
        %power law fluid;
        %p*U_mat(n-1,j+1)+q*U_mat(n,j+1)+r*U_mat(n+1,j+1)=s;
       
        if n==2
          A1(n-1,1)=   (1+ (r/2)*(f1+f2));
          A1(n-1,2)=   -r/2*f2;
          C1(n-1,1)=    U_mat(n,j)*(1-r)+r/2*(U_mat(n-1,j)+U_mat(n+1,j))+(r/2*f1*U_mat(n-1,j+1));
        elseif n==length(y_vec)-1
          A1(n-1,n-2)= (r/2*f1*U_mat(n-1,j+1));
          A1(n-1,n-1)= (1+ (r/2)*(f1+f2));
          C1(n-1,1)=    U_mat(n,j)*(1-r)+r/2*(U_mat(n-1,j)+U_mat(n+1,j))+(r/2*f2*U_mat(n+1,j+1));
        elseif n<=length(y_vec)-2
          A1(n-1,n-2)= (r/2*f1*U_mat(n-1,j+1));
          A1(n-1,n-1)= (1+ (r/2)*(f1+f2));
          A1(n-1,n)=   -r/2*f2;
          C1(n-1,1)=   U_mat(n,j)*(1-r)+r/2*(U_mat(n-1,j)+U_mat(n+1,j));
        end
          
    end
    U=TDMA(A1,C1);
    for i=1:length(U)
        U_mat(i+1,j+1)=U(i,1);
    end
    
end

for i=1:length(t_vec)
 v(i,:)=exp(-y_vec/2^0.5).*sin(t_vec(i)-y_vec/2^0.5);
 DU(i,:)= (-1/(2^0.5))*exp(-y_vec/(2^0.5)).*(-sin(t_vec(i)-y_vec/(2^0.5))+cos(t_vec(i)-y_vec/(2^0.5)));
end
DU=DU';
U_mat=v';
theta_mat=zeros(length(y_vec),length(t_vec));

%%%% Boundary Conditions
theta_mat(:,1)=0;
theta_mat(end,:)=0;
theta_mat(1,:)=1;

pr=1.5;
Ek=10;
for j=1:length(t_vec)-1
    A2=zeros(length(y_vec)-2,length(y_vec)-2);  
    C2=zeros(length(y_vec)-2,1);
    
    for n=2:length(y_vec)-1
        
        %FDM discretized equation for governing energy equation;
        %power law fluid;
        %p*theta_mat(n-1,j+1)+q*theta_mat(n,j+1)+r*theta_mat(n+1,j+1)=s;
       
        if n==2
          A2(n-1,1)=   (1+ r/pr);
          A2(n-1,2)=   -(r/2)*(1/pr);
          C2(n-1,1)=    theta_mat(n,j)*(1-r/pr)+ r/2*(1/pr)*(theta_mat(n-1,j)+theta_mat(n+1,j))+(r/2*(1/pr)*theta_mat(n-1,j+1))+ Ek*(K/mu_o)*(DU(n,j)^(n-1))*DU(n,j)^2;
        elseif n==length(y_vec)-1
          A2(n-1,n-2)= -(r/2*(1/pr)*theta_mat(n-1,j+1));
          A2(n-1,n-1)= (1+ r/pr);
          C2(n-1,1)=    theta_mat(n,j)*(1-r/pr)+ r/2*(1/pr)*(theta_mat(n-1,j)+theta_mat(n+1,j))+(r/2*(1/pr)*theta_mat(n+1,j+1))+ Ek*(K/mu_o)*(DU(n,j)^(n-1))*DU(n,j)^2;
        elseif n<=length(y_vec)-2
          A2(n-1,n-2)= -(r/2*(1/pr)*theta_mat(n-1,j+1));
          A2(n-1,n-1)= (1+ r/pr);
          A2(n-1,n)=   -(r/2)*(1/pr);
          C2(n-1,1)=   theta_mat(n,j)*(1-r/pr)+ r/2*(1/pr)*(theta_mat(n-1,j)+theta_mat(n+1,j))+ Ek*(K/mu_o)*(DU(n,j)^(n-1))*DU(n,j)^2;
        end
          
    end
    
    theta=TDMA(A2,C2); %%using TDMA to solve the linear equation;
    for i=1:length(U)
        theta_mat(i+1,j+1)=theta(i,1);
    end
    
end


%plot(U_mat(:,50),y_vec,'-r');
%hold on;
 
[tt,yy]=meshgrid(t_vec,y_vec);

mesh(yy,tt,U_mat)

xlabel('eta ');
ylabel('Time(t)');
zlabel('theta(U)');







