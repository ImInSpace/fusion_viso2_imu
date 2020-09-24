%% State declaration
syms x          %x-coordinate
syms y          %y-coordinate
syms th         %yaw
syms th_dot     %yaw rate
syms Ve         %lateral velocity
syms Vn         %longitudinal velocity
X=[x;y;th;Vn;Ve;th_dot];
Xs=[0;0;0;0;0;1];
N=length(X);

%% Input declaration
syms Pf          %Force to front tyres
syms Pr          %Force to rear tyres
syms d           %Steering angle
U=[Pf;Pr;d];
Us=[0;5000;5*3.14/180];

%% Parameter declaration
m=1292.2;       %Vehicle mass
I=2380.7;       %Vehicle inertia
a=1.006;        %CG to front axle
b=1.534;        %CG to rear axle
h=0.3;          %CG to ground
g=9.81;         %gravity     
mu=0.85;        %friction coeff
c=20000;        %cornering stiffness
p=[m I a b h g mu c].';
syms m I a b h g mu c
P=[m I a b h g mu c].';
%% Segel Lateral Force Model
ignore_lateral_force=0;
if ~ignore_lateral_force
    Fzf=(m*g*b-(Pf+Pr)*h)/(a+b); %normal load on front tyre
    Fzr=(m*g*a-(Pf+Pr)*h)/(a+b); %normal load on rear tyre

    af=d-(a*th_dot+ Ve)/Vn;  %front wheel slip angle
    ar=(b*th_dot-Ve)/Vn;    %rear wheel slip angle

    %asf=c*af/(mu*Fzf);
    %asr=c*ar/(mu*Fzr); 
    
    Fef=c*af;
    Fer=c*ar;
    %Fef=mu*Fzf*(asf-(asf^2)/3+asf^3/27)*sqrt(1-Pf^2/(mu^2*Fzf^2)+Pf^2/c^2); %Lateral force on front tyre
    %Fer=mu*Fzr*(asr-(asr^2)/3+asr^3/27)*sqrt(1-Pr^2/(mu^2*Fzr^2)+Pr^2/c^2); %Lateral force on rear tyre
else
    %If we assume that Fef always does maximum force to counter 
    Fef=Pr*d/(1+d);
    Fer=-Pr*d/(1+d);
end

%% ODE definition
F= [Vn*cos(th)-Ve*sin(th);     %x=X
    Vn*sin(th)+Ve*cos(th);     %y=Y
    th_dot;                    %th=psi
    (Pr-Fef*sin(d))/m+Ve*th_dot; %Vn=x_dot
    (Fef*cos(d)+Fer)/m-Vn*th_dot;%Ve=y_dot
    (a*Fef-b*Fer)/I;    %th_dot=psi_dot
  ]; %a=lf, b=lr

h=[Vn;     %x=X
   Ve;     %y=Y
   th_dot];

generate_c_functions(F,h,X,U)