%% State declaration
syms x          %x-coordinate
syms y          %y-coordinate
syms theta         %yaw
syms x_s         %
syms y_s         %stationary copy
syms theta_s        %
X=[x;y;theta;x_s;y_s;theta_s];
Xs=[0;0;0;0;0;0];
N=length(X);

%% Input declaration
syms V_n          %Normal Speed
syms V_e          %Lateral Speed
syms delta           %Steering angle

U=[V_n;V_e;delta];
Us=[1;5000;5*3.14/180];

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

%% ODE definition
F= [V_n*cos(theta)-V_e*sin(theta);
    V_n*sin(theta)+V_e*cos(theta); 
    V_n*delta/(a+b); 
    0;
    0;
    0
  ]; %a=lf, b=lr

h=[ cos(theta_s)*(x - x_s) - sin(theta_s)*(y - y_s);
    sin(theta_s)*(x - x_s) + cos(theta_s)*(y - y_s);
	theta - theta_s];


generate_c_functions;

