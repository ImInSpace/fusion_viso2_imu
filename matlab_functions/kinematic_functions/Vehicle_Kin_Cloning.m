%% State declaration
syms x          %x-coordinate
syms y          %y-coordinate
syms psi         %yaw
syms x_s         %
syms y_s         %stationary copy
syms psi_s        %
X=[x;y;psi;x_s;y_s;psi_s];
Xs=[0;0;0;0;0;0];
N=length(X);

%% Input declaration
syms V_n          %Normal Speed
syms V_e          %Lateral Speed
syms psi_dot    %Steering angle

U=[V_n;V_e;psi_dot];
Us=[1;5000;0];

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
F= [V_n*cos(psi)-V_e*sin(psi);
    V_n*sin(psi)+V_e*cos(psi); 
    psi_dot; 
    0;
    0;
    0
  ]; %a=lf, b=lr

h=[ cos(psi_s)*(x - x_s) + sin(psi_s)*(y - y_s);
   -sin(psi_s)*(x - x_s) + cos(psi_s)*(y - y_s);
	psi - psi_s];


generate_c_functions;

