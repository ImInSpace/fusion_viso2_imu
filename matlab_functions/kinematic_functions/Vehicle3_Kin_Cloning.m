%% State declaration
syms x          %x-coordinate
syms y          %y-coordinate
syms z         %yaw
syms q1         %
syms q2         %stationary copy
syms q3        %
syms q4        %
X=[x;y;z;q1;q2;q3;q4];
Xs=[0;0;0;0;0;0;0];
N=length(X);

%% Input declaration
syms vx          %Normal Speed
syms vy          %Lateral Speed
syms vz          %Steering angle

U=[V_n;V_e;theta_dot];
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
F= [V_n*cos(theta)-V_e*sin(theta);
    V_n*sin(theta)+V_e*cos(theta); 
    theta_dot; 
    0;
    0;
    0
  ]; %a=lf, b=lr

h=[ cos(theta_s)*(x - x_s) + sin(theta_s)*(y - y_s);
   -sin(theta_s)*(x - x_s) + cos(theta_s)*(y - y_s);
	theta - theta_s];


generate_c_functions;

