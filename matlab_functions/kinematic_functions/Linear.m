%% State declaration
syms x          %x-coordinate
syms y          %y-coordinate
syms x_dot     %
syms y_dot     %
X=[x;y;x_dot;y_dot];
Xs=[0;0;1;0];
N=length(X);

%% Input declaration
syms x_ddot          %Normal Speed
syms y_ddot          %Lateral Speed
U=[x_ddot;y_ddot];
Us=[0;0.1];

%% ODE definition
F= [x_dot
    y_dot; 
    x_ddot; 
    y_ddot
  ]; %a=lf, b=lr

h=[ x;
    y];


generate_c_functions;

