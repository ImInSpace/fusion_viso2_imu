%% State declaration
syms x          %x-coordinate
syms y          %y-coordinate
syms psi         %yaw
syms psi_dot     %yaw rate
syms V_e         %lateral V_elocity
syms V_n         %longitudinal V_elocity
X=[x;y;psi;V_n;V_e;psi_dot];
Xs=[0;0;0;0;0;1];

%% Input declaration
syms P_f          %Force to front tyres
syms P_r          %Force to rear tyres
syms delta           %Steering angle
U=[P_f;P_r;delta];
Us=[0;5000;5*3.14/180];

%% Parameter declaration
m=1292.2;       %V_ehicle mass
I=2380.7;       %V_ehicle inertia
a=1.006;        %CG to front axle
b=1.534;        %CG to rear axle
h=0.3;          %CG to ground
g=9.81;         %gravity     
mu=0.85;        %friction coeff
c=20000;        %cornering stiffness
Ps=[m I a b h g mu c].';
syms m I a b h g mu c
P=[m I a b h g mu c].';
%% Segel Lateral Force Model
ignore_lateral_force=0;
if ~ignore_lateral_force
    Fzf=(m*g*b-(P_f+P_r)*h)/(a+b); %normal load on front tyre
    Fzr=(m*g*a-(P_f+P_r)*h)/(a+b); %normal load on rear tyre

    af=delta-(a*psi_dot+ V_e)/V_n;  %front wheel slip angle
    ar=(b*psi_dot-V_e)/V_n;    %rear wheel slip angle

    %asf=c*af/(mu*Fzf);
    %asr=c*ar/(mu*Fzr); 
    
    Fef=c*af;
    Fer=c*ar;
    %Fef=mu*Fzf*(asf-(asf^2)/3+asf^3/27)*sqrt(1-P_f^2/(mu^2*Fzf^2)+P_f^2/c^2); %Lateral force on front tyre
    %Fer=mu*Fzr*(asr-(asr^2)/3+asr^3/27)*sqrt(1-P_r^2/(mu^2*Fzr^2)+P_r^2/c^2); %Lateral force on rear tyre
else
    %If we assume psiat Fef always does maximum force to counter 
    Fef=P_r*delta/(1+delta);
    Fer=-P_r*delta/(1+delta);
end

%% ODE definition
F= [V_n*cos(psi)-V_e*sin(psi);     %x=X
    V_n*sin(psi)+V_e*cos(psi);     %y=Y
    psi_dot;                    %psi=psi
    (P_r-Fef*sin(delta))/m+V_e*psi_dot; %V_n=x_dot
    (Fef*cos(delta)+Fer)/m-V_n*psi_dot;%V_e=y_dot
    (a*Fef-b*Fer)/I;    %psi_dot=psi_dot
  ]; %a=lf, b=lr

Fef_s=sym('F_ef');
Fer_s=sym('F_er');
F=subs(F,[Fef Fer],[Fef_s Fer_s])
h=[V_n;     %x=X
   V_e;     %y=Y
   psi_dot];

generate_c_functions;