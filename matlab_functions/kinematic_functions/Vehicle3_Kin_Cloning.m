%% State declaration
syms x y z %pos
pos=[x y z].';
syms psi theta phi       %euler yrp
yaw=psi;roll=theta;pitch=phi;
eul=[psi theta phi].';
syms   x_s y_s z_s %static pos
pos_s=[x_s y_s z_s].';
syms   psi_s theta_s phi_s       %static euler
yaw_s=psi_s;roll_s=theta_s;pitch_s=phi_s;
eul_s=[psi_s theta_s phi_s].';

X=[pos;eul;pos_s;eul_s];
Xs=zeros(size(X));
N=length(X);

%% Input declaration

syms vl_x vl_y vl_z      %local v
vl= [vl_x vl_y vl_z].'; 
syms wl_x wl_y wl_z      %local w
wl= [wl_x wl_y wl_z].';

U=[vl;wl];
Us=zeros(size(U));

%% Parameter declaration
p=[].';
P=[].';

%% ODE definition
R=my_eul2rotm(eul);
lw2euld=localw2euld(eul);

small_angles_approx=0
if small_angles_approx
    R=subs(R,[sin([pitch roll]) cos([pitch roll])],[pitch roll 1 1]);
    lw2euld=subs(lw2euld,[sin([pitch roll]) cos([pitch roll])],[pitch roll 1 1]);
end

euld2lw=simplify(inv(lw2euld),100);

R_s=subs(R,eul,eul_s);
lw2euld_s=subs(lw2euld,eul,eul_s);
euld2lw_s=subs(euld2lw,eul,eul_s);

F= sym(zeros(size(X)));
F(1:3)=R*vl;%dot
F(4:6)=lw2euld*wl;

h=[R_s.'*(pos-pos_s);
   euld2lw_s*(eul-eul_s)];%delta_pos
                          %TODO: delta_eul

                          
                          
%F=subs(F,[pitch roll pitch_s roll_s],[0 0 0 0]);
%h=subs(h,[pitch roll pitch_s roll_s],[0 0 0 0]);
generate_c_functions;




