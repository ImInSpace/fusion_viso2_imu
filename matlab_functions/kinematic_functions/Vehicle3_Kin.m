%% State declaration
syms pos_x pos_y pos_z %pos
pos=[pos_x pos_y pos_z].';
syms psi theta phi       %euler yrp
eul=[psi theta phi].';
syms vl_x vl_y vl_z      %local v
vl= [vl_x vl_y vl_z].'; 
syms wl_x wl_y wl_z      %local w
wl= [wl_x wl_y wl_z].';

X=[pos;eul;vl;wl];
Xs=zeros(size(X));
N=length(X);

%% Input declaration

U=sym([]);
Us=[];

%% Parameter declaration
p=[].';
P=[].';

%% ODE definition
F= sym(zeros(size(X))); %a=lf, b=lr
R=my_eul2rotm(eul);
lw2euld=localw2euld(eul);

small_angles_approx=0
if small_angles_approx
    R=subs(R,[sin([psi phi]) cos([psi phi])],[psi phi 1 1]);
    lw2euld=subs(lw2euld,[sin([psi phi]) cos([psi phi])],[psi phi 1 1]);
end

F(1:3)=R*vl;%pos_dot
F(4:6)=lw2euld*wl;

h=[ pos; eul];


generate_c_functions;

function rotm=my_eul2rotm(eul)
    psi=eul(1);theta=eul(2);phi=eul(3);
    cps=cos(psi);sps=sin(psi);
    cth=cos(theta);sth=sin(theta);
    cph=cos(phi);sph=sin(phi);
    rotm=[cps*cth cps*sph*sth-cph*sps sph*sps+cph*cps*sth;
        cth*sps cph*cps+sph*sps*sth cph*sps*sth-cps*sph;
        -sth cth*sph cph*cth];
end

function lw2euld=localw2euld(eul)
    psi=eul(1);theta=eul(2);phi=eul(3);
    cps=cos(psi);sps=sin(psi);
    cth=cos(theta);sth=sin(theta);
    cph=cos(phi);sph=sin(phi);
    lw2euld=[0 sph/cth cph/cth;
    0 cph -sph;
    1 sph*sth/cth cph*sth/cth];
end
