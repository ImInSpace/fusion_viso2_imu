close all

car1=stlread('car_1.stl');
V=car1.Points;
F=car1.ConnectivityList;
V(:,3)=-V(:,3);
V(:,1)=-V(:,1);
figure
p=patch('Vertices',V,'Faces',F,'FaceColor','r');
view(3)
axis equal
xlabel x
ylabel y
zlabel z
xlim([-3 3])
ylim([-3 3])
zlim([0 2])
camlight
q=[0 0 .5 .5].'; %using quaternions as JPL (natural order)
               %quaternion=q(4)*+q(1)**i+q(2)**j+q(3)**k
w=[0;0;1];
dt=0.1;
q=q/norm(q);

syms w_x w_y w_z
dts=sym('dt');
norms=sym('norm');
ws=[w_x;w_y;w_z];
ex=expm(.5*Omega(ws)*dts);
ex=simplify(ex);
ex=subs(ex,(- w_x^2 - w_y^2 - w_z^2)^(1/2),norms*i);
ex=simplify(ex);
exb=(eye(4)+.25*Omega(ws)*dts)*(eye(4)-.25*Omega(ws)*dts)^-1;
exb=simplify(exb)
try
    rosinit;
end
r=rosrate(pi/dt);
for it=1:ceil(2*pi/dt)
    %q=(eye(4)+.5*Omega(w)*dt)*q;%first order
    %q=(eye(4)-.5*Omega(w)*dt)^-1*q;%inverse
    %q=(eye(4)+.25*Omega(w)*dt)*(eye(4)-.25*Omega(w)*dt)^-1*q;%both
    %q=(eye(4)+.5*Omega(w)*dt+.25*(.5*Omega(w)*dt)^2)*q;%second order
    %q=expm(.5*Omega(w)*dt)*q;%true
    %q=q/norm(q);
    R=rot_mat(q);
    Vr=V*R.';
    p.set('Vertices',Vr)
    waitfor(r);
end
disp(norm(q))
close all



              
function r=quat_mult(q,p)
    r=[q(4)*p(1) + q(3)*p(2) - q(2)*p(3) + q(1)*p(4);
    	-q(3)*p(1) + q(4)*p(2) + q(1)*p(3) + q(2)*p(4);
    	q(2)*p(1) - q(1)*p(2) + q(4)*p(3) + q(3)*p(4);
    	-q(1)*p(1) - q(2)*p(2) - q(3)*p(3) + q(4)*p(4)];
end

function O=Omega(w)
    O=[0 w(3) -w(2) w(1);
        -w(3) 0 w(1) w(2);
        w(2) -w(1) 0 w(3);
        -w(1) -w(2) -w(3) 0];
end

function C=rot_mat(q)
    C=[q(1)^2-q(2)^2-q(3)^2+q(4)^2 2*(q(1)*q(2)+q(3)*q(4)) 2*(q(1)*q(3)-q(2)*q(4));
        2*(q(1)*q(2)-q(3)*q(4)) -q(1)^2+q(2)^2-q(3)^2+q(4)^2 2*(q(2)*q(3)+q(1)*q(4));
        2*(q(1)*q(3)+q(2)*q(4)) 2*(q(2)*q(3)-q(1)*q(4)) -q(1)^2-q(2)^2+q(3)^2+q(4)^2].';
end