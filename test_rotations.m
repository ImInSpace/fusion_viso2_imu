close all

car1=stlread('car_1.stl');
V=car1.Points;
F=car1.ConnectivityList;
V(:,3)=-V(:,3);
V(:,1)=-V(:,1);
figure
pat=patch('Vertices',V,'Faces',F,'FaceColor','r');
V=V.';
view(3)
axis equal
xlabel x
ylabel y
zlabel z
xlim([-3 3])
ylim([-3 3])
zlim([0 2])
camlight
eul=[0 0 15].'*pi/180; %ypr
p=[0 0 0].';
w=[0 0 pi].';
v=[0 0 0].';
dt=0.02;

try
    rosinit;
end
r=rosrate(1/dt);
for it=1:ceil(2/dt)
    euld=localw2euld(eul,w);
    eul=eul+euld*dt;
    R=eul2rotm(eul.');
    p=p+R*v*dt;
    pat.set('Vertices',(p+R*V).')
    waitfor(r);
end
disp(p)
close all

syms psi theta phi
ypr=[psi;theta;phi];
eulds=localw2euld(ypr,w);
latex_show(eulds)

function euld=globalw2euld(eul,w)
    psi=eul(1);theta=eul(2);phi=eul(3);
    cps=cos(psi);sps=sin(psi);
    cth=cos(theta);sth=sin(theta);
    cph=cos(phi);sph=sin(phi);
    euld=[cps*sth/cth sps*sth/cth 1;
          -sps cps 0;
          cps/cth sps/cth 0]*w;
end
function euld=localw2euld(eul,w)
    psi=eul(1);theta=eul(2);phi=eul(3);
    cps=cos(psi);sps=sin(psi);
    cth=cos(theta);sth=sin(theta);
    cph=cos(phi);sph=sin(phi);
    euld=[0 sph/cth cph/cth;
        0 cph -sph;
        1 sph*sth/cth cph*sth/cth]*w;
end
