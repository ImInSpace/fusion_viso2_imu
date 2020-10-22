close all

[car,V]=plot_car();
hold on
view(3)
axis equal
xlabel x
ylabel y
zlabel z
xlim([-3 3])
ylim([-3 3])
zlim([0 2])
eul=[0 0 0].'*pi/180; %ypr
p=[0 0 0].';
w=[0 0 1].';
v=[1 0 0].';
t=0;
dt=0.02;
ts=t;ps=p;euls=eul;
pl=plot3(ps(1,:),ps(2,:),ps(3,:));

try
    rosinit;
end
r=rosrate(1/dt);
for it=1:ceil(2/dt)
    t=t+dt;
    lw2euld=localw2euld(eul);
    eul=eul+lw2euld*w*dt;
    R=my_eul2rotm(eul);
    p=p+R*v*dt;
    ts(end+1)=t;ps(:,end+1)=p;euls(:,end+1)=eul;
    car.set('Vertices',(p+R*V.').')
    pl.set('XData',ps(1,:),'YData',ps(2,:),'ZData',ps(3,:));
    waitfor(r);
end

csvwrite('fake_logs\log_1.csv',[ps.' euls.']);
csvwrite('fake_logs\log_1_time.csv',ts.');
disp([p eul])
close all

syms psi theta phi
ypr=[psi;theta;phi];
eulds=localw2euld(ypr)*w;
%latex_show(eulds)

function rotm=my_eul2rotm(eul)
    psi=eul(1);theta=eul(2);phi=eul(3);
    cps=cos(psi);sps=sin(psi);
    cth=cos(theta);sth=sin(theta);
    cph=cos(phi);sph=sin(phi);
    rotm=[cps*cth cps*sph*sth-cph*sps sph*sps+cph*cps*sth;
        cth*sps cph*cps+sph*sps*sth cph*sps*sth-cps*sph;
        -sth cth*sph cph*cth];
end
function euld=localw2euld(eul)
    psi=eul(1);theta=eul(2);phi=eul(3);
    cps=cos(psi);sps=sin(psi);
    cth=cos(theta);sth=sin(theta);
    cph=cos(phi);sph=sin(phi);
    euld=[0 sph/cth cph/cth;
        0 cph -sph;
        1 sph*sth/cth cph*sth/cth];
end
