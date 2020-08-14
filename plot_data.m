close all
system('cmake-build-debug\test_fusion_viso2_imu.exe 3');
X=csvread("data.csv");
tx=X(:,1); ty=X(:,2);
ox=X(:,3); oy=X(:,4);
kx=X(:,5); ky=X(:,6);
figure()
hold on
scatter(tx(1),ty(1),30,'r');
plot(tx,ty,'r')
plot(ox,oy,'g')
plot(kx,ky,'b')
axis equal
legend start ground\_truth observations kalman

Qo=cov(ox-tx,oy-ty);
Qk=cov(kx-tx,ky-ty);
disp(det(Qk))
disp(det(Qo)/det(Qk));