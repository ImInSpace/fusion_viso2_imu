close all
system('cmake-build-debug\test_fusion_viso2_imu.exe');
X=csvread("data.csv");
figure()
hold on
plot(X(:,1),X(:,2),'r')
plot(X(:,3),X(:,4),'g')
plot(X(:,5),X(:,6),'b')
axis equal
legend ground\_truth observations kalman