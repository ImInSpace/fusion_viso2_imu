close all
system('cmake-build-debug\test_fusion_viso2_imu.exe 1');
X=csvread("data.csv");
figure()
hold on
scatter(X(1,1),X(1,2),30,'r');
plot(X(:,1),X(:,2),'r')
plot(X(:,3),X(:,4),'g')
plot(X(:,5),X(:,6),'b')
axis equal
legend start ground\_truth observations kalman