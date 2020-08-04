close all
system('cmake-build-debug\test_fusion_viso2_imu.exe');
X=csvread("data.csv");
figure()
hold on
plot(X(:,1),X(:,2))
plot(X(:,3),X(:,4))
legend observations kalman