system('cmake-build-debug\test_fusion_viso2_imu.exe 5 0 config.yaml');
X=csvread("exoter_logs\log_1_filtered.csv").';
config=yaml.ReadYaml('config.yaml');
configCase=config.testCase{5};
dt=configCase.dt;
N=size(X,2);
ps=zeros(3,N);
euls=zeros(3,N);
ps(:,1)=X(1:3,1);
euls(:,1)=X(4:6,1);
for i=2:N
    vl=X(7:9,i);wl=X(10:12,i);
    lw2euld=localw2euld(euls(:,i-1));
    euls(:,i)=euls(:,i-1)+lw2euld*wl*dt;
    R=my_eul2rotm(euls(:,i));
    ps(:,i)=ps(:,i-1)+R*vl*dt;
end
figure
hold on
plot3(X(1,:),X(2,:),X(3,:),'r')
plot3(ps(1,:),ps(2,:),ps(3,:),'b')
axis equal
xlabel x
ylabel y
zlabel z
title(sprintf("dt:%g, Qm:%g, Rm:%g",configCase.dt,configCase.Qm,configCase.obs{1}.Rm))
