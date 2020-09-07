close all
system('cmake-build-debug\test_fusion_viso2_imu.exe 3 0');
X=csvread("data.csv");
tx=X(:,1); ty=X(:,2);
ox=X(:,3); oy=X(:,4);
kx=X(:,5); ky=X(:,6);
Pk=reshape(X(:,end-3:end).',[2 2 size(X,1)]);

figure()
hold on
scatter(tx(1),ty(1),30,'r');
plot(tx,ty,'r')
plot(ox,oy,'g')
plot(kx,ky,'b')
axis equal
lgd = legend('start','ground\_truth','observations','kalman');
lgd.Location='northwest';
shadow=polyshape();
for i=1:size(X,1)
    shadow=union(shadow,draw_ellipse([kx(i),ky(i)],Pk(:,:,i)));
end
shplot = plot(shadow,'FaceColor','b',...
            'EdgeColor','b',...
            'LineStyle','--',...
            'FaceAlpha',.2,...
            'EdgeAlpha',.5);
legend start ground\_truth observations kalman kalman\_uncertainty
%{
Qo=cov(ox-tx,oy-ty);
Qk=cov(kx-tx,ky-ty);
disp(det(Qk))
disp(det(Qo)/det(Qk));
%}

dist_traveled=[0;cumsum(sqrt(diff(tx).^2+diff(ty).^2))];
dist_error=vecnorm([ox-tx oy-ty],2,2);
%figure()
%hold on
%plot(dist_traveled,dist_error)
dlm = fitlm(dist_traveled,dist_error,'Intercept',false);
fprintf('%.1f%%\n',dlm.Coefficients.Estimate*100)