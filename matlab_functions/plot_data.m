%close all
%figure()
%for seed=1:9
%addpath(genpath('third-party/yamlmatlab'))
%yamls=yaml.ReadYaml('config.yaml');
%yamls.seed=seed;
%yaml.WriteYaml('config2.yaml',yamls);
%system('cmake-build-debug\test_fusion_viso2_imu.exe 3 0 config2.yaml');
X=csvread("data.csv");
tx=X(:,1); ty=X(:,2);
nObs=X(1,3);
ox=X(:,4:2:(2+2*nObs)); oy=X(:,5:2:(3+2*nObs));
kx=X(:,4+2*nObs); ky=X(:,5+2*nObs);
Pk=reshape(X(:,end-3:end).',[2 2 size(X,1)]);

%subplot(3,3,seed)
hold on
scatter(tx(1),ty(1),30,'r');
legend_text={'start'};
plot(tx,ty,'r');
legend_text{end+1}='ground\_truth';
obs_colors='gcmy';
for i=1:nObs
    li=plot(ox(:,i),oy(:,i),obs_colors(i));
    if (i==2)
        li.LineStyle='-';
        li.Marker='.';
        li.MarkerSize=10;
    end
    legend_text{end+1}=['observation ' num2str(i)];
end
plot(kx,ky,'b')
legend_text{end+1}='kalman';
axis equal
if nObs==1
    obs_legend={'observation'};
end
lgd = legend(legend_text);
lgd.Location='Best';
%title(['Seed ' num2str(seed)])
%drawnow()
%end
return
shadow=polyshape();
for i=1:size(X,1)
    shadow=union(shadow,draw_ellipse([kx(i),ky(i)],Pk(:,:,i)));
end
shplot = plot(shadow,'FaceColor','b',...
            'EdgeColor','b',...
            'LineStyle','--',...
            'FaceAlpha',.2,...
            'EdgeAlpha',.5,...
            'PickableParts','none');

legend_text{end+1}='kalman\_uncertainty';
lgd = legend(legend_text);
lgd.Location='Best';
return
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
