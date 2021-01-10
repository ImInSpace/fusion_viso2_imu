clear 
%close all
load iterate_data\1000r.mat
fig=figure;
axis=gca;
hold on

names={'i',...%1
    'R1',...%2
    'R2',...%3
    'Obs1 drift [m]',...%4
    'Obs2 drift [m]',...%5
    'Velocioty drift [m]',...%6
    'Cloning drift [m]',...%7
    'Kalman uncertainty [m]',...%8
    'Cloning uncertainty [m]'};%9
short_names={'i',...%1
    'R1',...%2
    'R2',...%3
    'obs1',...%4
    'obs2',...%5
    'vel',...%6
    'clo',...%7
    'vel_unc',...%8
    'clo_unc'};%9

DATA=[(1:iter).' R1s R2s imu_drifts vis_drifts vel_kdrifts, clo_kdrifts, vel_kuncerts clo_kuncerts];
DATA=rmoutliers(DATA);
DATA=DATA(~any(isnan(DATA),2),:);
x=4;y=5;filt1=8;filt2=9;
%x=6;y=8;filt1=1;filt2=1;
x=7;y=9;filt1=1;filt2=1;

%filt1=7;filt2=9;x=7;y=9;
filt=DATA(:,filt1)>=DATA(:,filt2);

%x=imu_drifts;y=vis_drifts;
if filt1~=filt2
    s=scatter3(axis,DATA(:,x),DATA(:,y),DATA(:,filt1)-DATA(:,filt2),20,filt,'filled');
elseif filt1==1
    s=scatter(axis,DATA(:,x),DATA(:,y),20,'filled');
else
    s=scatter3(axis,DATA(:,x),DATA(:,y),DATA(:,filt1),20,DATA(:,filt1),'filled');
end
set(axis,'XScale','Log')
set(axis,'YScale','Log')
%set(axis,'ZScale','Log')
xlabel(names{x})
ylabel(names{y})
axis equal
drawnow()
xl=get(axis,'XLim');
yl=get(axis,'YLim');
rng=[min([xl yl]) max([xl yl])];
rng=[ 0.0077   18.5];
xlim(rng);
ylim(rng)
grid on
drawnow()
grid minor
for i=1:length(names)
    s.DataTipTemplate.DataTipRows(i) = dataTipTextRow(names{i},DATA(:,i));
end

if filt1==filt2
    name = ['iterate_data/' short_names{x} '_vs_' short_names{y} '.pdf'];
    exportgraphics(axis,name)
    return
end
colormap([0 0 1;0 0 1;0 0 0;1 0 0;1 0 0])
%{
cb = colorbar();
caxis([0 1])
cb.TickLabels = [];
ax = axes('Position', cb.Position,...
    'Color', 'none',...
    'XTick', [],...
    'YLim', [0 1],...
    'YTick', [0.1667 0.5 0.6667],...
    'YTickLabels',{[names{filt1} '<' names{filt2}], '', [names{filt1} '>=' names{filt2}]},...
    'YAxisLocation','right',...
    'FontSize', cb.FontSize);
ytickangle(90)
%}
%title('Kalman uncertainty')
%title('Kalman drift')
%legend imu vis

drawnow()
name = ['iterate_data/' short_names{filt1} '_vs_' short_names{filt2} '.pdf'];
exportgraphics(axis,name)