clear 
%close all
load iterate_data\1000r.mat
fig=figure;
axis=gca;
hold on

names={'i',...%1
    'R1',...%2
    'R2',...%3
    'IMU drift',...%4
    'Visodom drift',...%5
    'Vel Kdrift',...%6
    'Clo Kdrift',...%7
    'Vel Kuncert',...%8
    'Clo Kuncert'};%9

DATA=[(1:iter).' R1s R2s imu_drifts vis_drifts vel_kdrifts, clo_kdrifts, vel_kuncerts clo_kuncerts];

x=4;y=5;filt1=6;filt2=7;

%filt1=7;filt2=9;x=7;y=9;
filt=DATA(:,filt1)>=DATA(:,filt2);

%x=imu_drifts;y=vis_drifts;
s=scatter3(axis,DATA(:,x),DATA(:,y),DATA(:,filt1)-DATA(:,filt2),20,filt,'filled');
set(axis,'XScale','Log')
set(axis,'YScale','Log')
%set(axis,'ZScale','Log')
xlabel(names{x})
ylabel(names{y})

for i=1:length(names)
    s.DataTipTemplate.DataTipRows(i) = dataTipTextRow(names{i},DATA(:,i));
end
colormap([0 0 1;0 0 1;0 0 0;1 0 0;1 0 0])
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
%title('Kalman uncertainty')
%title('Kalman drift')
%legend imu vis

drawnow()
