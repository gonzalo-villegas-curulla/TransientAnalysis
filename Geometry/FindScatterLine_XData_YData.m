% Find all scatter plot objects in the current figure
scatterObjs = findobj(gca, 'Type', 'Scatter');
% scatterObjs = findobj(gca, 'Type', 'Line');

% Initialize arrays to hold all x and y data
allXData = [];
allYData = [];

% Loop through each scatter object to extract data
for i = 1:length(scatterObjs)
    allXData = [allXData; scatterObjs(i).XData(:)];
    allYData = [allYData; scatterObjs(i).YData(:)];
end

figure();
% plot(allXData(allYData~='nan'), allYData(allYData~='nan'), 'o');
% plot(allXData, log(allYData), 'o');
plot(allXData, allYData, 'o');

%%
fs = 51.2e3; dt = 1/fs;
% scatterObjs = findobj(gca, 'Type', 'Line');
objs = findobj('Type', 'Line');

% Initialize arrays to hold all x and y data
x = [];

x = [objs(8).YData',objs(7).YData'];

init = fix(0.0154297*fs);
ll = init:init+fix(0.030*fs);

figure; hold on;
for idx = 1:2
    plot(1e3*[0:length(ll)-1]/fs,x(ll,idx));
end
xlim([0 30]);
ylim([0 250]);
PTARG = 219.206;
plot([0,30],[1,1]*PTARG,'--k');

IDX20 = find(x(ll,1)>0.2*PTARG,1,'first') ;
IDX80 = find(x(ll,1)>0.8*PTARG,1,'first') ;
plot(1e3*dt*IDX20,0.2*PTARG, 'ok');
plot(1e3*dt*IDX80,0.8*PTARG, 'ok');

plot(1e3*dt*IDX20*[1,1],[0,0.2*PTARG],'--k');
plot(1e3*dt*IDX80*[1,1],[0,0.8*PTARG],'--k');


box on;
xlabel('Time [ms]');
ylabel('Pressure [Pa]');


%%
xlabel('Time [ms]');
ylabel('Pressure [Pa]');
grid on; box on;
yyaxis right;
plot(1e3*[0:length(ll)-1]*dt, x(ll,1));
ylabel('Key velocity [m/s]');

