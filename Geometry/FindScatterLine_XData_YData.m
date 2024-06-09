% Find all scatter plot objects in the current figure
% scatterObjs = findobj(gca, 'Type', 'Scatter');
scatterObjs = findobj(gca, 'Type', 'Line');

% Initialize arrays to hold all x and y data
allXData = [];
allYData = [];

% Loop through each scatter object to extract data
for i = 1:length(scatterObjs)
    allXData = [allXData; scatterObjs(i).XData(:)];
%     allYData = [allYData; scatterObjs(i).YData(:)];
end


% Display the extracted data
figure();
% plot(allXData(allYData~='nan'), allYData(allYData~='nan'), 'o');
plot(allXData, log(allYData), 'o');

%%


% scatterObjs = findobj(gca, 'Type', 'Line');
scatterObjs = findobj('Type', 'Line');

% Initialize arrays to hold all x and y data
x = [];

% Loop through each scatter object to extract data
for i = 1:length(scatterObjs)
%     allXData = [allXData; scatterObjs(i).XData(:)];
    x = [x, scatterObjs(i).YData(:)];
end

fs = 51.2e3; dt = 1/fs;
len = length(x);
ll = fix(0.220*fs):fix(0.520*fs);
% Display the extracted data
figure(12); clf; hold on;
for idx = 6 : -1:2
    plot(1e3*[0:length(ll)-1]*dt,x(ll,idx));
end
xlabel('Time [ms]');
ylabel('Pressure [Pa]');
grid on; box on;
yyaxis right;
plot(1e3*[0:length(ll)-1]*dt, x(ll,1));
ylabel('Key velocity [m/s]');

