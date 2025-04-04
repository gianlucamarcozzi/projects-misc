clearvars, clc %, close all
% Pulse lengths vary from 8ns to 24ns
figure(1)
clf
tL = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

filename = "../data/raw/misc-E-003001";
yy = readtable(filename + ".Wfm.csv");
yrm = yy{:, 1};
% ytm = yy{:, 2};
Param = readtable(filename + ".csv");
xx = linspace(Param{15, 2}, Param{16, 2}, Param{17, 2})*1e9;

figure(1)
clf
tL = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
plot(xx, yrm)
xlim(setaxlim(xx, 1))
% plot(xx, ytm)

idx(1) = 1;
for ii = 0:7
    idx(ii+ 2) = 175000 + ii*100000; 
    xline(xx(idx(ii + 2)))
end
idx(end + 1) = size(yrm, 1);
% legend('RM');
% savePath = '../images/RM_TM_pulseLengthDependence.png';
%
% for ii = 1:9
%     idx1 = 501:2750;
%     idx2 = (idx1(end) + 1):numel(xx);
%     if numel(idx1) ~= numel(idx2)
%         error("idx1 and idx2 should be same length")
%     end
%     nx = numel(idx1);
% end

nexttile
% clf
% Data correction
% Subtract constant background
yrm = yrm - mean(yrm(1:100000));  

for ii = 1:9
    rangeii = (idx(ii) + 1):idx(ii + 1);
    % Separate first pulse from second one
    p0{ii} = yrm(rangeii);
    xx0{ii} = xx(rangeii);
    plot(xx0{ii}, p0{ii})
    hold on
end

labelaxesfig(tL, 'Time / ns', '');

xlim(setaxlim(xx, 1))
nx = numel(xx0{1});

%% FFT

tStep = xx(2) - xx(1);  % ns
fSampl = 1/tStep;
nzf = 2^20;  % Zero filling
if nzf <= nx && nzf ~= 0
    warning("nzf <= nx. Continuing without zero-filling.")
    nzf = 0;
end
if nzf == 0
    ff = fSampl/(nx)*(-(nx)/2:(nx)/2 - 1);
else
    ff = fSampl/(nzf)*(-(nzf)/2:(nzf)/2 - 1);
end

ff = ff*1e3;  % MHz
figure(2)
viri = viridis(9);
clf
tiledlayout(3, 1, "TileSpacing", "compact", "Padding", "compact")
ax1 = nexttile;
for ii = 9:-1:1

    if nzf ~= 0  % Zero filling
        p0{ii}(nzf) = 0;  % Zero filling
        % p1{ii}(nzf) = 0;  % Zero filling
    else  % No zero filling (adjust values in memory in case)
        p0{ii} = p0{ii}(1:nx);
        % p1{ii} = p1{ii}(1:nx);
    end

    fp0{ii} = fft(p0{ii});
    
    plot(ff, abs(fftshift(fp0{ii})), 'o-', 'DisplayName', num2str(6 + ii*2), 'Color', viri(ii, :))
    hold on
end
% ax1 = gca;
xlimRatio = 0.0005;
xlim(ax1, setaxlim(ff, xlimRatio))
legend()
box(ax1, 'off')
nexttile

sload = load("../../zech-psi/data/processed/ZePSI-E-007015.mat");
bfield = sload.x{2} - (mean(sload.x{2}) + 2);
spc = sum(sload.y(120:260, :));

ffb = mhz2mt(ff)*10;

for ii = 5:-2:3
    transFunc = interp1(ffb, abs(fftshift(fp0{ii})), bfield, "linear");
    plot(bfield, transFunc.*spc', 'DisplayName', num2str(6 + ii*2), 'Color', viri(ii, :))
    hold on
    trapz(abs(transFunc.*spc'))
end

lim1 = get(gca, 'YLim');
ratioLim1 = lim1(2)/lim1(1);

yyaxis right
plot(bfield, spc, 'r');
xlim(setaxlim(ffb, xlimRatio))
% Make the zero match the other axis
lim2 = get(gca, 'YLim');
ratioLim2 = lim2(2)/lim2(1);
set(gca, "YLim", [lim2(1), lim2(2)*ratioLim1/ratioLim2])
legend()

nexttile  % Normalized
for ii = 5:-2:3
    transFunc = interp1(ffb, abs(fftshift(fp0{ii})), bfield, "linear");
    yplot = transFunc.*spc';
    yplot = yplot/max(yplot);
    plot(bfield, yplot, 'DisplayName', num2str(6 + ii*2), 'Color', viri(ii, :))
    hold on
end
plot(bfield, spc/max(spc), 'r');
xlim(setaxlim(ffb, xlimRatio))
legend()

ax2 = axes('Position', ax1.Position, 'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', 'Color', 'none');
yticks(ax2 , [])
xlim(ax2, setaxlim(ffb, xlimRatio))
xline([-10, 10])

%%
figure(4)
clf
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
for ii = 1:9
    nexttile
    plot(1:nzf, p0{ii})
end
% Secondary x-axis
% for ii = 1:9
    % ax2(ii) = axes('Position', ax1(ii).Position, 'XAxisLocation', 'top', ...
        % 'YAxisLocation', 'right', 'Color', 'none');
    % ffb = mhz2mt(ff)*10;
    % xlim(ax2(ii), setaxlim(ffb, 1))
    % yticks(ax2(ii) , [])
    % xlim(ax2(ii), setaxlim(ffb, 0.05))
    % xline([-10, 10])

    % line(ff2, abs(fftshift(fp0{ii})), 'Parent', ax2(ii), 'Color', 'none');
    % 
    % ax2(ii).XColor = 'r';
    % ax2(ii).XLabel.String = 'Secondary X-axis';
    % ax2(ii).XLabel.Color = 'r';
% end


%{
% Create secondary x-axis
ax1 = gca;
ax2 = axes('Position', ax1.Position, 'XAxisLocation', 'top', 'Color', 'none');
x2 = x1 * 2; % Example conversion for secondary x-axis values
line(x2, y1, 'Parent', ax2, 'Color', 'none');

% Link the two x-axes
linkaxes([ax1 ax2], 'x');

% Customize secondary x-axis appearance
ax2.XColor = 'r';
ax2.XLabel.String = 'Secondary X-axis';
ax2.XLabel.Color = 'r';

% Ensure primary axes are active for further modifications
axes(ax1);
%}
