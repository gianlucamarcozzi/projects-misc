%%
clearvars

%% IMPORT
generalFolder = "../data/raw/";
expName = "misc-E-002-004";
loadDir = dir(generalFolder + expName + "/*.DTA");
nMeas = numel(loadDir);

for ii = 1:nMeas
    filename = append(loadDir(ii).folder, '/', loadDir(ii).name);
    [x{ii}, y{ii}, param{ii}] = eprload(filename);
end

xAmp = table2array(readtable(append(...
    generalFolder, expName, '/', expName, '-param.txt'), "Range", "B1"));
% Check for nans
isnanxAmp = isnan(xAmp);
if sum(isnanxAmp)
    for ii = 1:numel(xAmp)
        if isnanxAmp(ii) == 1
            xAmp = xAmp(1:ii - 1);
            break
        end
    end
end

%---------------------------- FIT ----------------------------------------
expcos = @(xx, p) [ones(numel(xx), 1), exp(-xx/p(1)).*cos(2*pi*p(2)*xx)];
fitOpt = optimoptions('lsqnonlin','Display','off');

figure(1)
clf
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
w0 = linspace(1, sqrt(23), nMeas).^2*1e-3;
for ii = 1:nMeas
    xdata = x{ii};
    ydata = real(y{ii});
    ydata = ydata/max(ydata);

    fitmodel = @(p) expcos(xdata, p);
    p0 = [100, w0(ii)];
    
    [yfit{ii}, pfit{ii}, pci{ii}] = lsqnonlin2steps(...
        ydata, fitmodel, p0, fitOpt);

    freq(ii) = pfit{ii}(2)*1e3;  % MHz
    piPerf(ii) = 1/2/freq(ii)*1e3;  % ns
    dfreq(ii) = pci{ii}(2)*1e3;

    nexttile
    plot(xdata, ydata, 'o')
    hold on
    plot(xdata, yfit{ii})
    title(sprintf(...
        '%d: %.4f MHz, %.2f ns', ii, freq(ii), piPerf(ii)))
end

%-------------------------- PLOT ALL TOGETHER ----------------------------
figure(2)
clf
errorbar(xAmp, freq, dfreq, 'o')
labelaxesfig(gca, 'mpfu amplitude / %', 'Rabi freq / MHz')
xlim(setaxlim(xAmp, 1.05))

%% FIT SIGMOIDAL AND EXTRACT EQUIDISTANT PULSES
xx = linspace(min(xAmp), max(xAmp), 1000);
sigmoidalfun = @(xx, p) [...
    1./(1 + p(1)*exp(-p(2)*(xx - p(3)))), ones(numel(xx), 1) ...
    ];
fitmodel = @(p) sigmoidalfun(xAmp, p);

ydata = freq';
p0 = [5, 0.6, 48];
[yfit2, pfit2, pci2] = lsqnonlin2steps(ydata, fitmodel, p0, fitOpt);

sigmoidalfun0 = @(xx, p) p(4)./(1 + p(1)*exp(-p(2)*(xx - p(3)))) + p(5);
ympfu = sigmoidalfun0(xx, pfit2);

clf
plot(xAmp, freq, 'o', xAmp, yfit2, xx, ympfu)
hold on

ny = 25;
ygoal = linspace(ympfu(1), ympfu(end), ny);
ixgoal = zeros(ny, 1);
for ii = 1:ny
    [~, ixgoal(ii)] = min(abs(ympfu - ygoal(ii)));
end
xgoal = xx(ixgoal);

figure(3)
clf
plot(xgoal, ympfu(ixgoal), 'k', 'Marker', 'square', 'LineStyle', 'none')

%% SAVE TO FILE
fileID = fopen(append(expName, '_xAmpEquidistant.txt'), 'w');
for ii = 1:ny - 1
    fprintf(fileID, '%.3f, ', xgoal(ii));
end
fprintf(fileID, '%.3f', xgoal(end));
fclose(fileID);

