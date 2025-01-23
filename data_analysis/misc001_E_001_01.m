clearvars

%% IMPORT

generalFolder = "../data/raw/";
expName = "misc-E-002001";
measFolder = generalFolder + expName;

[x0, y0, Param0] = loadfolderelexsys(measFolder);
nMeas = numel(Param0);
x = x0{1};

temp = readtable(measFolder + "/misc-E-002001-param.txt");
pulseAmp = temp{:, 2};

%% 

disp('Fit...')
normModel = @(xx, p) p(1)*exp(-(xx - p(2)).^2/2/p(3)^2);
fitmodel = @(p) normModel(x, p);
p0 = [0, 200, 3];

% Smooth
for ii = 1:nMeas
    % y{ii} = datasmooth(y0{ii}, 10, 'savgol');
    y{ii} = y0{ii};
end
fitOpt = optimoptions('lsqnonlin','Display','off');
for ii = 1:nMeas
    % Real
    ydata = real(y{ii});
    p0(1) = max(abs(ydata));
    if -min(ydata) > max(ydata)
        p0(1) = -p0(1);
    end
    [pfit{ii, 1}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
        @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
    pci(ii, 1, :, :) = ...
        nlparci(pfit{ii, 1}, residual, 'jacobian', jacobian);
    [yfit{ii}, pfit{ii, 1}(1), pfit{ii, 1}(4)] = mldividefun(...
        fitmodel, ydata, pfit{ii, 1});

    % Imag
    ydata = imag(y{ii});
    p0(1) = max(abs(ydata));
    if -min(ydata) > max(ydata)
        p0(1) = -p0(1);
    end
    [pfit{ii, 2}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
        @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
    pci(ii, 2, :, :) = ...
        nlparci(pfit{ii, 2}, residual, 'jacobian', jacobian);
    [yfitImag, pfit{ii, 2}(1), pfit{ii, 2}(4)] = mldividefun(...
        fitmodel, ydata, pfit{ii, 2});
    yfit{ii} = yfit{ii} + 1i*yfitImag;
end

for ii = 1:nMeas
    % vre = max(real(yfit{ii}));
    % vim = min(imag(yfit{ii}));
    vre(ii) = trapz(real(yfit{ii}) - pfit{ii, 1}(4));
    vim(ii) = trapz(imag(yfit{ii} - 1i*pfit{ii, 2}(4)));
    
    phase(ii) = atan(vim(ii)/vre(ii));
    radius(ii) = sqrt(vre(ii)^2 + vim(ii)^2);
end
for ii = 1:19
    % phase(ii) = 0;  % Fit on the imaginary part is not reliable
end

yshift = -500;
figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact") 
for ii = 11:30
    nexttile
    plot(x, real(y{ii}), 'o-')
    hold on
    plot(x, real(yfit{ii}))
    plot(x, imag(y{ii}) + yshift, 'o-')
    plot(x, imag(yfit{ii}) + yshift)
end

idxs = [1, 4:13, 15:51];
figure(2)
clf
plot(pulseAmp, phase, 'o-')
hold on
plot(pulseAmp(idxs), phase(idxs), 'o-')
yyaxis right
plot(pulseAmp, radius, 'o-')

fileID = fopen('output.txt', 'w');
for ii = 1:numel(phase)
    fprintf(fileID, '%.2f, ', phase(ii));
end
fclose(fileID);

%% PLOTS

ii = 11;
yshift = -750;
yre = real(y{ii});
yim = imag(y{ii}) + yshift;
yfre = real(yfit{ii});
yfim = imag(yfit{ii}) + yshift;
% figure()
clf
plot(x, yre)
hold on
plot(x, yfre)
plot(x, yim)
plot(x, yfim)
xlim(setaxlim(x, 1))
ylim(setaxlim([yre, yim], 1.05))
labelaxesfig(gca, "Time / ns", "Intensity / a.u.")
legend("0deg ch.", "Fit 0deg", "90deg ch.", "Fit 90deg")
saveImageFolder = "../images/";
imageName = "misc-E-001_exampleEcho.png";
saveImagePath = saveImageFolder + imageName;
% saveas(gcf, saveImagePath)

%%

% figure()
% idxs = [1, 4:13, 15:51];
figure()
clf
plot(pulseAmp, phase, 'ko-')
hold on
ylabel("Phase / rad")
% plot(pulseAmp(idxs), phase(idxs), 'o-')
yyaxis right
plot(pulseAmp, radius, 'o-')
xlim(setaxlim(pulseAmp, 1.05))
% ylim(setaxlim([yre, yim], 1.05))
labelaxesfig(gca, "MPFU amplitude / %", "Intensity / a.u.")
lg = legend("Phase " + char(hex2dec('03C6')), ...
    "Intensity " + char(hex2dec('03C1')));
lg.Location = "east";
saveImageFolder = "../images/";
imageName = "misc-E-001_phaseAndIntensity_vs_mpfuAmp.png";
saveImagePath = saveImageFolder + imageName;
% saveas(gcf, saveImagePath)
