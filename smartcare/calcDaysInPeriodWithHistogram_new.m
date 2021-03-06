function [brIntChkptPat, brIntChkptSum] = calcDaysInPeriodWithHistogram(brIntChkptPat, brIntChkptSum, perioddata, study, meastype, period1, period2, twindow, type, cutoffd, plotsubfolder)

% calcDaysInPeriodWithHistogram - calculate the days on IV for pre and during
% study and plot histogram for each

twtext = sprintf('%dm', twindow);
pghght = 6;
pgwdth = 8.5;
plotsacross = 2;
plotsdown = 1;
npat = size(brIntChkptPat, 1);

name = sprintf('%s %s %s %s vs %s', study, meastype, twtext, period1, period2);
fprintf('1) %s\n', name);
[fig, pan] = createFigureAndPanelForPaper(name, pgwdth, pghght);
colname1 = sprintf('%s%s%s', meastype, twtext, period1);
colname2 = sprintf('%s%s%s', meastype, twtext, period2);
brIntChkptPat{:, {colname1}} = 0.0;
brIntChkptPat{:, {colname2}} = 0.0;

for p = 1:npat
    scid = brIntChkptPat.ID(p);
    drugd = brIntChkptPat.DTStart(p);
    fromd = brIntChkptPat.DTStart(p) - calmonths(twindow);
    tod   = brIntChkptPat.DTStart(p) + calmonths(twindow);
    if tod > cutoffd
        tod = cutoffd;
    end

    l1 = days(drugd-fromd);
    l2 = days(tod-drugd);

    fprintf('%2d: ID %d, Drug Therapy Date %11s ', p, scid, datestr(drugd, 1));

    pperioddata = perioddata(perioddata.ID == scid, :);

    if ismember(meastype, {'IVDays'})
        %normalise since this is no longer 6 months each
        fprintf(': %s ', period1);
        brIntChkptPat{p, {colname1}} = calcDaysOnIVForPatient(pperioddata,  fromd,  drugd)/l1;

        fprintf(': %s ', period2);
        brIntChkptPat{p, {colname2}} = calcDaysOnIVForPatient(pperioddata,  drugd,    tod)/l2;
    elseif ismember(meastype, {'EmergCont'})
        fprintf(': %s ', period1);
        brIntChkptPat{p, {colname1}} = calcEmergContForPatient(pperioddata, fromd,  drugd)/l1;

        fprintf(': %s ', period2);
        brIntChkptPat{p, {colname2}} = calcEmergContForPatient(pperioddata, drugd,    tod)/l2;
    else
        fprintf('Unknown meastype %s\n', meastype);
        return
    end

    fprintf('\n');
end

% plot histogram for each period
thisplot = 1;
edges = linspace(0,1,11);
ax = subplot(plotsdown, plotsacross, thisplot, 'Parent', pan);
if ismember(meastype, {'IVDays'})
    h = histogram(ax, brIntChkptPat{:, {colname1}}, 'BinEdges', edges);
else
    h = histogram(ax, brIntChkptPat{:, {colname1}});
end
ax.FontSize = 8;
ax.FontName = 'Arial';
title(ax, sprintf('%s %s %s', meastype, twtext, period1));
xlabel(ax, sprintf('%s %s %s', meastype, twtext, period1), 'FontSize', 8);
ylabel(ax, 'Count', 'FontSize', 8);

thisplot = 2;
ax = subplot(plotsdown, plotsacross, thisplot, 'Parent', pan);
if ismember(meastype, {'IVDays'})
    h = histogram(ax, brIntChkptPat{:, {colname2}}, 'BinEdges', edges);
else
    h = histogram(ax, brIntChkptPat{:, {colname2}});
end

ax.FontSize = 8;
ax.FontName = 'Arial';
title(ax, sprintf('%s %s %s', meastype, twtext, period2));
xlabel(ax, sprintf('%s %s %s', meastype, twtext, period2), 'FontSize', 8);
ylabel(ax, 'Count', 'FontSize', 8);

% save plot and close
if exist('fig', 'var')
    savePlotInDir(fig, name, plotsubfolder);
    close(fig);
end

% calc paired t-test and store results + demographics
[~, pval] = ttest(brIntChkptPat{:, {colname1}}, brIntChkptPat{:, {colname2}});
brIntChkptSum.DataType(type)       = {sprintf('%s%s T-Test', meastype, twtext)};
brIntChkptSum.n(type)              = npat;
brIntChkptSum.Period1Mean(type)    = mean(brIntChkptPat{:, {colname1}});
brIntChkptSum.Period2Mean(type)    = mean(brIntChkptPat{:, {colname2}});
brIntChkptSum.Period1StdErr(type)  =  std(brIntChkptPat{:, {colname1}}) / (npat ^ 0.5);
brIntChkptSum.Period2StdErr(type)  =  std(brIntChkptPat{:, {colname2}}) / (npat ^ 0.5);
brIntChkptSum.pVal(type)           = pval;
