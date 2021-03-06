clc; clear; close all;

tic
studynbr = 4;
study = 'BR';
studyfullname = 'Breathe';
basedir = setBaseDir();
subfolder = 'MatlabSavedVariables';
[datamatfile, clinicalmatfile, demographicsmatfile] = getRawDataFilenamesForStudy(study);
[physdata, offset] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);
[cdPatient, cdDrugTherapy, cdMicrobiology, cdAntibiotics, cdAdmissions, cdPFT, cdCRP, ...
    cdClinicVisits, cdOtherVisits, cdEndStudy, cdHghtWght] = loadAndHarmoniseClinVars(clinicalmatfile, subfolder, study);
toc

tic
% get the date scaling offset for each patient
patientoffsets = getPatientOffsets(physdata);

% extract clinical FEV1 measures and join with offsets to keep only those patients who
% have enough data (ie the patients left after outlier date handling
pclinicalfev = sortrows(cdPFT(:,{'ID', 'LungFunctionDate', 'FEV1'}), {'ID', 'LungFunctionDate'}, 'ascend');
pclinicalfev.Properties.VariableNames{'ID'} = 'SmartCareID';
pclinicalfev = innerjoin(pclinicalfev, patientoffsets);

% create a scaleddatenum to translate the study date to the same normalised
% scale as measurement data scaled date num
pclinicalfev.ScaledDateNum = datenum(pclinicalfev.LungFunctionDate) - offset - pclinicalfev.PatientOffset;

% extract study date and join with offsets to keep only those patients who
% have enough data (ie the patients left after outlier date handling
pstudydate = sortrows(cdPatient(:,{'ID', 'Hospital', 'StudyNumber', 'StudyDate'}), 'ID', 'ascend');
pstudydate.Properties.VariableNames{'ID'} = 'SmartCareID';
pstudydate = innerjoin(patientoffsets, pstudydate);

% create a scaleddatenum to translate the study date to the same normalised
% scale as measurement data scaled date num
pstudydate.ScaledDateNum = datenum(pstudydate.StudyDate) - offset - pstudydate.PatientOffset + 1;


% extract just the weight measures from smartcare data
pmeasuresfev = physdata(ismember(physdata.RecordingType,'FEV1Recording'),{'SmartCareID', 'ScaledDateNum', 'FEV'});
pmeasuresfev.Properties.VariableNames{'FEV'} = 'FEV1';

% store min and max to scale x-axis of plot display. Set min to -5 if less
% than, to avoid wasting plot space for the one patient with a larger delay
% between study date and active measurement period
mindays = min([pmeasuresfev.ScaledDateNum ; pstudydate.ScaledDateNum]);
if mindays < -5
    mindays = -5;
end
maxdays = max([pmeasuresfev.ScaledDateNum ; pstudydate.ScaledDateNum + 183]);

plotsacross = 3;
plotsdown = 5;
plotsperpage = plotsacross * plotsdown;

subfolder = sprintf('Plots/%s', study);
if ~exist(strcat(basedir, subfolder), 'dir')
    mkdir(strcat(basedir, subfolder));
end
toc

tic
% create plots for patients with differences home vs clinical
%fprintf('FEV Plots for diff values home vs clinical\n');
%patientlist = [82 ; 99 ; 175 ; 194 ; 200];
%filenameprefix = 'CalcClinicalVsHomeFEV1 - Different Values';
%figurearray = createAndSaveFEVPlots(patientlist, pmeasuresfev, pclinicalfev, pstudydate, ...
%    mindays, maxdays, plotsacross, plotsdown, plotsperpage, basedir, subfolder, filenameprefix);
%close all;
%toc
%tic
% create plots for patients with step function change in home measures
%fprintf('FEV Plots for step function change in home measuresl\n');
%patientlist = [54 ; 81 ; 138];
%filenameprefix = 'CalcClinicalVsHomeFEV1 - Step Function Change';
%figurearray = createAndSaveFEVPlots(patientlist, pmeasuresfev, pclinicalfev, pstudydate, ...
%    mindays, maxdays, plotsacross, plotsdown, plotsperpage, basedir, subfolder, filenameprefix);
%close all;
%toc
%tic
% create plots for patients with incorrect height
%fprintf('FEV Plots for incorrect height\n');
%patientlist = [179];
%filenameprefix = 'CalcClinicalVsHomeFEV1 - Incorrect Height';
%figurearray = createAndSaveFEVPlots(patientlist, pmeasuresfev, pclinicalfev, pstudydate, ...
%    mindays, maxdays, plotsacross, plotsdown, plotsperpage, basedir, subfolder, filenameprefix);
%close all;
%toc
%tic
% create plots for potential anomalous clinical FEV1 measures identified
%fprintf('FEV Plots for potential anomalous clinical measures\n');
%patientlist = [24 ; 66];
%filenameprefix = 'CalcClinicalVsHomeFEV1 - Outlier Clinical Values';
%figurearray = createAndSaveFEVPlots(patientlist, pmeasuresfev, pclinicalfev, pstudydate, ...
%    mindays, maxdays, plotsacross, plotsdown, plotsperpage, basedir, subfolder, filenameprefix);
%close all;
%toc
%tic
% create plots for potential anomalous Home FEV1 measures identified
%fprintf('FEV Plots for potential anomalous Home measures\n');
%patientlist = [46 ; 53 ; 78 ; 171];
%filenameprefix = 'CalcClinicalVsHomeFEV1 - Outlier Home Values';
%figurearray = createAndSaveFEVPlots(patientlist, pmeasuresfev, pclinicalfev, pstudydate, ...
%    mindays, maxdays, plotsacross, plotsdown, plotsperpage, basedir, subfolder, filenameprefix);
%close all;
%toc

tic
% create plots for all patients
fprintf('FEV Plots for all patients\n');
patientlist = unique(pmeasuresfev.SmartCareID);
filenameprefix = sprintf('%s-CalcClinicalVsHomeFEV1', study);
createAndSaveBreatheFEVPlots(patientlist, pmeasuresfev, pclinicalfev, pstudydate, ...
    mindays, maxdays, plotsacross, plotsdown, plotsperpage, subfolder, filenameprefix);

toc
