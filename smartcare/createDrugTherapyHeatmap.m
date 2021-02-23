function createDrugTherapyHeatmap(brPatient, brDrugTherapy, filter)
% set various parameters
study = 'BR';

if nargin<2
    filter = 0;
    subfolder = 'MatlabSavedVariables';
    clinicalmatfile = 'breatheclinicaldata';
    [brPatient, brDrugTherapy, ~, ~, ~, ~, ~, ...
    ~, ~, ~, ~, ~, ~, ~] = loadAndHarmoniseClinVars(clinicalmatfile, subfolder, study);
end

startd = datetime(2019, 02, 04);
cutoffd = datetime(2020, 11, 30); % cutoff date is last date the data was processed
mintwindow = 6; % minimum time window for analysis

%might want to do some filtering here - eg. filter study start (129->80) - ok but maybe don't need
if filter==1
    brPrePostPat   = brPatient;
    [brPrePostPat] = filterPrePostByStudyStart(brPrePostPat, mintwindow);
    a = unique(brDrugTherapy.ID(ismember(brDrugTherapy.ID, brPrePostPat.ID)));
    brDrugTherapy(~ismember(brDrugTherapy.ID, a),:) = [];
end

%set up mapping container
keySet = {'iva','sym','ork','tri'};
M = containers.Map(keySet,[1 2 3 4]);

%rename DrugTherapyType to the common ones
a = lower(brDrugTherapy.DrugTherapyType);
a(contains(a,'modulator')) = {'tri'};
a(contains(a,'iva')) = {'iva'};
a(contains(a,'sym')) = {'sym'};
a(contains(a,'ork')) = {'ork'};
a(contains(a,'tri')) = {'tri'};
b = size(find(~ismember(a, keySet)),1);
if b~=0 fprintf('%d Unidentified drug therapy',b), end

[d, id] = findgroups(a); counts = histcounts(d);
for i=1:size(id,1); fprintf('%s: %d   ',id{i},counts(i)); end

brDrugTherapy.DrugTherapyType = a;

%check comment
%brDrugTherapy(contains(brDrugTherapy.DrugTherapyComment,'date'),:) = [];
a = brDrugTherapy.ID(not(cellfun('isempty', brDrugTherapy.DrugTherapyComment)));
cmt = unique(a);

%scale DrugTherapyStartDate
ref = datenum(startd); %min(datenum(a))
dtd = datenum(brDrugTherapy.DrugTherapyStartDate) - ref + 1;
dtd(dtd<1)=1;
cld = datenum(brPatient.PatClinDate) - ref + 1;
std = datenum(brPatient.StudyDate) - ref + 1;

%generate mat
%0 else/not_on_DT NaN   1 iva   2 sym   3 ork   4 tri   5 overallstart   6 patstart   7 patclindate
ptID = unique(brDrugTherapy.ID);
ndays = datenum(cutoffd)-datenum(startd)+1;
dtArray = zeros(size(ptID,1),ndays+20);
for i=1:size(ptID,1)
    if ismember(ptID(i),cmt)
        dtArray(i,8)=8;
    end
    nd = find(arrayfun(@(x)isequal(x,ptID(i)),brDrugTherapy.ID))';
    ns = find(arrayfun(@(x)isequal(x,ptID(i)),brPatient.ID))';
    for n=nd
        dtArray(i,dtd(n)+10:end-10)=M(char(brDrugTherapy.DrugTherapyType(n)));
    end

    dtArray(i, 1+10) = 5;
    dtArray(i, std(ns(1))+10) = 6;
    dtArray(i, cld(ns(1))+10) = 7;
end

%heatmap
colors = [37 37 37; 239 243 255; 189 215 231; 107 174 214; 33 113 181; 252 174 145; 251 106 74; 203 24 29; 254 204 92]/256;
title = sprintf('%s-Heatmap of Drug Therapy Period', study);
[f, p] = createFigureAndPanel(title, 'portrait', 'a4');
h = heatmap(p, dtArray, 'Colormap', colors);
h.XLabel = 'Days since study';
h.YLabel = 'Patients';
h.CellLabelColor = 'none';
h.GridVisible = 'off';

%imagesc
img = imagesc(dtArray);
colormap(colors);
c = colorbar('Ticks',linspace(0,8,9),...
         'TickLabels',{'Not on DT', 'Ivacaftor', 'Symkevi', 'Orkambi', 'Trikafta', 'Study start', 'Pt start', 'Pt clinic','Comment'});

basedir = setBaseDir();
filename = sprintf('%s-HeatmapAllPatientsWithStudyPeriod', study);
subfolder = sprintf('Plots/%s', study);
if ~exist(strcat(basedir, subfolder), 'dir')
    mkdir(strcat(basedir, subfolder));
end
savePlotInDir(f, filename, subfolder);
close(f);

end
