totalTimerId = tic;

% Add tensor toolbox
% TODO Improve so that no hardcoded urls are required!!
try
   here = pwd; cd 'C:\Users\rleira\Dropbox\08-Proyecto ANII\matlabLibs\tensor_toolbox';addpath(pwd);cd(here);
catch ME
   here = pwd; cd 'D:\Dropbox\08-Proyecto ANII\matlabLibs\tensor_toolbox';addpath(pwd);cd(here);
   %here = pwd; cd '~/Dropbox/08-Proyecto ANII/matlabLibs/tensor_toolbox';addpath(pwd);cd(here); %mac
end

format long g

% Add path to utils so they can be used here
mfilepath = fileparts(which(mfilename));
addpath(fullfile(mfilepath,'../MatlabOptimizationUtils'));

disp('LOADING DATA...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------PARAMETERS -------------
format long g

if ~exist('cantRepeat', 'var')
   cantRepeat = 20;
end

if ~exist('lumCount', 'var')
   lumCount = 3;
end

if ~exist('CargarMatrices', 'var')
   CargarMatrices=1;
end

%objectiveImage = '../GeneralData/surfacesImages/piso_patio_fino2Pisos21824.bmp';
objectiveImage = '../GeneralData/surfacesImages/paredObjetivo.bmp';
if ~exist('positionsImage', 'var')
    positionsImage = '../GeneralData/surfacesImages/techo_flotante_no_black.bmp';
end
% -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------LOAD DATA---------------
if CargarMatrices
    disp('LOADING LUMINARIES MATRICES...');
    
    %load([fullfile(mfilepath,'../GeneralData/mats/TopBottom/Luminaires/') 'luminaireMatrix.mat'], 'luminaireMatrix.mat');%load('luminaireHemicubeMatrix.mat', 'luminaireHemicubeMatrix');
    load([fullfile(mfilepath,'../GeneralData/mats/TopBottom/Luminaires/') 'luminaireNameMap.mat']);%load('luminaireNameMap.mat', 'luminaireNameMap');
    load([fullfile(mfilepath,'../GeneralData/mats/TopBottom/Luminaires/') 'luminairePower.mat']);%load('luminairePower.mat', 'luminairePower');
    load([fullfile(mfilepath,'../GeneralData/mats/TopBottom/Luminaires/') 'orderedLuminarieListIndexes.mat']);%load('orderedLuminarieListIndexes.mat', 'orderedLuminarieListIndexes');
    
    load([fullfile(mfilepath,'../GeneralData/mats/TopBottom/PositionsNumero8/') 'positionsMatrix.mat'], 'positionsMatrix');
    load([fullfile(mfilepath,'../GeneralData/mats/TopBottom/PositionsNumero8/') 'positionsNameMap.mat'], 'positionsNameMap');
    %load([fullfile(mfilepath,'../GeneralData/mats/TopBottom/Positions/') 'positionsMatrix.mat'], 'positionsMatrix');
    %load([fullfile(mfilepath,'../GeneralData/mats/TopBottom/Positions/') 'positionsNameMap.mat'], 'positionsNameMap');
    %load([fullfile(mfilepath,'../GeneralData/CentroTechoFlotante/OptimizationEngine2/') 'positionsMatrix.mat'], 'positionsMatrix');%load('positionsMatrix.mat', 'positionsMatrix');
    %load([fullfile(mfilepath,'../GeneralData/CentroTechoFlotante/OptimizationEngine2/') 'positionsNameMap.mat'], 'positionsNameMap');%load('positionsNameMap.mat', 'positionsNameMap');
    
    % ORDER AUXILIAR MATRICES
    % Order luminaireNameMap to make it consistent with the luminaireHemicubeMatrix
    orderedLuminaireNameMap = luminaireNameMap(orderedLuminarieListIndexes);
    % Order luminairePower to make it consistent with the luminaireHemicubeMatrix
    luminairePower = luminairePower(orderedLuminarieListIndexes);
    maxPower = max(luminairePower);
    %Order luminaireNameMap with matrix order
    luminaireNameMap = luminaireNameMap(orderedLuminarieListIndexes);
    
    disp('LUMINARIES MATRICES LOADED SUCCESSFULLY...');
else
    disp('LUMINARIES MATRICES SEEM TO BE ALREADY LOADED, skipping step...');
end

% -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------CALCULATE---------------
% Get the triangles on the source surface
%[positionsMapMatrix selectedTriangles selectedIndexesOrdered] = loadViewSurface('../GeneralData/surfacesImages/techo_flotante_no_black.bmp', 'BMP', 70, 21824);
[positionsMapMatrix selectedTriangles selectedIndexesOrdered] = loadViewSurface('../GeneralData/surfacesImages/JustIndexes/numero8_9088/numero8_techo_flotante.bmp', 'BMP', 70, 9088);

% Get the triangles on the objective
%[objectivePositionsMapMatrix objectiveSelectedTriangles objectiveSelectedIndexesOrdered] = loadViewSurface(objectiveImage, 'BMP', 1, 21824);
objectiveSelectedIndexesOrdered = load('../GeneralData/surfacesImages/JustIndexes/numero8_9088/selectedIndexesOrdered_piso');
objectiveSelectedIndexesOrdered = objectiveSelectedIndexesOrdered.selectedIndexesOrdered;

% Transform positions file names to indexes
positionsMap = regexp(positionsNameMap, 'H_(\d+).bmp', 'tokens', 'once');
positionsMap = vertcat(positionsMap{:});
positionsMap = sprintf('%s ', positionsMap{:});
selectedIndexesOrdered = sscanf(positionsMap, '%u') + 1;

prefijo = 'numero8'; %prefijo = 'patio_fino2Pisos21824';   %es el prefijo de los archivos que contienen U, V, R, ?rea.

%<<<<<<<<< PAR?METROS DE EJECUCI?N <<<<<<<<<<<<<<<<<<

nvars = 6;            % count of variables
allVars = nvars * lumCount; % count of all vars on each solution
if ~exist('MaxNroIterIf', 'var')
   MaxNroIterIf = 25000; %cu?ntos c?lculos de radiosidad realiza? 
end
cppit = 10;

deltavalues = [0.1 0.15 0.3 0.5 0.7 1.0]; 
if lumCount == 1
    %lmstepvalues = 1:nvars;
    lmstepvalues = 1:nvars;
else
    %lmstepvalues = 1:nvars:allVars;
    lmstepvalues = 1:3:allVars;
end

[aux maxdeltavalues] = size(deltavalues);
[aux maxlmstepvalues] = size(lmstepvalues);

[y x] = size(positionsMapMatrix);
[l aux] = size(orderedLuminarieListIndexes);

% Calculate the max value a neighbor can have
variableMaxNeighbor = [y;x;l;360;180;360];

if ~exist('oneEach', 'var')
   oneEach = [1;1;1;45;45;90];
end

if ~exist('oneEachMinVal', 'var')
   variableMinVal = [1;1;1;0;0;0];
end

variableMaxNeighborLumCount = [];
oneEachLumCount = [];
variableMinValLumCount = [];
for i=1:lumCount
    variableMaxNeighborLumCount = [variableMaxNeighborLumCount variableMaxNeighbor'];
    oneEachLumCount = [oneEachLumCount;oneEach];
    variableMinValLumCount = [variableMinValLumCount;variableMinVal];
end
variableMaxNeighborLumCountT = variableMaxNeighborLumCount';

% Fix to make tensor not use zeros as indexes (we later sum the min value for it)
variableTensorFix = variableMinValLumCount';
nonZeroIndexes = find(variableTensorFix);
zeroIndexes = find(~variableTensorFix);
variableTensorFix(nonZeroIndexes) = 0; % if values are nonzero then we don sum
variableTensorFix(zeroIndexes) = 1; % if values are nonzero then we don sum

if ~exist('weightArray', 'var')
   weightArray = [1 1 10 8 4 4]
   %weightArray = [1 1 3 1 1 1]
end



weightArray = [ones(1, weightArray(1)) ones(1, weightArray(2))*2 ones(1, weightArray(3))*3 ones(1, weightArray(4))*4 ones(1, weightArray(5))*5 ones(1, weightArray(6))*6];
[aux weightArraySize] = size(weightArray);

weightArrayAll = [];
for i=1:nvars:allVars
    weightArrayAll = [weightArrayAll (weightArray+(i-1))];
end

[aux weightArrayAllSize] = size(weightArrayAll);

%target = load('../GeneralData/targets/WithRots/target_1Lum_920_6_180_90_180_ExpV.mat');

%target = load('../GeneralData/surfacesImages/target.mat');
%target = load('../GeneralData/targets/WithRots/target_1Lum_667_18_360_180_90.mat');
target = load('../GeneralData/targets/WithRots/target_2Lum_315_921_20_35_180_45_90_135_360_180_X.mat'); % Test Eduardo Girona, cambiar objetivo y ver si mejora
%target = load('../GeneralData/targets/WithRots/target_2Lum_1149_436_13_27_225_315_180_45_360_90.mat');
%target = load('../GeneralData/targets/WithRots/target_3Lum_688_938_66_11_46_13_180_45_90_180_45_135_90_270_360.mat');
target = target.B;
% target = load('../GeneralData/targets/focal1Lum/B.mat');
% target = target.B;
% target = target(:, 1);

if CargarMatrices
    %here = pwd; cd '..\OptimizationEngine\hemicubesGenerator\';
    %here = pwd; cd '../OptimizationEngine/hemicubesGenerator/';
    here = pwd; cd '../GeneralData/lowRankData/numero8_9088/';
    
    disp('LOADING LOW-RANK RADIOSITY MATRICES...');
    CargarMatrices=0;
    V=single(load([prefijo '_V.dat'],'-ascii'));
    [k,n]=size(V);
    VI=zeros(1,n);
    for i=1:n,VI(i)=find(V(:,i)==1);end
    U=single(load([prefijo '_U.dat'],'-ascii'));
    R=(load([prefijo '_ros.dat'],'-ascii'));
    RU=zeros(n,k,'single');
    for i=1:n, RU(i,:)=U(i,:)*sum(R(i,:))/3;end 

    VRU=zeros(k,'single');
    for i=1:n, VRU(VI(i),:)=VRU(VI(i),:)+RU(i,:); end
    Y=-RU*inv(eye(k)-(VRU));
    clear U;
    area=single(load([prefijo '_area.dat'],'-ascii'));
    Yp = Y(objectiveSelectedIndexesOrdered, :);
    Rmean = mean(R')';
    disp('LOW-RANK RADIOSITY MATRICES LOADED SUCCESSFULLY...');
    cd(here);
else
    disp('LOW-RANK RADIOSITY MATRICES SEEM TO BE ALREADY LOADED, skipping step...');
end
triangleCount = n;

% Calculate cx and cgamma fixed hemicubes, top and bottom
[Hcx, Hcgamma, TopHcx, TopHcgamma, pixelDeltaAngle] = generateCXAndCGammaHemicubes;
[xV, yV, zV] = sph2cart(deg2rad(hemicubeToVector(Hcx)), deg2rad(hemicubeToVector(Hcgamma) - 90), 1);
bottomHemicubeCXAndCGammaAsCartesian = [xV';yV';zV'];
[xV, yV, zV] = sph2cart(deg2rad(hemicubeToVector(TopHcx)), deg2rad(hemicubeToVector(TopHcgamma) - 90), 1);
topHemicubeCXAndCGammaAsCartesian = [xV';yV';zV'];

matlabFilesPath = fullfile(mfilepath,'../GeneralData/matlab/');
matlabFileNamePostFix = '_matlab2.mat';
% -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%REORGANIZE VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem specific
problemObj.objectiveSelectedIndexesOrdered = objectiveSelectedIndexesOrdered;
problemObj.selectedIndexesOrdered = selectedIndexesOrdered;
problemObj.target = target;
problemObj.triangleCount = triangleCount;
problemObj.lumCount = lumCount;
problemObj.nvars = nvars;
problemObj.allVars = allVars;
problemObj.MaxNroIterIf = MaxNroIterIf;
problemObj.cppit = cppit;
problemObj.luminarieLocation = '../OptimizationEngine/hemicubesGenerator/lumOneByOne/';
problemObj.luminarieMatlabsLocation = '../GeneralData/matlab/';
problemObj.weightArray = weightArray;
problemObj.oneEach = oneEach;
problemObj.oneEachLumCount = oneEachLumCount;
problemObj.variableMinValLumCount = variableMinValLumCount;
problemObj.variableTensorFix = variableTensorFix;
problemObj.weightArrayAllSize = weightArrayAllSize;
problemObj.objectiveImage = objectiveImage;
problemObj.positionsImage = positionsImage;
problemObj.deltavalues = deltavalues;
problemObj.maxdeltavalues = maxdeltavalues;
problemObj.lmstepvalues = lmstepvalues;
problemObj.maxlmstepvalues = maxlmstepvalues;
problemObj.hemicubeSize = 512;
problemObj.variableMaxNeighborLumCount = variableMaxNeighborLumCount;
problemObj.variableMaxNeighbor = variableMaxNeighbor;

%Restrictions
restrictionObj.maxPower = 100;

% Hemicube data
hemicube.Hcx = Hcx;
hemicube.Hcgamma = Hcgamma;
hemicube.TopHcx = TopHcx;
hemicube.TopHcgamma = TopHcgamma;

% Hemicube data as vector
hemicubeAsVector.bottomHemicubeCXAndCGammaAsCartesian = bottomHemicubeCXAndCGammaAsCartesian;
hemicubeAsVector.topHemicubeCXAndCGammaAsCartesian = topHemicubeCXAndCGammaAsCartesian;
hemicubeAsVector.pixelDeltaAngle = hemicubeToVector(pixelDeltaAngle)';
hemicubeAsVector.pixelDeltaAngle = [hemicubeAsVector.pixelDeltaAngle hemicubeAsVector.pixelDeltaAngle];

positions.hemicubesAsVectors = positionsMatrix;
positions.positionsNameMap = positionsNameMap;
positions.black = 16777216;


% Low-Rank specific
lowRankObj.Rmean = Rmean;
lowRankObj.RmeanMatrix = diag(Rmean);
lowRankObj.Yp = Yp;
lowRankObj.Y = Y;
lowRankObj.V = V;
lowRankObj.areaPatches = area;
lowRankObj.areaSelectedPatches = lowRankObj.areaPatches(problemObj.objectiveSelectedIndexesOrdered);
lowRankObj.sumAreaSelectedPatches = sum(lowRankObj.areaSelectedPatches);

% Positions specific
positionObj.positionsMapMatrix = positionsMapMatrix;
positionObj.positionsMatrix = positionsMatrix;

% Luminarie specific
luminarieObj.orderedLuminarieListIndexes = orderedLuminarieListIndexes;
luminarieObj.orderedLuminaireNameMap = orderedLuminaireNameMap;
luminarieObj.maxPower = maxPower;
luminarieObj.luminairePower = luminairePower;
% - Luminaire files path
luminarieObj.pathToMatlab2Files = matlabFilesPath;
% - Luminaire files postfix
luminarieObj.matlabFileNamePostFix = matlabFileNamePostFix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ALL DATA LOADED SUCCESSFULLY...');
totalTime = toc(totalTimerId);
X = ['TIME CONSUMED LOADING DATA: ', num2str(totalTime), ' seconds'];
disp(X);
