% positionsPath - path containing the position hemicube images
% luminarePath  - path containing the luminarie .matlab files
%function [emision selectedLuminary positionsAux res totalTime] = optimizeLuminariesN(cantIter, problemObj, restrictionObj, luminarieObj)

res = ones(1, 4 + (problemObj.nvars - 1) * problemObj.lumCount);

if ~exist('cantRepeat', 'var')
    cantRepeat = 1;    %cantidad de veces que se desea ejecutar
end

if exist('cantIter', 'var')
    cantRepeat = cantIter;
end

areaPatch = (0.3*0.3)/2; % triangles of half a meter

MENOR = Inf;   %%% Esta variable contiene el mejor valor de la funci?n a optimizar (hasta el momento).
NroIterIf = 0;

[y x] = size(positionObj.positionsMapMatrix);
[l aux] = size(luminarieObj.orderedLuminarieListIndexes);

MaxNroIterIf = min(problemObj.MaxNroIterIf, prod(problemObj.variableMaxNeighbor));

problemObj.variableMaxNeighborLumCount = [];
for i=1:problemObj.lumCount
    problemObj.variableMaxNeighborLumCount = [problemObj.variableMaxNeighborLumCount problemObj.variableMaxNeighbor'];
end
variableMaxNeighborLumCountT = problemObj.variableMaxNeighborLumCount';

Nsolu=0;

MENOR_REP = ones(cantRepeat, 1) * Inf;
tic;
totalTimerId = tic;
for ii=1:cantRepeat
    
    disp(['iteracion ' int2str(ii)])

    % Generate start solution, when lum count = 1
    if (problemObj.lumCount == 1)
        VExtrOpt = (floor(rand(problemObj.nvars, 1).*(problemObj.variableMaxNeighbor./problemObj.oneEach)) + 1).*problemObj.oneEach;
        variableMaxNeighborForNLum = problemObj.variableMaxNeighbor;
    else
        VExtrOpt = [];
        variableMaxNeighborForNLum = [];
        % Generate start solution, when lum count > 1
        for lum=1:problemObj.lumCount
            VExtrOpt = [VExtrOpt; ((floor(rand(problemObj.nvars, 1).*(problemObj.variableMaxNeighbor./problemObj.oneEach)) + 1).*problemObj.oneEach)];
            variableMaxNeighborForNLum = [variableMaxNeighborForNLum; problemObj.variableMaxNeighbor];
        end
    end
    
    % Generate tabu matrix
    tabu = sptensor(problemObj.variableMaxNeighborLumCount);
    tabu(VExtrOpt') = 1;tabu(VExtrOpt') = 0; % HORRIBLE HACK!!!! THIS is made so that sptensor is correctly initialized and ready to use after this point!!
    tabuFilled = 0;
    % Generate indexes to travel neighborhoods
    iii=1;jjj=1;nextjjj = 1;
    
    MENOR=Inf; %%% Esta variable contiene el mejor valor de la funci?n a optimizar (hasta el momento).%76;
    NroIter=0; NroIterIf=0; 


    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ProbarVariasVeces=0;
    h=1; 
    delta = problemObj.deltavalues(iii);
    lmstep = problemObj.lmstepvalues(jjj);
    
    singleTimerId = tic;
    bestIt = 1;
    
    timesInvalidPosition = 0;
    while NroIterIf < MaxNroIterIf && ~tabuFilled
        %tic; %delete me
        
        NroIter=NroIter+1;
        
        % Get the current solution
        VExtr = VExtrOpt;
        
%-------------------------------------------
% GET RANDOM NEW SOLUTION WITH NEIGHBORGHOOD
%-------------------------------------------

        % Get a random sort of the variables
        [aux, sortedRandIndexes] = sort(rand(1, problemObj.weightArrayAllSize));
        randomWeightArrayAll = weightArrayAll(sortedRandIndexes);
        % https://www.mathworks.com/matlabcentral/newsreader/view_thread/263654
        [aux, v] = unique(randomWeightArrayAll);
        v = randomWeightArrayAll(sort(v));
       
        %Get a random value for the selected variables (chosen randomly)
        % new value is the existing value +/- the neighborgh size that is 
        %calculted by multipling the max neighborgh size by delta
        random = rand(problemObj.allVars, 1);
        plusMinus = sign(rand(problemObj.allVars, 1)-0.5);
        
        %VExtr(v(1:lmstep)) = VExtr(v(1:lmstep)) + plusMinus(1:lmstep).*(floor((variableMaxNeighborLumCountT(v(1:lmstep)).*problemObj.oneEach(v(1:lmstep))).*random((1:lmstep))*delta) + 1);
        selectedVars = v(1:lmstep);
        
        variableMaxNeighborLumCountTOnRange = variableMaxNeighborLumCountT(selectedVars)./problemObj.oneEachLumCount(selectedVars);
        radius = floor(variableMaxNeighborLumCountTOnRange*delta) + 1;
        VExtrOnRange = VExtr(selectedVars)./problemObj.oneEachLumCount(selectedVars);
        s_minus_r = VExtrOnRange - radius;
        s_plus_r = VExtrOnRange + radius;

        right = min(s_plus_r, variableMaxNeighborLumCountTOnRange);
        left = max(s_minus_r, problemObj.variableMinValLumCount(selectedVars));

        VExtr(selectedVars) = round((right - left).*random(1:lmstep) + left).*problemObj.oneEachLumCount(selectedVars);
        
        %VExtr(selectedVars) = VExtr(selectedVars) + (plusMinus(1:lmstep).*(floor((variableMaxNeighborLumCountT(selectedVars)./problemObj.oneEachLumCount(selectedVars)).*random((1:lmstep))*delta) + 1)).*problemObj.oneEachLumCount(selectedVars);
        
%------------------------------------------- 

%------------------ DO NOT ALLOW INVALID POSITIONS ------------------------
% this is so ugly...
        invalidPosition = false;
        for VExtrIndex=1:problemObj.lumCount
            lowerLimit = (((VExtrIndex - 1) * problemObj.nvars ) + 1);
            VExtrAux = VExtr(lowerLimit:lowerLimit + (problemObj.nvars - 1));
            % Calculate the emision merging both the view hemicube with the
            % radiance hemicube taken from the luminaire
            position = positionObj.positionsMapMatrix(VExtrAux(1), VExtrAux(2));
            invalidPosition = position == positions.black;
            if (invalidPosition) 
                break;
            end
        end
        if (invalidPosition)
            timesInvalidPosition = timesInvalidPosition + 1;
            if (timesInvalidPosition>problemObj.cppit*(lmstep+1)*(delta/problemObj.deltavalues(1)))
                % Move to next neighborghood
                [iii jjj nextjjj] = nextMatrixPoint(iii, jjj, nextjjj, problemObj.maxdeltavalues, problemObj.maxlmstepvalues);
                delta = deltavalues(iii);
                lmstep = problemObj.lmstepvalues(jjj);
                h = 1;
                timesInvalidPosition = 0;
            end
            continue;
        end
        timesInvalidPosition = 0;
%--------------------------------------------------------------------------

        % New variables must be inside the [problemObj.variableMinValLumCount problemObj.variableMaxNeighbor] range
        if (min(VExtr-problemObj.variableMinValLumCount)>=0) && (min(variableMaxNeighborForNLum-VExtr)>=0)
               
            if (ProbarVariasVeces)
                if (h > problemObj.cppit)
                    iii = 1;
                    jjj = 1;
                    nextjjj = 1;
                    ProbarVariasVeces = 0;
                    lmstep = problemObj.lmstepvalues(jjj);
                    delta = problemObj.deltavalues(iii);
                    h = 0;
                end
            end

            h=h+1;
            if ~(ProbarVariasVeces) && h>problemObj.cppit*(lmstep+1)*(delta/problemObj.deltavalues(1)), 
                %En diagonal, siempre intenta ir abajo a la izquierda si no
                %puede vuelve a iniciar en el primer elemento disponible de
                %la primer columna que no fue consultada (o la ultima
                %columna de la matriz en caso que no hayan mas)
                [iii jjj nextjjj] = nextMatrixPoint(iii, jjj, nextjjj, problemObj.maxdeltavalues, problemObj.maxlmstepvalues);
                delta = deltavalues(iii);
                lmstep = problemObj.lmstepvalues(jjj);
                %[delta, lmstep,h]
                h = 1;
            end
            
            if (~tabu(VExtr'))
                tabu(VExtr'+problemObj.variableTensorFix) = 1;
                
                %%%%%%%%%% C?lculo de E %%%%%%%%%%%%%%
                E = zeros(n,1);
%-------------------------------------------
% GENERATE E FOR EACH SUB SET OF SOLUTIONS
%-------------------------------------------
                % REMOVE ME
                %VExtr = [1;1;260;180;45;90]; %original
                %VExtr = [1;1;1512;180;45;90]; %original
                %VExtr = [1;1;1512;90;45;90]; % rotado
                %VExtr = [1;1;1508;0;0;0]; % rotado
                %VExtr = [2;6;1422;135;90;90]; %[1422 37 135 90 90] pos 37
                %VExtr = [3;10;1334;225;90;90;12;1;1347;45;90;90]; %Unif Num8
                %VExtr = [9;1;1424;45;90;90;3;10;456;225;90;90;15;9;1377;315;90;90]; %Unif Num8
                % REMOVE ME end

                %VExtr = [1;1;1325;339;139;330;1;1;1517;188;85;139;1;1;619;196;32;343];
                %VExtr = [1;1;1517;339;139;330;1;1;1517;188;85;139;1;1;619;196;32;343];
                for VExtrIndex=1:problemObj.lumCount
                    lowerLimit = (((VExtrIndex - 1) * problemObj.nvars ) + 1);
                    VExtrAux = VExtr(lowerLimit:lowerLimit + (problemObj.nvars - 1));
                    
                    % Calculate the emision merging both the view hemicube with the
                    % radiance hemicube taken from the luminaire
                    position = find(problemObj.selectedIndexesOrdered==positionObj.positionsMapMatrix(VExtrAux(1), VExtrAux(2)));
                    
                    % REMOVE ME
                    %position = 1;
                    % REMOVE ME end
                    
                    fileName = [problemObj.luminarieMatlabsLocation char(luminarieObj.orderedLuminaireNameMap{VExtrAux(3)}) '_matlab2.mat'];
                    if strcmp(char(luminarieObj.orderedLuminaireNameMap{VExtrAux(3)}), 'EMPTY')
                        luminPolarCurveMatrix = zeros(361, 181);
                    else
                        luminPolarCurveMatrix = load(fileName);
                        luminPolarCurveMatrix = luminPolarCurveMatrix.matlabVersion2;
                    end
                    
                    z1 = VExtrAux(4);
                    y2 = VExtrAux(5);
                    z3 = VExtrAux(6);
                                      
                    lumE = calculateFirstReflectionWithHemicubesAsVectors(...
                        positionObj.positionsMatrix(position, :), ...
                        luminPolarCurveMatrix, ...
                        problemObj.hemicubeSize, ...
                        hemicubeAsVector.bottomHemicubeCXAndCGammaAsCartesian, ...
                        hemicubeAsVector.topHemicubeCXAndCGammaAsCartesian, ...
                        hemicubeAsVector.pixelDeltaAngle, ...
                        z1, ...
                        y2, ...
                        z3, ...
                        problemObj.triangleCount ...
                    );
                    
                    E = E + lumE(1:problemObj.triangleCount);
                    
                    positionsAux(VExtrIndex) = position;
                    selectedLuminary(VExtrIndex) = VExtrAux(3);
                    z1Angles(VExtrIndex) = z1;
                    y2Angles(VExtrIndex) = y2;
                    z3Angles(VExtrIndex) = z3;
                end
%-------------------------------------------

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                %<<<<<<<<< OBJETIVO A MINIMIZAR <<<<<<<<<<
                %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                % TODO no minimization objective
                E = E(1:problemObj.triangleCount);
                E = (E.*lowRankObj.Rmean)./lowRankObj.areaPatches;
                NroIterIf=NroIterIf+1;
                
                %---------------------------------------
                % Find Uniformity/Dispersion:
                %---------------------------------------
%                 B = -lowRankObj.Yp*(lowRankObj.V*E) + E(problemObj.objectiveSelectedIndexesOrdered);
%                 BB = zeros(triangleCount, 1);
%                 BB(problemObj.objectiveSelectedIndexesOrdered) = B; 
%                 objetivo = std(B)/mean(B);
                %---------------------------------------
                
                %---------------------------------------
                % Find Weighted Uniformity/Dispersion:
                %---------------------------------------
                B = -lowRankObj.Yp*(lowRankObj.V*E) + E(problemObj.objectiveSelectedIndexesOrdered);
                weightedMean = (sum(B.*lowRankObj.areaSelectedPatches))/lowRankObj.sumAreaSelectedPatches;
                weightedSTD = sqrt((sum(lowRankObj.areaSelectedPatches.*((B-weightedMean).^2)))/lowRankObj.sumAreaSelectedPatches);
                objetivo = weightedSTD/weightedMean;
                %---------------------------------------
                
                %---------------------------------------
                % Optimize Power Find min and max luxes:
                %---------------------------------------
%                 B = -lowRankObj.Yp*(lowRankObj.V*E) + E(problemObj.objectiveSelectedIndexesOrdered);
%                 BB = zeros(problemObj.triangleCount, 1);
%                 BB(problemObj.objectiveSelectedIndexesOrdered) = B;
%                 objetivo = sum(luminarieObj.luminairePower(selectedLuminary));
%                 % Penalty
%                 c1 = ((max(50 - (min(B)/areaPatch), 0)^2)*100);
%                 if (c1 > 0)
%                     c1 = c1 + luminarieObj.maxPower * problemObj.lumCount; % c1 >= luminarieObj.maxPower
%                 end
%                 c2 = ((max((max(B)/areaPatch) - 750, 0)^2)*100);
%                 if (c2 > 0)
%                     c2 = c2 + luminarieObj.maxPower * problemObj.lumCount; % c2 >= luminarieObj.maxPower
%                 end
%                 
%                 objetivo = objetivo + c1 + c2;
                %---------------------------------------

                %---------------------------------------
                % Find exact solution:
                %---------------------------------------
%                 B = -lowRankObj.Y*(lowRankObj.V*E) + E;
%                 diffVect = (B - problemObj.target);
%                 objetivo=sqrt(sum(diffVect.*diffVect));
                %---------------------------------------
                
                %---------------------------------------
                % _Maximize_ sum of luxes on a surface:
                % Restrictions: max luminaries power
                %---------------------------------------
%                 B = -lowRankObj.Yp*(lowRankObj.V*E) + E(problemObj.objectiveSelectedIndexesOrdered);
%                 objetivo = sum(B/areaPatch);
%                 % Penalty
%                 c1 = -1*((min(restrictionObj.maxPower - sum(luminarieObj.luminairePower(selectedLuminary)), 0)^2)*100);
%                 if (c1 < 0)
%                     c1 = c1 + luminarieObj.maxPower * problemObj.lumCount; % c1 >= luminarieObj.maxPower
%                 end
%                 objetivo = -1*(c1 + objetivo/(size(problemObj.objectiveSelectedIndexesOrdered,1)));
                %---------------------------------------
                               
                %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    
                if (objetivo<=MENOR)
                    VExtrOpt = VExtr;              
                    if objetivo < MENOR
                        ProbarVariasVeces=1;
                        h=0; 
                        Efound=E;
                        MENOR=objetivo;
                        Nsolu=Nsolu+1;
                        EEE=[Efound Efound Efound];
                        emision = EEE;
                        orderedPos = VExtr(3);
                        [NroIterIf, objetivo, selectedLuminary, positionsAux, z1Angles, y2Angles, z3Angles]
                        
                        bestPositions = positionsAux;
                        bestLuminaires = selectedLuminary;
                        bestIt = NroIterIf;
                        bestAngles = [z1Angles, y2Angles, z3Angles];
                        
                        for VExtrIndex=1:problemObj.lumCount
                            aux = luminarieObj.orderedLuminaireNameMap(selectedLuminary(VExtrIndex));
                            aux{1}
                            positionsNameMap(positionsAux(VExtrIndex))
                        end
                        
                        %---------------------------------------
                        % Optimize Power Find min and max luxes: display mean
                        %---------------------------------------
%                         mean(B)
                        %---------------------------------------
                        
                        
                        % Save results to disc
%                         save EEE EEE -ascii
                        
                        % Save radiosity of scene plus traingle where
                        % luminaire is placed
%                         patchIndexCell = regexp(positionsNameMap(positionsAux), 'H_(\d+).bmp', 'tokens', 'once');
%                         BBB = double(BB);
%                         BBB(str2double(patchIndexCell{1})) = 1;
%                         BBB = [BBB BBB BBB];
%                         save BBB BBB -ascii
                        
                        if (objetivo == 0)
                            tabuFilled = 1;
                        end
                    end
                end
            end
        end
        %toc; %delete me
    end
    singleTimer = toc(singleTimerId);
    %%%%%%%%%%%%%Hallo la soluci?n real%%%%%%%%%%%%%%%%%%%%%%%
    VE=zeros(k,1);
    for i=1:n, VE(VI(i))=VE(VI(i))+Efound(i); end
    Creal=-Y*(VE);
    VE=zeros(k,1);
    %%%%%%%%%%%%%Hall? la soluci?n real%%%%%%%%%%%%%%%%%%%%%%%
    
    beep

    Nsolu=Nsolu+1;
    
    MENORopt=MENOR;
    
    MENOR_REP(ii, 1) = MENOR;

    disp('end iteration: ')
    disp(ii)
    res(ii, :) = [NroIterIf, bestIt, singleTimer, MENOR_REP(ii), bestLuminaires, bestPositions, bestAngles]
end
toc;
totalTime = toc(totalTimerId);

%hist(MENOR_REP)
MENOR_REP
disp('media   mediana   variancia   desviacion estandard   Luminaria   Position')
show = [mean(MENOR_REP) median(MENOR_REP) var(MENOR_REP) std(MENOR_REP) bestLuminaires bestPositions];
disp(show)
NroIterIf
save optimizeLuminariesNResult res

