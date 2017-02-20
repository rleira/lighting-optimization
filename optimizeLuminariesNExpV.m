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

MENOR = Inf;   %%% Esta variable contiene el mejor valor de la función a optimizar (hasta el momento).
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
    
    MENOR=Inf; %%% Esta variable contiene el mejor valor de la función a optimizar (hasta el momento).
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
    
    while NroIterIf < MaxNroIterIf && ~tabuFilled  

        % Get the current solution
        VExtr = VExtrOpt;
        
%-------------------------------------------------
% INIT VARIABLES FOR HOLE NEIGHBORHOOD EXPLORATION
%-------------------------------------------------   
        hMax = floor(problemObj.cppit * (lmstep + 1) * (delta / problemObj.deltavalues(1)));
        Edyn = [];
        variablesForCandidate = [];%zeros(problemObj.nvars, hMax);
        halleMENOR = 0;
        
        % Iterations are incremented with the amount of candidates in the
        % neighborhood
        NroIter = NroIter + hMax;
        
        candidateIndex = 1;

        % Calculate all the possible candidates for the neighborhood
        for h = 1:hMax
            % Iterate looking for plausible solutions until one is found
            foundPlausibleSol = false;
            triesToFound = 0;
            triesSumand = 1;
            % Do not increment if h = 1 cause we want at least one candidate selected
%             if (h <= 1)
%                 triesSumand = 0;
%             end
            
            while ~foundPlausibleSol && (triesToFound < hMax)
                
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
                
%                 % start DELETE
%                 if NroIterIf == 0
%                     VExtr = [4;5;667;360;180;90]; % DELETE ME!!!!!
%                 end
%                 % end DELETE
                
                %VExtr(selectedVars) = VExtr(selectedVars) + (plusMinus(1:lmstep).*(floor((variableMaxNeighborLumCountT(selectedVars)./problemObj.oneEachLumCount(selectedVars)).*random((1:lmstep))*delta) + 1)).*problemObj.oneEachLumCount(selectedVars);
                
                % Check if the solution is within the ranges for each variable
                % New variables must be inside the [problemObj.variableMinValLumCount problemObj.variableMaxNeighbor] range
                if ((min(VExtr-problemObj.variableMinValLumCount)>=0) && (min(variableMaxNeighborForNLum-VExtr)>=0) && ~tabu(VExtr'))
                    foundPlausibleSol = true;
                    tabu(VExtr'+problemObj.variableTensorFix) = 1;
                end
                triesToFound = triesToFound + triesSumand;
            end
            
            % In small neighborgs or when tabu has lots of elements it can
            % happen that we don´t find a plausible solution
            if (~foundPlausibleSol)
                continue;
            end
%-------------------------------------------
% GENERATE E FOR EACH SUB SET OF SOLUTIONS
%-------------------------------------------           
            
            for VExtrIndex=1:problemObj.lumCount
                E = zeros(n,1);
                
                lowerLimit = (((VExtrIndex - 1) * problemObj.nvars ) + 1);
                VExtrAux = VExtr(lowerLimit:lowerLimit + (problemObj.nvars - 1));

                % Calculate the emision merging both the view hemicube with the
                % radiance hemicube taken from the luminaire
                position = find(problemObj.selectedIndexesOrdered==positionObj.positionsMapMatrix(VExtrAux(1), VExtrAux(2)));

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

                % Acum E for the h candidate
                E = E + lumE(1:problemObj.triangleCount);
                
                positionsAux(VExtrIndex, candidateIndex) = position;
                selectedLuminary(VExtrIndex, candidateIndex) = VExtrAux(3);
                z1Angles(VExtrIndex, candidateIndex) = z1;
                y2Angles(VExtrIndex, candidateIndex) = y2;
                z3Angles(VExtrIndex, candidateIndex) = z3;
            end
            
            % Save emission into matrix
            E = E(1:problemObj.triangleCount);
            E = (E.*lowRankObj.Rmean)./lowRankObj.areaPatches;
            Edyn(:, candidateIndex) = E;
            % Save variables used for the h candidate
            variablesForCandidate(:, candidateIndex) = VExtr;
            candidateIndex = candidateIndex + 1;
        end
        
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        %<<<<<<<<<<< MIN/MAX OBJECTIVE <<<<<<<<<<<
        %<<<<<<< Do matrix-matrix product <<<<<<<<
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        hMaxCalculated = size(Edyn,2);
        % In case we cannot have a single candidate for this neighborg
        % search in the following
        if hMaxCalculated < 1
            [iii jjj nextjjj] = nextMatrixPoint(iii, jjj, nextjjj, problemObj.maxdeltavalues, problemObj.maxlmstepvalues);
            delta = problemObj.deltavalues(iii);
            lmstep = problemObj.lmstepvalues(jjj);
            continue;
        end
        
        
        NroIterIf = NroIterIf + hMaxCalculated;
        objetives = zeros(hMaxCalculated, 1);

        %---------------------------------------
        % Find Uniformity/Dispersion:
        %---------------------------------------
%         B = -lowRankObj.Yp * (lowRankObj.V * Edyn) + Edyn(:, problemObj.objectiveSelectedIndexesOrdered);
%         BB = zeros(triangleCount, hMax);
%         BB(problemObj.objectiveSelectedIndexesOrdered, :) = B; 
%         
%         for h = 1:hMaxCalculated
%             objectives(h) = std(B(:, h))/mean(B(:, h));
%         end
        %---------------------------------------

        %---------------------------------------
        % Optimize Power Find min and max luxes:
        %---------------------------------------
%         B = -lowRankObj.Yp * (lowRankObj.V * Edyn) + Edyn(:, problemObj.objectiveSelectedIndexesOrdered);
%         BB = zeros(problemObj.triangleCount, hMax);
%         BB(problemObj.objectiveSelectedIndexesOrdered, :) = B;
%         
%         for h = 1:hMaxCalculated
%             objetives(h) = sum(luminarieObj.luminairePower(selectedLuminary));
%             % Penalty
%             c1 = ((max(50 - (min(B(:, h))/areaPatch), 0)^2)*100);
%             if (c1 > 0)
%                 c1 = c1 + luminarieObj.maxPower * problemObj.lumCount; % c1 >= luminarieObj.maxPower
%             end
%             c2 = ((max((max(B(:, h))/areaPatch) - 750, 0)^2)*100);
%             if (c2 > 0)
%                 c2 = c2 + luminarieObj.maxPower * problemObj.lumCount; % c2 >= luminarieObj.maxPower
%             end
% 
%             objectives(h) = objectives(h) + c1 + c2;
%         end
        %---------------------------------------

        %---------------------------------------
        % Find exact solution:
        %---------------------------------------
        B = -lowRankObj.Y * (lowRankObj.V * Edyn) + Edyn;
        
        for h = 1:hMaxCalculated
            diffVect = (B(:, h) - problemObj.target);
            objectives(h) = sqrt(sum(diffVect.*diffVect));
        end
        %---------------------------------------

        %---------------------------------------
        % _Maximize_ sum of luxes on a surface:
        % Restrictions: max luminaries power
        %---------------------------------------
%         B = -lowRankObj.Yp * (lowRankObj.V * Edyn) + Edyn(:, problemObj.objectiveSelectedIndexesOrdered);
%         
%         for h = 1:hMaxCalculated
%             objetives(h) = sum(B(:, h)/areaPatch);
%             % Penalty
%             c1 = -1*((max(restrictionObj.maxPower - sum(luminarieObj.luminairePower(selectedLuminary(:, h))), 0)^2)*100);
%             if (c1 < 0)
%                 c1 = c1 + luminarieObj.maxPower * problemObj.lumCount; % c1 >= luminarieObj.maxPower
%             end
%             objectives(h) = -1*(c1 + objectives(h));
%         end
        %---------------------------------------

        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        % Evaluate each candidate from the neighborhood
        for h = 1:hMaxCalculated
            
            objetivo = objectives(h);
            
            if (objetivo <= MENOR)
                VExtrOpt = variablesForCandidate(:, h);              
                if objetivo < MENOR
                    ProbarVariasVeces = 1;
                    halleMENOR = 1;
                    Efound = Edyn(:, h);
                    MENOR = objetivo;
                    Nsolu = Nsolu + 1;
                    EEE = [Efound Efound Efound];
                    emision = EEE;
                    orderedPos = variablesForCandidate(3, h);
                    [(NroIterIf - hMaxCalculated + h), objetivo, selectedLuminary(:, h)', positionsAux(:, h)', z1Angles(:, h)', y2Angles(:, h)', z3Angles(:, h)']

                    bestPositions = positionsAux(:, h)';
                    bestLuminaires = selectedLuminary(:, h)';
                    bestIt = NroIterIf - hMaxCalculated + h;
                    bestAngles = [z1Angles(:, h)', y2Angles(:, h)', z3Angles(:, h)'];

                    for VExtrIndex=1:problemObj.lumCount
                        aux = luminarieObj.orderedLuminaireNameMap(selectedLuminary(VExtrIndex, h));
                        aux{1}
                        positionsNameMap(positionsAux(VExtrIndex, h))
                    end

                    if (objetivo == 0 || objetivo < 0.0001)
                        tabuFilled = 1;
                    end
                end
            end
        end

        
        % Reset neighborhood
        if (ProbarVariasVeces) && ~(halleMENOR)
            ProbarVariasVeces = 0;
            iii = 1;jjj = 1;nextjjj = 1;
            lmstep = problemObj.lmstepvalues(jjj);
            delta = problemObj.deltavalues(iii);
            continue;
        end
        
        % Update neighborhood
        if ~(ProbarVariasVeces)
            %En diagonal, siempre intenta ir abajo a la izquierda si no
            %puede vuelve a iniciar en el primer elemento disponible de
            %la primer columna que no fue consultada (o la ultima
            %columna de la matriz en caso que no hayan mas)
            [iii jjj nextjjj] = nextMatrixPoint(iii, jjj, nextjjj, problemObj.maxdeltavalues, problemObj.maxlmstepvalues);
            delta = problemObj.deltavalues(iii);
            lmstep = problemObj.lmstepvalues(jjj);
            %[delta, lmstep,h]
        end
    end
    singleTimer = toc(singleTimerId);
    %%%%%%%%%%%%%Hallo la solución real%%%%%%%%%%%%%%%%%%%%%%%
    VE=zeros(k,1);
    for i=1:n, VE(VI(i))=VE(VI(i))+Efound(i); end
    Creal=-Y*(VE);
    VE=zeros(k,1);
    %%%%%%%%%%%%%Hallé la solución real%%%%%%%%%%%%%%%%%%%%%%%
    
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

