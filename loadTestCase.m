try
    sleepTime = 3600;
    here = pwd;
    cantRuns = 30;
    % Load queue of tests to run
    ended = false;
    while (~ended)
        load('testCases.mat', 'testCases');
        if length(testCases)
            % Get test that needs to be executed
            test = testCases(1, :);

            % Update pool of tests removing the one that is being executed
            testCases = testCases(2:end, :);
            save('testCases.mat', 'testCases');

            % Set parameters for tests
            lumCount = test(1)
            cantRepeat = test(2)
            MaxNroIterIf = test(3)
            try 
                Algorithm = test(4)
            catch
                Algorithm = 0;
                disp('No algorithm found! defaulting to 0')
            end
            
            
            ended = true;
        else
            disp(['No test found! sleeping:' int2str(sleepTime)])
            pause(sleepTime);
        end    
    end
    disp('Finalizing....')
catch
  disp('Sorry some unhandled excpetion occurred!!!... closing matlab...');
  exit;
end