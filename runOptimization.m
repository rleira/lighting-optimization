try
    loadTestCase;
    
    loadData;

    if (~Algorithm)
        disp('Running optimizeLuminariesN ...')
        optimizeLuminariesN;
    end
    
    if (Algorithm == 1)
        disp('Running optimizeLuminariesNExpV ...')
        optimizeLuminariesNExpV;
    end
    
catch
  disp('Sorry some unhandled excpetion occurred!!!... closing matlab...');
  exit;
end