function t = register(moving,fixed, transformType)
  % Rigidly register two monomodal images in imregconfig sense
  %
  % t = function register(I)
  %   I, J - n x m dimensional images to be registered
  %   t - output from  imregtform from I to J

  [opt,met] = imregconfig('monomodal');
  %opt.MaximumStepLength = 10^(-4);
  %opt.MaximumIterations = 100;
  %factor = 100;
  %opt.MinimumStepLength = opt.MinimumStepLength/factor;
  %opt.MaximumStepLength = opt.MaximumStepLength/factor;
  %opt.MaximumIterations = opt.MaximumIterations*factor;
  %levels = round(log2(max(size(I)))/3);
  opt.MaximumIterations = 1000;
  levels = 1;

  %t = imregtform(I, J,'rigid',opt,met,'PyramidLevels',levels,'InitialTransform',t,'DisplayOptimization',false);
  t = imregtform(moving, fixed, transformType, opt, met, 'PyramidLevels',levels,'DisplayOptimization',true);