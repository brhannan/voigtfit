function [estimates, model] = voigtfit(xdata, ydata, initGuess, peakBounds, varargin)
% VOIGTFIT Fit data to a Voigt profile model.
%
% [estimates, model] = voigtfit(x, y, initGuess, peakBounds) fits the x and y
% data to one or more Voigt profile models, returning the Voigt profile
% parameters in the vector estimates and the fit function handle, model.
% The Voigt model fit is initialized with the parameters in initGuess. The
% background is fit to a 3rd order polynomial by excluding the data contained
% within the upper and lower bounds specified by peakBounds.
%
% [estimates, model] = voigtfit(x, y, initGuess, peakBounds, bkgdFitOrder) fits
% the background to a polynomial of order bkgdFitOrder.
%
% Inputs:
%   xdata       A vector of x data.
%   ydata       A vector of y data.
%   initGuess   A row vector containing initial guesses for the nonlinear fit
%               parameters. There are 3 parameters for each peak:
%               peak center value, gamma, and sigma. For multiple peaks,
%               use the format
%               [peak_1, gamma_1, sigma_1, peak_2, gamma_2, sigma_2, ...].
%               For more info on Voigt profile parameters, see
%               https://en.wikipedia.org/wiki/Voigt_profile.
%   peakBounds  A 1x2 vector of the form [LB, UB]. The user must select upper
%               and lower bounds for the region containing the peak(s).
%               The background is fit to a polynomial by excluding this data.
%   bkgdOrder   Optional input. The order of the polynomial used for the
%               background fit. The default value is 2.
% Outputs:
%   estimates   A row vector containing the final estimates for nonlinear
%               parameters. The elements are organized identically to the
%               initGuess input parameter.
%   model       A handle to the model function.

% Resources:
% MathWorks curve fitting tutorial
% mathworks.com/help/matlab/math/example-curve-fitting-via-optimization.html
% Steven G. Johnsons's Faddeeva package
% http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package

% Set optional input value.
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:voigtfit:TooManyInputs', ...
        'accepts at most 1 optional input');
end
optargs = {2}; % Set default for optional input.
optargs(1:numvarargs) = varargin;
[bkgdOrder] = optargs{:};

X_TOL_STOP = 1e-3; % Termination tolerance for x.
MAX_EVALS = 1e8; % optimset max number of function evaluations.

if size(xdata,1) < size(xdata,2)
    xdata = xdata';
end
if size(ydata,1) < size(ydata,2)
    ydata = ydata';
end

assert(~rem(numel(initGuess), 3), ...
    'The number of elements in initGuess must be a multiple of 3.');
assert(numel(peakBounds) == 2, ...
    'The number of elements in peakBounds must equal 2.');
assert(peakBounds(1) < peakBounds(2), ...
    'peakBounds must have the form [lowerBound, upperBound].');
assert(size(xdata,2)==1 && size(ydata,2)==1, ...
    'X and Y inputs must be vectors.');
assert(~rem(numel(initGuess), 3), ...
    ['Number of parameters passed to model must be equal to the product of ' ...
    'the number of peaks and the number of nonlinear parameters to fit ' ...
    'to each peak (3).']);

% Fit the background to a polynomial.
isPeak = xdata < peakBounds(2) & xdata > peakBounds(1);
bkgdPoly = polyfit(xdata(~isPeak), ydata(~isPeak), bkgdOrder);
bkgdFit = polyval(bkgdPoly, xdata);

ydata = ydata - bkgdFit;
peaksCell = calc_model_from_nlparams(initGuess, xdata, ydata);
% Plot data and model for first guess parameters.
mycmap = lines(numel(initGuess)/3);
figure(1);
    clf;
    hData = plot(xdata, ydata + bkgdFit, 'k', 'DisplayName', 'Data');
    hold on;
    hBkgd = plot(xdata, bkgdFit, 'm--', 'DisplayName', 'Background');
    hold off;
    hold on;
    hFit = plot(xdata, peaksCell{1} + bkgdFit, 'r:', 'LineWidth', 2, ...
                'DisplayName', 'Fit');
    hold off;
    xlim([min(xdata), max(xdata)]);
    title('Inital model data');
    set(gca, 'XDir', 'reverse'); % Reverse x (XPS plot convention).
    set(gcf, 'color', 'w');
    legend([hData, hBkgd, hFit]); %, peakHandlesArray]);
    grid on;
    pause(1.5);

% Fit data to model with fminsearch. Plot model state at each iteration.
model = @fit_voigt;
outputFcn = @(x, optimvalues, state) fitoutputfun(x, optimvalues, state,...
    xdata, ydata, bkgdFit, hFit);
options = optimset('OutputFcn', outputFcn, 'TolX', X_TOL_STOP, ...
        'MaxFunEvals', MAX_EVALS);
estimates = fminsearch(model, initGuess, options);
    function [sse, fittedCurve] = fit_voigt(nlParamVec)
        A = zeros(numel(xdata), numel(nlParamVec)/3);
        for k = 1:numel(nlParamVec)/3
            A(:,k) = voigt(xdata, nlParamVec(3*k-2), nlParamVec(3*k-1), ...
                 nlParamVec(3*k))';
        end
        c = abs(A\ydata);
        fittedCurve = A*c;
        errorVec = sum(fittedCurve - ydata);
        sse = sum(errorVec.^2);
    end
end % main


function vp = voigt(xx, peakVal, gam, sig)
% VOIGT returns Voigt profile data.
% xx            Vector of x values.
% peakVal       Peak center value.
% gammaVal      Gamma value.
% sig           Sigma value.
z = arrayfun(@(q) ((q-peakVal)+i*gam)/(sqrt(2)*sig), xx);
vp = (1/(sig*sqrt(2*pi))) * real(Faddeeva_w(z)); % Get Voigt from Faddeeva fn.
vp = vp./max(vp);
end


function stop = fitoutputfun(nlpVec, optimvalues, state, x, y, bkgd, ...
                                hTotFit)
% FITOUTPUT Output function used by FITDEMO
% Copyright 1984-2004 The MathWorks, Inc.
% Modified by BH to display additional curves + background fit.
stop = false;
peaksCell = calc_model_from_nlparams(nlpVec, x, y);
switch state
    case 'init'
        set(hTotFit, 'ydata', peaksCell{1} + bkgd);
        drawnow;
        title('Voigt profile fit results');
    case 'iter'
        set(hTotFit, 'ydata', peaksCell{1} + bkgd);
        drawnow;
    case 'done'
        hold off;
end % switch
end


function peaksCell = calc_model_from_nlparams(nlpVec, x, y)
% Outputs fit results in the cell array peaksCell. The 1st element is the
% total fit. When >1 peaks are fit to the data, the n+1th element of peaksCell
% corresponds to the n+1th peak.
peaksCell = cell(1, numel(nlpVec)/3);
A = zeros(numel(x), numel(nlpVec)/3);
for k = 1:numel(nlpVec)/3
    A(:,k) = voigt(x, nlpVec(3*k-2), nlpVec(3*k-1), nlpVec(3*k));
end
c = A\y;
peaksCell{1} = A*c; % 1st element of output is the total fit.
% Get z array for each individual peak if >1 peaks are fit to data.
if numel(nlpVec) > 3
    for nPeak = 1:numel(nlpVec)/3
        cNow = A(:,nPeak)\y;
        peaksCell{nPeak+1} = A(:,nPeak)*cNow; % Elements >1 are single peaks.
    end
end
end
