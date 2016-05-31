% voigtfit_demo.m
% Simulate Voigt profile fit for 2 sets of simulated data. One has 1 peak, the 
% other has 2.

x = linspace(25,125,1000);
params1 = [65, 0.1, 1];
params2 = [75, 0.2, 0.9];
data1 = myvoigt(x, params1(1), params1(2), params1(3));
data2 = myvoigt(x, params2(1), params2(2), params2(3));
% Simulate Gaussian white noise plus a background emission spectrum.
noise = 0.05*randn(size(x)) + 0.1*polyval([0.02,0.01], x);
y1 = data1 + noise; % Signals y1 and y2 are Voigt peaks plus noise, background.
y2 = y1 + data2;

% Fit signal y1.
initGuess1 = [64, 0.5, 0.1];
[estimates1, model1] = voigtfit(x, y1, initGuess1, [60, 70]);
disp('Single peak fit results [peak1, gamma1, sigma1]');
disp(estimates1);
pause(2); % Display single peak fit and pause.

% Fit signal y2.
initGuess2 = [62, 0.5, 0.5, 79, 0.5, 0.5];
[estimates2, model2] = voigtfit(x, y2, initGuess2, [60, 80]);
disp('Double peak fit results [peak1, gamma1, sigma1, peak2, gamma2, sigma2]');
disp(estimates2);