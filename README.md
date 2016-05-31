# voigtfit.m

Fits data to a Voigt profile. Voigt models are commonly used to fit XPS 
spectra.

[estimates, model] = voigtfit(x, y, initGuess, peakBounds) fits the x and y 
data to one or more Voigt profile models, returning the Voigt profile 
parameters in the vector estimates and the fit function handle, model.
The Voigt model fit is initialized with the parameters in initGuess. The 
background is fit to a 3rd order polynomial which is fit to all data except 
for a user-defined region contained within the upper and lower bounds given 
in peakBounds (see below).

[estimates, model] = voigtfit(x, y, initGuess, peakBounds, bkgdFitOrder) fits 
the background to a polynomial of order bkgdFitOrder.

Inputs:

  xdata - A column vector containing x data.

  ydata - A column vector containing y data.

  initGuess - A row vector containing initial guesses for the nonlinear fit 
              parameters. There are 3 parameters for each peak: 
              peak center value, gamma, and sigma. For multiple peaks, 
              use the format 
              [peak_1, gamma_1, sigma_1, peak_2, gamma_2, sigma_2, ...].
              For more info on Voigt profile parameters, see
              https://en.wikipedia.org/wiki/Voigt_profile.

  peakBounds - A 1x2 vector of the form [LB, UB]. The user must select upper 
              and lower bounds for the region containing the peak(s).
              The background is fit to a polynomial by excluding this data.

  bkgdOrder - Optional input. The order of the polynomial used for the 
              background fit. The default value is 2.

Outputs:

  estimates - A row vector containing the final estimates for nonlinear 
              parameters. The elements are organized identically to the 
              initGuess input parameter.

  model - A handle to the model function.
  

I used Steven G. Johnsons's Faddeeva package, which needs to be in the 
search path. It can be found at 
http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package.
