function vp = myvoigt(xx, peakVal, gam, sig)
% MYVOIGT returns Voigt profile data.
% xx            Vector of x values.
% peakVal       Peak center value.
% gammaVal      Gamma value.
% sig           Sigma value.
z = arrayfun(@(q) ((q-peakVal)+i*gam)/(sqrt(2)*sig), xx);
vp = (1/(sig*sqrt(2*pi))) * real(Faddeeva_w(z)); % Get Voigt from Faddeeva fn.
vp = vp./max(vp);
end