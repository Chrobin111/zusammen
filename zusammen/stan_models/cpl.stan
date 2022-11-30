// CPL normalization using Gamma functions
real ggrb_int_cpl(real alpha, real ec, real emin, real emax) {

  real i1 = gamma_q(2 + alpha, emin / ec) * tgamma(2 + alpha);
  real i2 = gamma_q(2 + alpha, emax / ec) * tgamma(2 + alpha);

  return -square(ec) * (i2 - i1);

}



// alternative Cut-off Power Law
vector cpl_alt(vector energy, real norm, real ec, real alpha) {

  // real piv = 100.;

  return norm * pow(energy / ec, alpha) .* exp(-energy/ec);

}



// energy-differentiated flux
vector differential_flux_alt( vector energy, real norm, real ec, real alpha) {

  return cpl_alt(energy, norm, ec, alpha);

}



// flux
vector integral_flux_alt(vector ebounds_lo, vector ebounds_hi, vector ebounds_add, vector ebounds_half, real norm, real ec, real alpha) {

  return (ebounds_add
	  .* (differential_flux_alt(ebounds_lo, norm, ec, alpha)
	      + 4 * differential_flux_alt(ebounds_half, norm, ec, alpha)
	      + differential_flux_alt(ebounds_hi, norm, ec, alpha)));

}



// Cut-off Power Law
vector cpl(vector energy, real norm, real ec, real alpha) {

  // real piv = 100.;

  return norm * pow(energy, alpha) .* exp(-energy/ec);

}



real cpl_indi(real energy, real K, real alpha, real ec) {

  // real piv = 100.;

  return K * pow(energy, alpha) * exp(-energy/ec);

}



real cpl_flux_integrand(real x, real xc, real[] theta, real[] x_r, int[] x_i) {

  real out = x * cpl_indi(x, theta[1], theta[2], theta[3]);

  return out;

}



// energy-differentiated flux
vector differential_flux( vector energy, real norm, real ec, real alpha) {

  return cpl(energy, norm, ec, alpha);

}



// flux
vector integral_flux(vector ebounds_lo, vector ebounds_hi, vector ebounds_add, vector ebounds_half, real norm, real ec, real alpha) {

  return (ebounds_add
	  .* (differential_flux(ebounds_lo, norm, ec, alpha)
	      + 4 * differential_flux(ebounds_half, norm, ec, alpha)
	      + differential_flux(ebounds_hi, norm, ec, alpha)));

}