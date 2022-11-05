functions {
#include cpl.stan
}

data {
    int N;
    array[N] real alpha;
    array[N] real ec;
    array[N] real energy_flux;
}

transformed data {
    real erg2kev = 6.24151e8;
    real x_r[0];
    int x_i[0];
    real K;
    vector[3] ebounds_lo = [100.,105.,110.]';
    vector[3] ebounds_hi = [105.,110.,115.]';
    vector[3] ebounds_add = (ebounds_hi - ebounds_lo) / 6.;
    vector[3] ebounds_half = (ebounds_hi + ebounds_lo) / 2.;

    for (i in 1:N){
        for (j in 1:N){
            for (k in 1:N) {
                print("alpha = ", alpha[i]);
                print("ec = ", ec[j]);
                print("flux = ", energy_flux[k]);

                K = erg2kev * energy_flux[k] * inv(ggrb_int_cpl(alpha[i], ec[j], 10, 10e4));
                print("Exakt: ", K);

                print("Integral Flux: ", integral_flux(ebounds_lo, ebounds_hi, ebounds_add, ebounds_half, K, ec[j], alpha[i]));

                // array[3] real theta = {1., alpha[i], ec[j]};
                // print("Numerisch: ", erg2kev * energy_flux[k] * inv(integrate_1d(cpl_flux_integrand, 10., 1.e4, theta, x_r, x_i)));

                print("");
            }
        }
    }
}