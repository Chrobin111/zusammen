functions {
#include cpl.stan

    real differential_flux_scalar (real energy, real xc, real[] theta, real[] x_r, int[] x_i) {
        real K = theta[1];
        real ec = theta[2];
        real alpha = theta[3];
        return K * pow(energy/100., alpha) * exp(-energy/ec);
    }
}

data{
  int<lower=1> N_intervals; // number of intervals
  int max_n_echan; // number of maximum energy side channels
  int max_n_chan; // maximum number of channels

  array[N_intervals] int<lower=0> N_dets; // number of detectors per data type
  array[N_intervals, max(N_dets)] int<lower=0> N_chan; // number of channels in each detector
  array[N_intervals,  max(N_dets)] int<lower=0> N_echan; // number of energy side channels in each detector

  array[N_intervals, max(N_dets)] real exposure; // exposure time

  array[N_intervals, max(N_dets)] vector[max_n_echan] ebounds_hi;
  array[N_intervals, max(N_dets)] vector[max_n_echan] ebounds_lo;

  array[N_intervals, max(N_dets)] matrix[max_n_chan, max_n_echan] response;

  vector[N_intervals] K;
  vector[N_intervals] ec;
  vector[N_intervals] alpha;
}


transformed data{
    array[N_intervals, max(N_dets)] vector[max_n_echan] ebounds_add;
    array[N_intervals, max(N_dets)] vector[max_n_echan] ebounds_half;

    for (n in 1:N_intervals) {
        for (m in 1:N_dets[n]) {
            ebounds_half[n, m, :N_echan[n, m]] = 0.5*(ebounds_hi[n, m, :N_echan[n, m]]+ebounds_lo[n, m, :N_echan[n, m]]);
            ebounds_add[n, m, :N_echan[n, m]] = (ebounds_hi[n, m, :N_echan[n, m]] - ebounds_lo[n, m, :N_echan[n, m]])/6.0;
        }
    }

    array[max(N_dets)] vector[max_n_chan] expected_model_counts;

    for (n in 1:N_intervals) {
        for (m in 1:N_dets[n]) {
            expected_model_counts[m, : N_chan[n,m]] = (
                response[n, m,:N_chan[n,m], :N_echan[n,m]] * integral_flux(
                    ebounds_lo[n, m, :N_echan[n, m]],
                    ebounds_hi[n, m, :N_echan[n, m]],
                    ebounds_add[n, m, :N_echan[n, m]],
                    ebounds_half[n, m, :N_echan[n, m]],
                    K[n],
                    ec[n],
                    alpha[n]
                )
            ) * exposure[n,m];
        }
    }

    print(expected_model_counts);


    real x_r[0];
    int x_i[0];
    array[max(N_dets)] vector[max_n_chan] expected_model_counts_int;

    for (n in 1:N_intervals) {
        for (m in 1:N_dets[n]) {
            array[3] real theta = {K[n], ec[n], alpha[n]};

            vector[N_echan[n,m]] integrated_flux;

            for (l in 1:N_echan[n,m]){
                integrated_flux[l] = integrate_1d(
                    differential_flux_scalar,
                    ebounds_lo[n, m, l],
                    ebounds_hi[n, m, l],
                    theta,
                    x_r,
                    x_i
                );
            }
            expected_model_counts[m, : N_chan[n,m]] = (
                response[n, m,:N_chan[n,m], :N_echan[n,m]] * integrated_flux
            ) * exposure[n,m];
        }
    }

    print(expected_model_counts_int);

}