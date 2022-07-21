data{

    int<lower=0> N; // number of GRBs
    int<lower=0> Nunknown; // number of GRBs with unknown redshift
    int<lower=0> Ngroups; // number of groups TODO

    vector[N] z_known; // known redshifts
    vector[N] Nz;
    int group[N];
    int unknown_group[N];
    int known[N]; // array of booleans whether redshift is known
    real maxSlope; // maximum value for gamma

    vector[N] FE_obs; // observed flux
    vector[N] Ep_obs; // observed peak energy
    vector[N] Ep_sig; // sigma of peak energy
    vector[N] FE_sig; // sigma of flux

}



parameters {

    // latent variables
    vector<lower=-2,upper=5>[N] Ep_true; // true peak energy
    vector[N] FE_true; // true peak energy

    vector<lower=50>[Ngroups] Nrest; // GC normalization
    vector<lower=0>[Ngroups] gamma; // exponent
    vector<lower=0.0,upper=15>[Nunknown] z; // unknown redshifts
    vector<lower=0>[Ngroups] int_scatter_sq; // intrinsic scatter squared

    // hyperpriors
    real<lower=0> gamma_mu_meta;
    real Nrest_mu_meta;
    real<lower=0> gamma_sig_meta;
    real<lower=0> Nrest_sig_meta;

}



transformed parameters {

    vector<lower=0>[Ngroups] int_scatter; // intrinsic scatter
    for(i in 1:Ngroups){
        int_scatter[i] = sqrt(int_scatter_sq[i]); // before: <-
    }

}



model{

    gamma_sig_meta ~ cauchy(0., 2.5);
    Nrest_sig_meta ~ cauchy(0., 2.5);
    gamma_mu_meta ~ normal(0, maxSlope);
    Nrest_mu_meta ~ normal(52, 5);
    gamma ~ normal(gamma_mu_meta, gamma_sig_meta);
    Nrest ~ normal(Nrest_mu_meta, Nrest_sig_meta);
    int_scatter_sq ~ cauchy(0, 2.5)
    z ~ uniform(0, 15);
    Ep_true ~ uniform(-2, 5);
    Ep_obs ~ normal(Ep_true, Ep_sig);

    for(i in 1:N){

        if (known[i] = 0){

            FE_true[i] ~ normal(
                Nrest[group[i]] - ( 1.099 + 2 * log10(DL(z[unknown_group[i]])) ) // DL = luminosity distance?
                    + gamma[group[i]] * ( log10(1+z[unknown_group[i]]) + Ep_true[i]-2 ),
                int_scatter[group[i]]
            );

        } else {

            FE_true[i] ~ normal(
                Nrest[group[i]] - Nz[i] + gamma[group[i]] * ( log10(1+z_known[i]) + Ep_true[i] - 2 ),
                int_scatter[group[i]]
            );

        }

    }

    FE_obs ~ normal(FE_true,FE_sig);

}