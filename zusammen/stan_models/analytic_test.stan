functions {
#include cpl.stan
}

transformed data {
    real x_r[0];
    int x_i[0];

    array[3] real alpha = { -0.8, -1, -1.2 };
    array[3] real ec = { 1e1, 1e2, 1e3 };

    for (i in 1:3){
        for (j in 1:3){
                print("alpha = ", alpha[i]);
                print("ec = ", ec[j]);

                print("Exakt: ", ggrb_int_cpl(alpha[i], ec[j], 10., 1.e4));

                array[3] real theta = {1., alpha[i], ec[j]};
                print("Numerisch: ", integrate_1d(cpl_flux_integrand, 10., 1.e4, theta, x_r, x_i));

                print("");
        }
    }
}