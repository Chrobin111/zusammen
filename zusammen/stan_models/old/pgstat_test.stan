functions {
#include pgstat.stan
}

data {
    int N;
    vector[N] observed_counts;
    vector[N] background_counts;
    vector[N] background_errors;
    vector[N] expected_model_counts;
    array[N] int idx_background_zero;
    array[N] int idx_background_nonzero;
    int N_bkg_zero;
    int N_bkg_nonzero;
}

transformed data {
    print(
        pgstat(
            observed_counts, background_counts, background_errors, expected_model_counts, idx_background_zero[ :N_bkg_zero], idx_background_nonzero[ :N_bkg_nonzero]
        )
    );
    print(num_elements(idx_background_zero));
    print(num_elements(idx_background_nonzero));
}