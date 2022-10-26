data {
    int N;
    vector[N] alpha;
    vector[N] log_ec;
    vector[N] log_Nrest;
    vector[N] dl;
    vector[N] gamma;
    vector[N] z;
}

transformed data {
    vector[N] log_epeak;
    vector[N] epeak;
    vector[N] log_energy_flux;
    vector[N] energy_flux;

    log_epeak = log10(2+alpha) + log_ec;

    log_energy_flux = log_Nrest - (1.099 + 2 * log10(dl)) + gamma .* (log10(1 + z) + log_epeak - 2);
    energy_flux = pow(10, log_energy_flux);

    print(log_epeak, log_energy_flux);
    print(epeak, energy_flux);
}