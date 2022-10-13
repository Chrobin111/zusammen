real x(real z, real Om) {
    return (1 - Om) / (Om * (1 + z) * (1 + z) * (1 + z));
}


real my_Phi(real my_x) {
    return (1 + 1.32 * my_x + 0.4415 * my_x * my_x + 0.02656 * my_x * my_x * my_x) / ( 1 + 1.392 * my_x + 0.5121 * my_x * my_x + 0.03944 * my_x * my_x * my_x);
}


real dl_func(real z) {
    real H0 = 2.2e-18; // 1/s
    real Om = 0.27;
    real c = 299792458; // m/s

    return (2 * c / H0) * ((1 + z) / sqrt(Om)) * (my_Phi(x(0, Om)) - my_Phi(x(z, Om)) / sqrt(1 + z)) * 100;
}