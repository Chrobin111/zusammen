# Files
## Stan models
### alpha_correlation
- data:
    - $N_\mathrm{intervals}$
    - $N_\mathrm{chan,max}$
    - $N_\mathrm{echan,max}$
    - $N_\mathrm{det}$
    - $N_\mathrm{chan}$
    - $N_\mathrm{echan}$
    - grb_id
    - $N_\mathrm{GRB}$
    - ebounds_hi
    - ebounds_lo
    - $N_\mathrm{obs}$
    - $N_\mathrm{bg}$
    - BG errors
    - idx_background_zero
    - idx_background_nonzero
    - $N_\mathrm{BG,zero}$
    - $N_\mathrm{BG,nonzero}$
    - exposure
    - response
    - mask
    - $N_\text{channels used}$
    - $dl$
    - $z$
    - $N_\text{gen spectra}$
    - $E_\mathrm{model}$
    - $N_\mathrm{correlation}$
    - model_correlation
- transformed data:
    - ebounds_half
    - ebounds_add
    - $N_\text{total channels}$
- parameters: 
    - $\alpha$
    - $\beta$
    - $\log E_\mathrm{peak}$
    - $\gamma$
    - $\delta$
- transformed parameters:
    - $\gamma$
    - $\delta$
    - $E_\mathrm{peak}$,
    - $\log L$
    - $L$
    - expected model counts
- model:
    - $\alpha \sim \mathcal{N}(-1,0.5)$
    - $\beta \sim \mathcal{N}(-3,1)$
    - $\log E_\mathrm{peak} \sim \mathcal{N}(2,1)$
    - $\gamma_\mu \sim \mathcal{N}(3,1)$
    - $\delta_\mu \sim \mathcal{N}(-4,4)$
    - $\gamma_\mathrm{off} \sim \mathcal{N}(0,1)$
    - $\delta_\mathrm{off} \sim \mathcal{N}(0,1)$
    - $\gamma_\sigma \sim \mathcal{N}(0,5)$
    - $\delta_\sigma \sim \mathcal{N}(0,5)$
    - log_like / target: ???
- generated quantities: ???

### band_grb
- ggrb_into_pl:
- ggrb_int_cpl:
- band_precalculation:
- differential_flux:
- integral_flux: 

### big_fuck_correlation
- data: usual
- transformed data:
    - ebounds_half
    - ebounds_add
    - $N_\text{total channels}$
    - log_zp1_grb
    - log_dl2_grb
- parameters:
    - $\alpha$
    - $\beta$
    - $\log E_\mathrm{peak}$
    - $\gamma$
    - $\delta$
    - $\phi$
- transformed parameters:
    - $\gamma$
    - $\delta$
    - $\phi$
    - $E_\mathrm{peak}$
    - $\log L$
    - $L$
    - expected model counts
- model:
    - $\alpha \sim \mathcal{N}(-1,0.5)$
    - $\beta \sim \mathcal{N}(-3,1)$
    - $\log E_\mathrm{peak} \sim \mathcal{N}(2,1)$
    - $\gamma_\mu \sim \mathcal{N}(0,5)$
    - $\phi_\mu \sim \mathcal{N}(0,5)$
    - $\delta_\mu \sim \mathcal{N}(0,4)$
    - $\gamma_\mathrm{off} \sim \mathcal{N}(0,1)$
    - $\phi_\mathrm{off} \sim \mathcal{N}(0,1)$
    - $\delta_\mathrm{off} \sim \mathcal{N}(0,1)$
    - $\gamma_\sigma \sim \mathcal{N}(0,5)$
    - $\phi_\sigma \sim \mathcal{N}(0,5)$
    - $\delta_\sigma \sim \mathcal{N}(0,5)$
    - log_like / target: ???
- generated quantities: ???

### correlation
- data: usual
- transformed data:
    - ebounds_half
    - ebounds_add
    - $N_\text{total channels}$
- parameters:
    - $\alpha$
    - $\beta$
    - $E_\mathrm{peak}$
    - $\gamma$
    - $\delta$
- transformed parameters:
    - $\gamma$
    - $\delta$
    - $\log L$
    - expected model counts
- model:
    - $\alpha \sim \mathcal{N}(-1,0.5)$
    - $\beta \sim \mathcal{N}(-3,1)$
    - $E_\mathrm{peak} \sim \mathcal{N}(500,500)$
    - $\gamma_\mu \sim \mathcal{N}(3,1)$
    - $\delta_\mu \sim \mathcal{N}(-4,4)$
    - $\gamma_\mathrm{off} \sim \mathcal{N}(0,1)$
    - $\delta_\mathrm{off} \sim \mathcal{N}(0,1)$
    - $\gamma_\sigma \sim \mathcal{N}(0,5)$
    - $\delta_\sigma \sim \mathcal{N}(0,5)$
    - log_like / target: ???
- generated quantities: ???

### cpl_interval_fold:
- partial_log_like: 

### cpl_simple_chunked:
- data: usual
- transformed data:
    - ebounds_half
    - ebounds_add
    - $N_\text{total channels}$
    - keV to erg
    - erg to keV
    - $x_r$
    - $x_i$
- parameters:
    - $\alpha$
    - $\log E_C$
    - $\log L$
- transformed parameters:
    - $E_C$
    - $L$
    - $K = L \cdot \left( \int_{10}^{1000} dx\ x \cdot K \cdot (E/100)^\alpha \cdot \mathrm{e}^{-E/100} \right)^{-1}$
- model:
    - $\log L \sim \mathcal{N}(0,1)$
    - $\alpha \sim \mathcal{N}(-1,0.5)$
    - $\log E_C \sim \mathcal{N}(2,1)$
    - target: ???
- generated quantities: ???

### cpl
- ggrb_int_cpl:
- cpl: $N \cdot (E/100)^\alpha \cdot \mathrm{e}^{-E/E_C}$
- cpl_indi: $K \cdot (E/100)^\alpha \cdot \mathrm{e}^{-E/E_C}$
- cpl_flux_integrand: $x \cdot \mathrm{cpl_indi}(x,\theta)$
- differential_flux: $L_\mathrm{diff}(E,N,E_C,\alpha) = \mathrm{cpl}(E,N,E_C,\alpha)$
- integral_flux: $E_\text{bounds,add} \cdot L_\mathrm{diff}(E_\text{bounds,lo}, N, E_C, \alpha) + 4 \cdot L_\mathrm{diff}(E_\text{bounds,half}, N, E_C, \alpha) + L_\mathrm{diff}(E_\text{bounds,hi}, N, E_C, \alpha)$

### flux_ep_correlation
- data: usual
- transformed data:
    - ebounds_half
    - ebounds_add
    - $N_\text{total channels}$
    - log_zp1_grb
    - log_dl2_grb
- parameters:
    - $\alpha$
    - $\beta$
    - $E_\mathrm{peak}$
    - $\gamma$
    - $\delta$
- transformed parameters:
    - $\gamma$
    - $\delta$
    - $\log L$
    - expected model counts
- model:
    - $\alpha \sim \mathcal{N}(-1,0.5)$
    - $\beta \sim \mathcal{N}(-3,1)$
    - $\log E_\mathrm{peak} \sim \mathcal{N}(2,1)$
    - $\gamma_\mu \sim \mathcal{N}(1.5,1)$
    - $\delta_\mu \sim \mathcal{N}(0,4)$
    - $\gamma_\mathrm{off} \sim \mathcal{N}(0,1)$
    - $\delta_\mathrm{off} \sim \mathcal{N}(0,1)$
    - $\gamma_\sigma \sim \mathcal{N}(0,5)$
    - $\delta_\sigma \sim \mathcal{N}(0,5)$
    - log_like / target: ???
- generated quantities: ???

### gbm_csv
- data: usual
- transformed data:
    - ebounds_half
    - ebounds_add
    - out: extract non-zero values from response matrix

### gbm_sum
- data: usual
- transformed data:
    - ebounds_half
    - ebounds_add
    - $N_\text{total channels}$
- parameters:
    - $\alpha$
    - $\beta$
    - $\log E_\mathrm{peak}$
    - $\log L$
- transformed parameters:
    - expected model counts
- model:
    - $\alpha \sim \mathcal{N}(-1,0.5)$
    - $\beta \sim \mathcal{N}(-3,1)$
    - $\log E_\mathrm{peak} \sim \mathcal{N}(2,1)$
    - $\log L \sim \mathcal{N}(-6,2)$
    - log_like / target: ???
- generated quantities: ???

### gbm_sum
- data: usual
- transformed data:
    - ebounds_half
    - ebounds_add
    - $N_\text{total channels}$
- parameters:
    - $\alpha$
    - $\beta$
    - $E_\mathrm{peak}$
    - $L$
- transformed parameters:
    - expected model counts
- model:
    - $\alpha \sim \mathcal{N}(-1,0.5)$
    - $\beta \sim \mathcal{N}(-3,1)$
    - $E_\mathrm{peak} \sim \mathcal{N}(500,500)$
    - $L \sim \mathcal{N}(10^{-6},10^{-2})$
    - log_like / target: ???
- generated quantities: ???

### old stuff
- ggrb_int_cpl:
- band_precalculation: compare!!

### pgstat
- background_model:
- pgstat:

### simple_cpl
- data: usual
- transformed data:
    - ebounds_half
    - ebounds_add
    - $N_\text{total channels}$
- parameters:
    - $\alpha$
    - $\log E_\mathrm{peak}$
    - $\log L$
- transformed parameters:
    - $E_\mathrm{peak}$
    - expected model counts
- model:
    - $\alpha \sim \mathcal{N}(-1,0.5)$
    - $\log E_\mathrm{peak} \sim \mathcal{N}(2,1)$
    - $\log L \sim \mathcal{N}(-7,1)$
    - log_like / target: ???
- generated quantities: empty

### simple_fit
- data: usual
- transformed data:
    - ebounds_half
    - ebounds_add
    - $N_\text{total channels}$
    - log_zp1_grb
    - log_dl2_grb
- parameters:
    - $\alpha$
    - $\beta$
    - $\log E_\mathrm{peak}$
    - $\log L$
- transformed parameters:
    - $E_\mathrm{peak}$
    - expected model counts
- model:
    - $\alpha \sim \mathcal{N}(-1,0.5)$
    - $\beta \sim \mathcal{N}(-3,1)$
    - $\log E_\mathrm{peak} \sim \mathcal{N}(2,1)$
    - $\log L \sim \mathcal{N}(-7,1)$
    - log_like / target: ???
- generated quantities: ???



## Synthetic populations
### aux_samplers
Samples decay time, duration, $E_\mathrm{peak}$ and $L$



## utils
### ogip2stan
- GRBDatum: single GRB with data
    - can export to: plugin, hdf5
    - can read from: ogip/FITS, hdf5
- GRBInterval: time interval consisting of all the detectors
    - can export to: plugin, hdf5
    - can read from: dict/YAML (reads ogip/FITS), hdf5
- GRBData: all intervals from a single GRB
    - can export to: plugin, hdf5
    - can read from: dict/YAML (reads ogip/FITS), hdf5
- DataSet: 
    - can export to: plugin, hdf5, Stan dict
    - can read from: dict/YAML (reads ogip/FITS), hdf5

### sim2fits
- GRBProcessor: retrieves light curves from GRB and exports ogip
- AnalysisBuilder: imports survey and processes GRBS
    - can export to: YAML



## examples
### corr_cpl
TODO

### single_grb
- generates data:
    1. defines samplers for parameters
    2. gets samplers from utils/aux_samplers
    3. observes samplers
    4. writes to data/single_grb.h5
- simulates single GRB to data/test_grb.h5
- simulates universe o single_grb to data/survey.h5