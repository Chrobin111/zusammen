# Files
## Stan functions
### band_grb
- ggrb_into_pl:  
    $pl(\alpha,\beta,E_c,E_\mathrm{min},E_\mathrm{max}) = \begin{cases}
        ((\alpha - \beta)^{\alpha-\beta} \cdot \frac{\mathrm{e}^{\beta-\alpha}}{E_c^\beta})/(2\beta) \cdot (E_\mathrm{max}^{2 + \beta} - E_\mathrm{min}^{2 + \beta}), & \beta \neq 2 \\
        (\alpha - \beta)^{\alpha-\beta} \cdot \frac{\mathrm{e}^{\beta-\alpha}}{E_c^\beta} \cdot \log\frac{E_\mathrm{max}}{E_\mathrm{min}}, & \mathrm{else}
    \end{cases}$
- ggrb_int_cpl:  
    $cpl(\alpha,E_c,E_\mathrm{min},E_\mathrm{max}) = -E_c^2 \cdot (\Gamma(2 + \alpha, E_\mathrm{max}/E_c) - \Gamma(2 + \alpha, E_\mathrm{min}/E_c)) \cdot \Gamma(2 + \alpha)$
- band_precalculation:  
    norm $= \begin{cases}
    F / (c(\alpha,E_c,E_\mathrm{min},E_\mathrm{split}) + p(\alpha,\beta,E_c,E_\mathrm{split},E_\mathrm{max})), & E_\mathrm{min} \leq E_\mathrm{split} \leq E_\mathrm{max} \\
    F / pl(\alpha,\beta,E_c,E_\mathrm{split},E_\mathrm{max}), & E_\mathrm{split} < E_\mathrm{min}
    \end{cases}$  
    $pre = (\alpha-\beta)^{\alpha-\beta} \cdot \mathrm{e}^{\beta-\alpha}$
- differential_flux:  
    $F_\mathrm{diff}^i = \mathrm{norm} \cdot o^i$  
    $o^i = \begin{cases}
        (E/E_c)^\alpha \cdot \mathrm{e}^{-E/E_c}, & E < E\mathrm{split} \\
        pre \cdot (E/E_c)^\beta, & \mathrm{else}
    \end{cases}$
- integral_flux:  
    $F_\mathrm{int} = E_\mathrm{bounds,add} \cdot (F_\mathrm{diff}(E_\mathrm{bounds,lo},\dots) + 4 \cdot F_\mathrm{diff}(E_\mathrm{bounds,half},\dots) + F_\mathrm{diff}(E_\mathrm{bounds,hi},\dots))$

### cpl_interval_fold:
- partial_log_like:  
    $N_\mathrm{exp} = R \cdot F_\mathrm{int}(E_\mathrm{bounds,lo},E_\mathrm{bounds,hi},E_\mathrm{bounds,add},E_\mathrm{bounds,half},K,E_c,\alpha) \cdot t_\mathrm{exposure}$  
    $L = \sum_i PG(N_\mathrm{obs}^i,N_\mathrm{back}^i,\sigma^i,N_\mathrm{exp}^i,idx_0^i,idx^i)$

### cpl
- ggrb_int_cpl:
    $cpl(\alpha,E_c,E_\mathrm{min},E_\mathrm{max}) = -E_c^2 \cdot (\Gamma(2 + \alpha, E_\mathrm{max}/E_c) - \Gamma(2 + \alpha, E_\mathrm{min}/E_c)) \cdot \Gamma(2 + \alpha)$
- cpl: $N \cdot (E/100)^\alpha \cdot \mathrm{e}^{-E/E_C}$
- cpl_indi: $K \cdot (E/100)^\alpha \cdot \mathrm{e}^{-E/E_C}$
- cpl_flux_integrand: $x \cdot$ cpl_indi $(x,\theta)$
- differential_flux: $L_\mathrm{diff}(E,N,E_C,\alpha) = \mathrm{cpl}(E,N,E_C,\alpha)$
- integral_flux: $E_\text{bounds,add} \cdot L_\mathrm{diff}(E_\text{bounds,lo}, N, E_C, \alpha) + 4 \cdot L_\mathrm{diff}(E_\text{bounds,half}, N, E_C, \alpha) + L_\mathrm{diff}(E_\text{bounds,hi}, N, E_C, \alpha)$

### pgstat
- background_model: $b = 0.5 \cdot \sqrt{MB^2 - 2 \sigma^2 \cdot (MB - 2 N_\mathrm{obs}) + \sigma^2} + N_\mathrm{back} - N_\mathrm{exp} - \sigma^2$
- pgstat:  
$L(idx \neq 0) = -\frac{(BM - N_\mathrm{back})^2}{2 \sigma^2} + N_\mathrm{obs} \cdot \log(BM + N_\mathrm{exp}) - BM - N_\mathrm{exp} + \log\Gamma(idx + 1) - 0.5 \cdot \log(2\pi) - \log(\sigma)$  
$L(idx=0) = N_\mathrm{obs} \log(N_\mathrm{exp} - N_\mathrm{exp} + \log\Gamma(idx + 1))$

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

### cpl_simple_chunked:
- data: usual
- transformed data:
    - ebounds_half
    - ebounds_add
    - $N_\text{total channels}$
    - all_N
    - keV to erg
    - erg to keV
    - $E_\mathrm{min}$
    - $E_\mathrm{max}$
    - $x_r$
    - $x_i$
- parameters:
    - $\alpha$
    - $\log E_C$ (observed energy)
    - $\log F$
- transformed parameters:
    - $E_C$
    - $F$
    - $K = L \cdot \left( \int_{10}^{1000} dx\ x \cdot K \cdot (E/100)^\alpha \cdot \mathrm{e}^{-E/100} \right)^{-1}$
- model:
    - $\log y \sim \mathcal{N}(0,1)$
    - $\log y_\sigma \sim \mathcal{N}(0,1)$
    - $\alpha \sim \mathcal{N}(-1,0.5)$
    - $\log E_C \sim \mathcal{N}(2,1)$
    - target: sums over partial_log_like

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