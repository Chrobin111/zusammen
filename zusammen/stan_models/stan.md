# Files
## alpha_correlation
- data: usual
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

## band_grb
- ggrb_into_pl:
- ggrb_int_cpl:
- band_precalculation:
- differential_flux:
- integral_flux: 

## big_fuck_correlation
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

## correlation
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

## cpl_interval_fold:
- partial_log_like: 

## cpl_simple_chunked:
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

## cpl
- ggrb_int_cpl:
- cpl: $N \cdot (E/100)^\alpha \cdot \mathrm{e}^{-E/E_C}$
- cpl_indi: $K \cdot (E/100)^\alpha \cdot \mathrm{e}^{-E/E_C}$
- cpl_flux_integrand: $x \cdot \mathrm{cpl_indi}(x,\theta)$
- differential_flux: $L_\mathrm{diff}(E,N,E_C,\alpha) = \mathrm{cpl}(E,N,E_C,\alpha)$
- integral_flux: $E_\text{bounds_add} \cdot L_\mathrm{diff}(E_\text{bounds_lo}, N, E_C, \alpha) + 4 \cdot L_\mathrm{diff}(E_\text{bounds_half}, N, E_C, \alpha) + L_\mathrm{diff}(E_\text{bounds_hi}, N, E_C, \alpha)$

## flux_ep_correlation