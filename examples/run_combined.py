import numpy as np

# plt.style.use("mike")
import warnings

warnings.simplefilter("ignore")
warnings.filterwarnings("ignore")

import os, sys

parent_dir = os.path.abspath("..")
if parent_dir not in sys.path:
    sys.path.append(parent_dir)
from zusammen.stan_models.stan_model import get_model
from zusammen import AnalysisBuilder, DataSet
from zusammen.spectral_plot import display_posterior_model_counts

from threeML import update_logging_level

import arviz as av


update_logging_level("FATAL")


# ds = DataSet.from_yaml("test_proc_sig_10_red.yml")
ds = DataSet.from_hdf5_file("sgrb.h5")

# cmdstanpy.set_cmdstan_path("/home/bjorn/general_sw/cmdstan-2.29.0/")

m = get_model("cpl_simple_chunked_combined")

# m.clean_model()

m.build_model()  # opt_exp=True)

data = ds.to_stan_dict()

n_threads = 2
n_chains = 2

fit = m.model.sample(
    data=data,
    parallel_chains=n_chains,
    chains=n_chains,
    inits={
        "alpha": -1.0 * np.ones(data["N_intervals"]),
        "log_ec": 2 * np.ones(data["N_intervals"]),
        # "log_K": -1 * np.ones(data["N_intervals"]),
        "log_energy_flux_mu_raw": 0,
        "log_energy_flux_sigma": 1,
        "log_energy_flux_raw": np.zeros(data["N_intervals"]),
        # "log_energy_flux": -6 * np.ones(data["N_intervals"])},
    },
    threads_per_chain=n_threads,
    seed=1234,
    iter_warmup=1000,
    iter_sampling=500,
    max_treedepth=15,
    adapt_delta=0.9,
    # show_progress=True,
    show_console=True
    # output_dir=".stan/",
)

res = av.from_cmdstanpy(fit)
res.to_netcdf("inference_data/analytic.nc")
