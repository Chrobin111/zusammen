import os
import sys

import arviz as av
import numpy as np
from threeML import update_logging_level

update_logging_level("FATAL")

parent_dir = os.path.abspath("..")
if parent_dir not in sys.path:
    sys.path.append(parent_dir)
from zusammen import DataSet
from zusammen.stan_models.stan_model import get_model

inference_folder = "inference/"
data_folder = "simulation/"
data_name = "data_2_sig_5"

model_name = "cpl_simple_chunked_gc"
inference_name = "simulated_2_sig_5_1000"

ds = DataSet.from_hdf5_file(data_folder + data_name + ".h5")
data = ds.to_stan_dict()

m = get_model(model_name)
m.clean_model()
m.build_model(opt_exp=True)

n_threads = 12
n_chains = 2
n_warmup = 2000
n_sampling = 1000

fit = m.model.sample(
    data=data,
    chains=n_chains,
    parallel_chains=n_chains,
    threads_per_chain=n_threads,
    inits={
        "alpha": -1 * np.ones(data["N_intervals"]),
        "log_ec": 2 * np.ones(data["N_intervals"]),
        "gamma_sig_meta": 1,
        "log_Nrest_sig_meta": 1,
        "gamma_mu_meta": 1.5,
        "log_Nrest_mu_meta": 52,
        "gamma": 1.5 * np.ones(data["N_grbs"]),
        "log_Nrest": 52 * np.ones(data["N_grbs"]),
        # 'gamma': 1.5,
        # 'log_Nrest': 52,
    },  # type: ignore
    iter_warmup=n_warmup,
    iter_sampling=n_sampling,
    max_treedepth=15,
    adapt_delta=0.99,
    step_size=0.1,
    show_progress=True,
    refresh=1,
)

res = av.from_cmdstanpy(fit)
res.to_netcdf(inference_folder + inference_name + ".nc")
