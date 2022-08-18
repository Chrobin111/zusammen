from cosmogrb.universe.survey import Survey

import os, sys

parent_dir = os.path.abspath("..")
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from zusammen.stan_models.stan_model import get_model
from zusammen import AnalysisBuilder, DataSet

import arviz as av


survey = Survey.from_file("data/survey.h5")
ab = AnalysisBuilder(survey, use_bb=True)
ab.write_yaml("test_proc.yml")
ds = DataSet.from_yaml("test_proc.yml")


m = get_model("cpl_simple_chunked_gc")
m.clean_model()
m.build_model()


data = ds.to_stan_dict()

n_threads = 5
n_chains = 2

fit = m.model.sample(
    data=data,
    parallel_chains=n_chains,
    chains=n_chains,
    # inits= {'alpha':-1.},
    threads_per_chain=n_threads,
    seed=1234,
    iter_warmup=2000,
    iter_sampling=1000,
    max_treedepth=12,
    show_progress=True,
)


res = av.from_cmdstanpy(fit)
res.to_netcdf("inference_data/test_inference_data_energy_flux_6.nc")
