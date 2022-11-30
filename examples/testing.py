import os
import sys

from cosmogrb.universe.survey import Survey

parent_dir = os.path.abspath("..")
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

import arviz as av

from zusammen import AnalysisBuilder, DataSet
from zusammen.stan_models.stan_model import get_model

# survey = Survey.from_file("data/survey.h5")
# ab = AnalysisBuilder(survey, use_bb=True)
# ab.write_yaml("test_proc.yml")
ds = DataSet.from_yaml("test_proc.yml")


m = get_model("cpl_simple_chunked")
m.clean_model()
m.build_model(opt_exp=True)


data = ds.to_stan_dict()

n_threads = 16
n_chains = 2

fit = m.model.sample(
    data=data,
    parallel_chains=n_chains,
    chains=n_chains,
    # inits= {'alpha':-1.},
    threads_per_chain=n_threads,
    seed=1234,
    iter_warmup=1000,
    iter_sampling=500,
    max_treedepth=15,
    adapt_delta=0.9,
    show_progress=True,
)


res = av.from_cmdstanpy(fit)
res.to_netcdf("inference_data/testing_wo_gc.nc")
