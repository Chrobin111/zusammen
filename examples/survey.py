import numpy as np

import popsynth
from popsynth.aux_samplers.trunc_normal_aux_sampler import TruncatedNormalAuxSampler
from popsynth.aux_samplers.lognormal_aux_sampler import LogNormalAuxSampler

import os, sys

parent_dir = os.path.abspath("..")
if parent_dir not in sys.path:
    sys.path.append(parent_dir)
from zusammen.synthetic_populations.aux_samplers import (
    TDecaySampler,
    DurationSampler,
    LumSampler,
    EpeakObsSampler,
)
from corr_cpl.corr_cpl_grb import GBMGRB_CORR_CPL
from corr_cpl.corr_cpl_universe import GBM_CORR_CPL_Universe


pop_gen = popsynth.populations.SFRPopulation(
    r0=2 / 20, rise=1.0, decay=4.0, peak=1.5, r_max=7.0, a=0.1
)

t90 = LogNormalAuxSampler(name="t90", observed=False)
t90.mu = np.log10(10)
t90.tau = 0.25

ep = LogNormalAuxSampler(name="ep", observed=False)
ep.mu = np.log(300)
ep.tau = 0.4

alpha = TruncatedNormalAuxSampler(name="alpha", observed=False)
alpha.lower = -1.5
alpha.upper = 0.0
alpha.mu = -1
alpha.tau = 0.25

ep_tau = TruncatedNormalAuxSampler(name="ep_tau", observed=False)
ep_tau.lower = -2
ep_tau.upper = -1
ep_tau.mu = -1.5
ep_tau.tau = 0.25

nrest = LogNormalAuxSampler(name="nrest", observed=False)
nrest.mu = np.log(1e52)
nrest.tau = 0.1

gamma = TruncatedNormalAuxSampler(name="gamma", observed=False)
gamma.mu = 1.5
gamma.tau = 0.5
gamma.lower = 1.0
gamma.upper = 2.0

duration = DurationSampler()
obs_lum = LumSampler()
ep_obs = EpeakObsSampler()

duration.set_secondary_sampler(t90)
obs_lum.set_secondary_sampler(ep)
obs_lum.set_secondary_sampler(nrest)
obs_lum.set_secondary_sampler(gamma)
ep_obs.set_secondary_sampler(ep)

pop_gen.add_observed_quantity(duration)
pop_gen.add_observed_quantity(obs_lum)
pop_gen.add_observed_quantity(ep_obs)
pop_gen.add_observed_quantity(alpha)
pop_gen.add_observed_quantity(ep_tau)

population = pop_gen.draw_survey()
population.writeto("data/single_grb.h5")


universe = GBM_CORR_CPL_Universe("data/single_grb.h5", save_path="data")
universe.go(client=None)
universe.save("data/survey.h5")
