from typing import Tuple

import arviz as av
import numpy as np


class PPC:
    def __init__(self, data: dict, res: av.InferenceData):
        self.N_intervals = data["N_intervals"]
        self.N_grbs = data["N_grbs"]
        self.chains = res.posterior.gamma.shape[0]
        self.draws = res.posterior.gamma.shape[1]

        self.alpha = np.zeros((self.N_intervals, self.chains, self.draws))
        self.log_ec = np.zeros((self.N_intervals, self.chains, self.draws))
        self.K = np.zeros((self.N_intervals, self.chains, self.draws))
        self.grb_id = data["grb_id"]

        self.cbounds_lo = data["cbounds_lo"]
        self.cbounds_hi = data["cbounds_hi"]
        self.ebounds_lo = data["ebounds_lo"]
        self.ebounds_hi = data["ebounds_hi"]

        self.response = data["response"]
        self.exposure = data["exposure"]
        self.N_chan = data["N_chan"]
        self.N_echan = data["N_echan"]

        for id in range(self.N_intervals):
            self.alpha[id] = res.posterior.alpha.values[:, :, id]
            self.log_ec[id] = res.posterior.log_ec.values[:, :, id]
            self.K[id] = res.posterior.K.values[:, :, id]

    @staticmethod
    def cpl(E, K, ec, alpha):
        return K * (E / ec) ** alpha * np.exp(-E / ec)

    @staticmethod
    def integral_flux(ebounds_lo, ebounds_hi, K, ec, alpha):
        ebounds_add = (ebounds_hi - ebounds_lo) / 6
        ebounds_half = (ebounds_hi + ebounds_lo) / 2
        return ebounds_add * (
            PPC.cpl(ebounds_lo, K, ec, alpha)
            + 4 * PPC.cpl(ebounds_half, K, ec, alpha)
            + PPC.cpl(ebounds_hi, K, ec, alpha)
        )

    def ppc(
        self,
        chain: int,
        draws: int,
        grb: int,
        interval: int,
        detector: int,
        interval_in_grb=False,
    ) -> Tuple[np.ndarray, np.ndarray]:

        assert chain < self.chains
        assert draws <= self.draws
        assert grb < self.N_grbs

        if interval_in_grb:
            interval_final = np.where(self.grb_id == grb + 1)[0][interval]
        else:
            interval_final = interval
        cenergies = np.append(
            self.cbounds_lo[interval_final, detector],
            self.cbounds_hi[interval_final, detector, -1],
        )

        rng = np.random.default_rng()
        rand = rng.choice(self.draws, size=draws, replace=False)
        ppc_alpha = np.array([self.alpha[interval_final, chain, i] for i in rand])
        ppc_log_ec = np.array([self.log_ec[interval_final, chain, i] for i in rand])
        ppc_K = np.array([self.K[interval_final, chain, i] for i in rand])

        ppc_expected_model_counts = np.zeros(
            (draws, self.N_chan[interval_final, detector])
        )

        for s in range(draws):
            ppc_expected_model_counts[s, : self.N_chan[interval_final, detector]] = (
                self.response[
                    interval_final,
                    detector,
                    : self.N_chan[interval_final, detector],
                    : self.N_echan[interval_final, detector],
                ]
                @ PPC.integral_flux(
                    self.ebounds_lo[
                        interval_final,
                        detector,
                        : self.N_echan[interval_final, detector],
                    ],
                    self.ebounds_hi[
                        interval_final,
                        detector,
                        : self.N_echan[interval_final, detector],
                    ],
                    ppc_K[s],
                    10 ** ppc_log_ec[s],
                    ppc_alpha[s],
                )
            ) * self.exposure[interval_final, detector]

        ppc_sampled_counts = np.random.poisson(ppc_expected_model_counts)

        return cenergies, ppc_sampled_counts

    @staticmethod
    def ppc_summary(ppc_sampled_counts):
        ppc_sampled_counts_mu = np.mean(ppc_sampled_counts, axis=0)
        ppc_sampled_counts_hdi_1s = np.zeros((2, ppc_sampled_counts.shape[1]))
        ppc_sampled_counts_hdi_2s = np.zeros((2, ppc_sampled_counts.shape[1]))

        for i, counts in enumerate(ppc_sampled_counts.T):
            ppc_sampled_counts_hdi_1s.T[i] = av.hdi(counts, 0.6828)
            ppc_sampled_counts_hdi_2s.T[i] = av.hdi(counts, 0.9545)

        return (
            ppc_sampled_counts_mu,
            ppc_sampled_counts_hdi_1s,
            ppc_sampled_counts_hdi_2s,
        )
