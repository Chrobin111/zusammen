import collections

from pathlib import Path

import numpy as np
import yaml

from cosmogrb.io.gbm_fits import grbsave_to_gbm_fits
from cosmogrb.universe.survey import Survey
from cosmogrb.utils.file_utils import if_directory_not_existing_then_make
from threeML import TimeSeriesBuilder
from threeML.utils.interval import Interval
from typing import Union


class GRBProcessor(object):
    def __init__(
        self,
        gbm_grb,
        n_nai_to_use: int = 3,
        use_bb: bool = False,
        sig_min: Union[float, None] = None,
        all_above_limit: bool = False,
    ):
        """
        :param gbm_grb:
        :param n_nai_to_use:
        :returns:
        :rtype:

        """

        self._grb_save = gbm_grb

        assert n_nai_to_use > 0, "yo use some detectors"
        self._n_nai_to_use: int = int(n_nai_to_use)

        self._use_bb: bool = use_bb

        self._sig_min: Union[float, None] = sig_min

        self._all_above_limit: bool = all_above_limit

        self._config_dict = collections.OrderedDict()

        self._config_dict["z"] = float(self._grb_save.z)

        if_directory_not_existing_then_make(self._grb_save.name)

        # gets the light curves we want
        self._setup_order_by_distance()

        self._create_fits_files()

        self._threeml_process()

    def _setup_order_by_distance(self):

        # now we will go through the lightcurves
        # collect

        angular_distances = []
        bgo_angular_distance = 1000
        bgo_det = ""
        lc_names = []

        for name, det in self._grb_save.items():

            lc = det["lightcurve"]

            if name.startswith("n"):

                lc_names.append(str(name))

                angular_distances.append(lc.extra_info["angle"])

            else:

                if lc.extra_info["angle"] < bgo_angular_distance:
                    bgo_angular_distance = lc.extra_info["angle"]
                    bgo_det = str(name)

        angular_distances = np.array(angular_distances)
        lc_names = np.array(lc_names)

        # now attach the lc_names sorted by
        # the detector distance to the GRB

        idx = angular_distances.argsort()

        self._lc_names = list(lc_names[idx][: self._n_nai_to_use])
        self._lc_names.append(bgo_det)
        self._lc_names = [str(x) for x in self._lc_names]

    def _create_fits_files(self):

        self._fits_files = grbsave_to_gbm_fits(
            self._grb_save,
            destination=self._grb_save.name,
            detectors=self._lc_names,
        )

    def _threeml_process(self):

        self._config_dict["dir"] = str(Path(self._grb_save.name).absolute())

        det_dic = {}

        for i, name in enumerate(self._lc_names):

            if name.startswith("n"):

                selection = "10-900"

            else:

                selection = "250-10000"

            det_dic[name] = selection

            ts = TimeSeriesBuilder.from_gbm_tte(
                name=name,
                tte_file=self._fits_files[name]["tte"],
                rsp_file=self._fits_files[name]["rsp"],
                poly_order=0,
                verbose=False,
            )

            ts.set_background_interval(
                "-20--5",
                f"{self._grb_save.duration + 5}-{self._grb_save.duration + 20}",
            )

            # for now do nothing else
            if self._use_bb:

                if i < 1:
                    ts.create_time_bins(
                        -25,
                        self._grb_save.duration + 1,
                        method="bayesblocks",
                        p0=0.1,
                    )
                    bins_to_use = ts
                    above_limit = np.ones(
                        (len(self._lc_names), len(ts.bins)), dtype=bool
                    )
                else:
                    ts.read_bins(bins_to_use)

                intervals = ts.bins

                # check for non-significant intervals
                if self._sig_min is not None:
                    sig = ts.significance_per_interval
                    above_limit[i] = sig > self._sig_min

                ts.write_pha_from_binner(
                    file_name=Path(self._grb_save.name) / name,
                    start=-25,
                    stop=self._grb_save.duration + 1,
                    # inner=True,
                    force_rsp_write=True,
                    overwrite=True,
                )

                if len(intervals) > 1:
                    fig = ts.view_lightcurve(use_binner=True)

                    fig.savefig(
                        f"{Path(self._grb_save.name) / name}_lc.png",
                        bbox_inches="tight",
                    )

                intervals_all = np.arange(len(intervals))

            else:

                ts.set_active_time_interval(f"0-{self._grb_save.duration}")

                plugin = ts.to_spectrumlike()

                plugin.write_pha(
                    filename=Path(self._grb_save.name) / name,
                    force_rsp_write=True,
                    overwrite=True,
                )

        if self._all_above_limit:
            interval_ids = np.argwhere(np.all(above_limit, axis=0)).flatten()
        else:
            interval_ids = np.argwhere(np.any(above_limit, axis=0)).flatten()

        self._config_dict["interval_ids"] = interval_ids.tolist()

        self._config_dict["detectors"] = det_dic

        for k, v in self._fits_files[name].items():
            Path(v).unlink()

    @property
    def yaml_params(self):

        return self._config_dict


class AnalysisBuilder(object):
    def __init__(
        self,
        survey_file,
        use_all: bool = False,
        use_bb: bool = False,
        sig_min: Union[float, None] = None,
        all_above_limit: bool = False,
        intervals_min: int = 1,
    ):

        if isinstance(survey_file, str):

            self._survey = Survey.from_file(survey_file)

        else:

            assert isinstance(survey_file, Survey)

            self._survey = survey_file

        if intervals_min < 1:
            print(
                "An event has to have at least one interval.\ninterval_min has been increased to 1."
            )
            intervals_min = 1

        self._config_dict = collections.OrderedDict()
        for k, v in self._survey.items():

            process = GRBProcessor(
                v.grb, use_bb=use_bb, sig_min=sig_min, all_above_limit=all_above_limit
            )

            if len(process.yaml_params["interval_ids"]) > intervals_min:
                self._config_dict[k] = process.yaml_params

    def write_yaml(self, file_name: str) -> None:
        """TODO describe function

        :param file_name:
        :type file_name: str
        :returns:

        """

        with open(file_name, "w") as f:
            yaml.dump(self.yaml_params, stream=f, default_flow_style=False)

    @property
    def yaml_params(self):
        return self._config_dict
