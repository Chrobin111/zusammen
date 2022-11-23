import os
from collections import OrderedDict

import numpy as np
import yaml
from threeML import FermiGBMBurstCatalog, TimeSeriesBuilder


class AnalysisBuilder:
    def __init__(self, download_name, data_directory, sig_min, interval_min) -> None:
        self._gbm_catalog = FermiGBMBurstCatalog()
        self._yaml_dict = OrderedDict()

        self._dload = dict()
        with open(data_directory + "/" + download_name, "r") as f:
            self._dload = yaml.load(f, Loader=yaml.SafeLoader)

        self._data_directory = data_directory
        self._sig_min = sig_min
        self._interval_min = interval_min

        self._threeml_process()

    def _threeml_process(self):
        for grb, dload in self._dload.items():

            selection = dict()
            source_time = dload["source_time"]

            self._gbm_catalog.query_sources(grb)
            det_info = self._gbm_catalog.get_detector_information()[grb]

            for i, det in enumerate(det_info["detectors"]):

                if det.startswith("n"):
                    selection[str(det)] = "10-900"
                else:
                    selection[str(det)] = "250-10000"

                ts = TimeSeriesBuilder.from_gbm_tte(
                    det,
                    tte_file=dload[det]["tte"].replace(".gz", ""),
                    rsp_file=dload[det]["rsp"],
                    poly_order=0,
                    verbose=False,
                )

                ts.set_background_interval(*det_info["background"]["full"].split(","))

                if i < 1:
                    ts.create_time_bins(
                        *source_time,
                        method="bayesblocks",
                        p0=0.1,
                    )
                    bins_to_use = ts
                    above_limit = np.ones(
                        (len(det_info["detectors"]), len(ts.bins)), dtype=bool
                    )
                else:
                    ts.read_bins(bins_to_use)

                sig = ts.significance_per_interval
                above_limit[i] = sig > self._sig_min

                ts.write_pha_from_binner(
                    file_name=self._data_directory + "/" + grb + "/" + det,
                    start=source_time[0],
                    stop=source_time[1],
                    force_rsp_write=True,
                    overwrite=True,
                )

                interval_ids = np.argwhere(np.any(above_limit, axis=0)).flatten()

            if len(interval_ids) < self._interval_min:
                self._yaml_dict[grb] = OrderedDict()

                self._yaml_dict[grb]["z"] = float(dload["z"])
                self._yaml_dict[grb]["dir"] = os.path.abspath(
                    self._data_directory + "/" + grb
                )
                self._yaml_dict[grb]["interval_ids"] = interval_ids.tolist()
                self._yaml_dict[grb]["detectors"] = {i: j for i, j in selection.items()}

    def write_yaml(self, file_name: str) -> None:
        with open(file_name, "w") as f:
            yaml.dump(self._yaml_dict, stream=f, default_flow_style=False)
