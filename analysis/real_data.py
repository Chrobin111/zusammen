import os
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Optional

import numpy as np
import yaml
from threeML import FermiGBMBurstCatalog, TimeSeriesBuilder

real_data_dir = os.path.abspath("real_data")
if real_data_dir not in sys.path:
    sys.path.append(real_data_dir)
from get_data import get_data, get_number_of_intervals


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


class SynchSample:
    def __init__(
        self,
        base_path,
        grb_names,
        sig_min: Optional[float] = None,
        interval_min: int = 1,
    ):
        if isinstance(base_path, str):
            self._base_path = Path(base_path)
        else:
            self._base_path = base_path
        self._grb_names = grb_names
        self._sig_min = sig_min
        self._interval_min = interval_min
        self._grbs = OrderedDict()

        self._threeml_process()

    def _threeml_process(self):
        for grb in self._grb_names:
            n_intervals = get_number_of_intervals(self._base_path, grb)
            intervals = []
            all_detectors = []

            for interval in range(n_intervals):
                spectra = get_data(self._base_path, grb, interval)
                detectors_per_interval = []

                significant = False
                for spectrum in spectra:
                    detectors_per_interval.append(spectrum._name)
                    if spectrum.significance >= self._sig_min:
                        significant = True
                if significant:
                    intervals.append(interval)
                    all_detectors.append(detectors_per_interval)

            if len(intervals) < self._interval_min:
                continue

            for dets1, dets2 in zip(all_detectors[1:], all_detectors):
                assert dets1 == dets2

            self._grbs[grb] = OrderedDict()
            self._grbs[grb]["dir"] = str((self._base_path / grb).absolute())
            self._grbs[grb]["interval_ids"] = intervals
            self._grbs[grb]["detectors"] = OrderedDict()
            for det in all_detectors[0]:
                active = ""
                if det.startswith("n"):
                    active = "10-900"
                elif det.startswith("b"):
                    active = "250-10000"
                self._grbs[grb]["detectors"][det] = active

    def write_yaml(self, file_name: str) -> None:
        with open(file_name, "w") as f:
            yaml.dump(self._grbs, stream=f, default_flow_style=False)
