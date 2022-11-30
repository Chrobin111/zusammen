from pathlib import Path

import astropy.io.fits as fits
from threeML.plugins.OGIPLike import OGIPLike


def get_number_of_intervals(base_path: Path, grb: str):

    grb_path = base_path / grb

    dets = [x.stem for x in grb_path.glob("*.pha") if "_bak" not in x.name]

    with fits.open(grb_path / f"{dets[0]}.pha") as f:

        n_intervals = int(f[1].header["NAXIS2"])

    return n_intervals


def get_data(base_path: Path, grb: str, interval: int):

    grb_path = base_path / grb

    dets = [x.stem for x in grb_path.glob("*.pha") if "_bak" not in x.name]

    plugins = []

    for det in dets:

        o = OGIPLike(
            det,
            observation=grb_path / f"{det}.pha",
            background=grb_path / f"{det}_bak.pha",
            response=grb_path / f"{det}.rsp",
            spectrum_number=interval + 1,
        )

        if det.startswith("n"):

            o.set_active_measurements("8-900")

        elif det.startswith("b"):

            o.set_active_measurements("250-30000")

        else:

            o.set_active_measurements("50000-100000")

        plugins.append(o)

    return plugins
