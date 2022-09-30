import os
import pkg_resources

import cmdstanpy

_available_models = {}

_available_models["cpl_simple_chunked"] = "cpl_simple_chunked.stan"
_available_models["cpl_simple_chunked_centered"] = "cpl_simple_chunked_centered.stan"
_available_models["cpl_simple_chunked_gc"] = "cpl_simple_chunked_gc.stan"
_available_models[
    "cpl_simple_chunked_gc_vectorized"
] = "cpl_simple_chunked_gc_vectorized.stan"
_available_models["cpl_simple_chunked_gc_cauchy"] = "cpl_simple_chunked_gc_cauchy.stan"
_available_models[
    "cpl_simple_chunked_gc_relaxed"
] = "cpl_simple_chunked_gc_relaxed.stan"
_available_models[
    "cpl_simple_chunked_combined"
] = "cpl_simple_chunked_combined.stan"


class StanModel(object):
    def __init__(self, name, stan_file):

        self._name = name
        self._stan_file = pkg_resources.resource_filename(
            "zusammen", os.path.join("stan_models", stan_file)
        )

        self._model = None

    def build_model(
        self, recompile: bool = False, opt: bool = False, opt_exp: bool = False
    ):
        """
        build the stan model

        :returns:
        :rtype:

        """

        compile = "force" if recompile else True
        cpp_options = dict(STAN_THREADS=True)
        stanc_options = dict()
        if opt:
            stanc_options["O1"] = True
        if opt_exp:
            stanc_options["Oexperimental"] = True

        self._model = cmdstanpy.CmdStanModel(
            stan_file=self._stan_file,
            model_name=self._name,
            compile=compile,
            cpp_options=cpp_options,
            stanc_options=stanc_options,
        )

    @property
    def model(self):
        return self._model

    def clean_model(self):
        """
        Clean the model bin file
        to allow for compiling

        :returns:
        :rtype:

        """

        bin = self._stan_file.replace(".stan", "")
        if os.path.exists(bin):
            os.remove(bin)
            print(f"Removed binary {bin}.")
        else:
            print(f"Binary {bin} does not exist.")


def get_model(model_name):
    """
    Retrieve the stan model

    :param model_name:
    :returns:
    :rtype:

    """

    assert (
        model_name in _available_models
    ), f"please chose {','.join(x for x in _available_models.keys()) }"

    return StanModel(model_name, _available_models[model_name])
