import astromodels.functions.numba_functions as nb_func
import astropy.units as astropy_units
from astromodels.functions.function import Function1D, FunctionMeta


class Cutoff_powerlaw_prime(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A power law multiplied by an exponential cutoff

    latex : $ K~\frac{x}{xc}^{index}~\exp{-x/xc} $

    parameters :

        K :

            desc : Normalization (differential flux at the pivot value)
            initial value : 1.0
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e3
            delta : 0.1

        index :

            desc : Photon index
            initial value : -2
            min : -10
            max : 10

        xc :

            desc : Cutoff energy
            initial value : 10.0
            transformation : log10
            min: 1.0

    """

    def _set_units(self, x_unit, y_unit):
        # The index is always dimensionless
        self.index.unit = astropy_units.dimensionless_unscaled

        self.xc.unit = x_unit

        # The normalization has the same units as the y
        self.K.unit = y_unit

    # noinspectionq PyPep8Naming
    def evaluate(self, x, K, index, xc):

        if isinstance(x, astropy_units.Quantity):
            index_ = index.value
            K_ = K.value
            xc_ = xc.value
            x_ = x.value

            unit_ = self.y_unit

        else:
            unit_ = 1.0
            K_, x_, index_, xc_ = K, x, index, xc

        result = nb_func.cplaw_eval(x_, K_, xc_, index_, xc_)

        return result * unit_
