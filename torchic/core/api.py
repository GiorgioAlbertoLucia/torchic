from torchic.core.dataset import (
    Dataset,
)

from torchic.core import histogram
from torchic.core.histogram import (
    AxisSpec,
    HistLoadInfo
)

from torchic.core.fitter import (
    Fitter,
    fit_TH1,
    fit_by_slices,
)

__all__ = [
    'Dataset',
    'AxisSpec',
    'HistLoadInfo',
    'histogram',
    'Fitter',
    'fit_TH1',
    'fit_by_slices',
]