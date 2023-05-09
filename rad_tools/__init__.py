r"""
All tools from the package.


|version|
"""

__version__ = "0.6.0"


from . import score
from .score import *
from . import io
from .io import *
from . import exchange
from .exchange import *
from . import kpoints
from .kpoints import *
from . import dos
from .dos import *
from .routines import *

__all__ = []
__all__.extend(score.__all__)
__all__.extend(io.__all__)
__all__.extend(exchange.__all__)
__all__.extend(kpoints.__all__)
__all__.extend(dos.__all__)
__all__.extend(routines.__all__)


if __name__ == "__main__":
    local_variables = dict(locals())
    for i in local_variables:
        print(i)
