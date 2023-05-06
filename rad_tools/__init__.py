r"""
All tools from the package.


|version|
"""

__version__ = "0.6.0"


from . import score
from .score import *
from . import io
from .io import *

__all__ = []
__all__.extend(score.__all__)
__all__.extend(io.__all__)


if __name__ == "__main__":
    local_variables = dict(locals())
    for i in local_variables:
        print(i)
