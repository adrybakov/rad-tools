r"""
All tools from the package.


|version|
"""

__version__ = "0.5.29"


from . import score
from .score import *

__all__ = []
__all__.extend(score.__all__)


if __name__ == "__main__":
    local_variables = dict(locals())
    for i in local_variables:
        print(i)
