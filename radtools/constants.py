from math import pi

__all__ = ["todegrees", "toradians"]

RED = "#FF4D67"
GREEN = "#58EC2E"
ORANGE = "#F7CB3D"
BLUE = "#274DD1"
PURPLE = "#DC5CFF"

# Constants are usually defined in uppercase
# but these two are intentionally defined in lowercase

todegrees = 180.0 / pi
"""Convert radians to degrees.
"""

toradians = pi / 180.0
"""Convert degrees to radians.
"""
