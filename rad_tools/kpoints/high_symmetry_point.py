from math import tan, cos, sin


class HighSymmetryPoints:
    r"""
    Generator of the high symmetry k-points.

    Parameters
    ----------
    kpoints : dict
        Distionary of the high symmetry k-points and coordinate
        (fractions of the reciprocal vectors).

        .. code-block:: python

            kpoints = {"name" : [xb1, xb2, xb3]}
    """

    def __init__(self,
                 kpoints=None,
                 path=None) -> None:
        self._kpoints = kpoints
        self._path = path

    def cub(self):
        r"""Cubic (CUB, cP)"""

        kpoints = {
            "Gamma": [0, 0, 0],
            "M": [1/2, 1/2, 0],
            "R": [1/2, 1/2, 1/2],
            "X": [0, 1/2, 0]}

        path = [
            ["Gamma", "X", "M", "Gamma", "R", "X"],
            ["M", "R"]]

        return HighSymmetryPoints(kpoints, path)

    def fcc(self):
        r"""Face-centered cubic (FCC, cF)"""

        kpoints = {
            "Gamma": [0, 0, 0],
            "K": [3/8, 3/8, 3/4],
            "L": [1/2, 1/2, 1/2],
            "U": [5/8, 1/4, 5/8],
            "W": [1/2, 1/4, 3/4],
            "X": [1/2, 0, 1/2]}

        path = [
            ["Gamma", "X", "W", "K", "Gamma", "L", "U", "W", "L", "K"],
            ["U", "X"]]

        return HighSymmetryPoints(kpoints, path)

    def bcc(self):
        r"""Body-centered cubic (BCC, cl)"""

        kpoints = {
            "Gamma": [0, 0, 0],
            "H": [1/2, -1/2, 1/2],
            "P": [1/4, 1/4, 1/4],
            "N": [0, 0, 1/2]}

        path = [["Gamma", "H", "N", "Gamma", "P", "H"],
                ["P", "N"]]

        return HighSymmetryPoints(kpoints, path)

    def tet(self):
        r"""Tetragonal (TET, tP)"""

        kpoints = {
            "Gamma": [0, 0, 0],
            "A": [1/2, 1/2, 1/2],
            "M": [1/2, 1/2, 0],
            "R": [0, 1/2, 1/2],
            "X": [0, 1/2, 0],
            "Z": [0, 0, 1/2]}

        path = [
            ["Gamma", "X", "M", "Gamma", "Z", "R", "A", "Z"],
            ["X", "R"],
            ["M", "A"]]

        return HighSymmetryPoints(kpoints, path)

    def bct1(self, a, c):
        r"""Body-centered tetragonal (BCT, tI),

        c < a"""

        eta = (1 + c**2 / a**2) / 4

        kpoints = {
            "Gamma": [0, 0, 0],
            "M": [-1/2, 1/2, 1/2],
            "N": [0, 1/2, 0],
            "P": [1/4, 1/4, 1/4],
            "X": [0, 0, 1/2],
            "Z": [eta, eta, -eta],
            "Z1": [-eta, 1 - eta, eta]}

        path = [
            ["Gamma", "X", "M", "Gamma", "Z", "P", "N", "Z1", "M"],
            ["X", "P"]]

        return HighSymmetryPoints(kpoints, path)

    def bct2(self, a, c):
        r"""Body-centered tetragonal (BCT, tI),

        c > a"""

        eta = (1 + c**2 / a**2) / 4
        zeta = a**2 / (2 * c**2)
        kpoints = {
            "Gamma": [0, 0, 0],
            "N": [0, 1/2, 0],
            "P": [1/4, 1/4, 1/4],
            "Sigma": [-eta, eta, eta],
            "Sigma1": [eta, 1 - eta, -eta],
            "X": [0, 0, 1/2],
            "Y": [-zeta, zeta, 1/2],
            "Y1": [1/2, 1/2, -zeta],
            "Z": [1/2, 1/2, -1/2]}

        path = [
            ["gamma", "X", "Y", "Sigma", "Gamma",
                "Z", "Sigma1", "N", "P", "Y1", "Z"],
            ["X", "P"]]

        return HighSymmetryPoints(kpoints, path)

    def orc(self):
        r"""Orthorhombic (ORC, oP)"""

        kpoints = {
            "Gamma": [0, 0, 0],
            "R": [1/2, 1/2, 1/2],
            "S": [1/2, 1/2, 0],
            "T": [0, 1/2, 1/2],
            "U": [1/2, 0, 1/2],
            "X": [1/2, 0, 0],
            "Y": [0, 1/2, 0],
            "Z": [0, 0, 1/2]}

        path = [
            ["Gamma", "X", "S", "Y", "Gamma", "Z", "U", "R", "T", "Z"],
            ["Y", "T"],
            ["U", "X"],
            ["S", "R"]]

        return HighSymmetryPoints(kpoints, path)

    def orcf1(self, a, b, c):
        r"""Face-centered orthorhombic (ORCF, oF),

        \dfrac{1}{a^2} > \dfrac{1}{b^2} + \dfrac{1}{c^2}"""

        eta = (1 + a**2 / b**2 + a**2 / c**2) / 4
        zeta = (1 + a**2 / b**2 - a**2 / c**2) / 4

        kpoints = {
            "Gamma": [0, 0, 0],
            "A": [1/2, 1/2 + zeta, zeta],
            "A1": [1/2, 1/2 - zeta, 1 - zeta],
            "L": [1/2, 1/2, 1/2],
            "T": [1, 1/2, 1/2],
            "X": [0, eta, eta],
            "X1": [1, 1 - eta, 1 - eta],
            "Y": [1/2, 0, 1/2],
            "Z": [1/2, 1/2, 0]}

        path = [["Gamma", "Y", "T", "Z", "Gamma", "X", "A1", "Y"],
                ["T", "X1"],
                ["X", "A", "Z"],
                ["L", "Gamma"]]

        return HighSymmetryPoints(kpoints, path)

    def orcf2(self, a, b, c):
        r"""Face-centered orthorhombic (ORCF, oF),

        \dfrac{1}{a^2} < \dfrac{1}{b^2} + \dfrac{1}{c^2}"""

        eta = (1 + a**2 / b**2 + a**2 / c**2) / 4
        delta = (1 + b**2 / a**2 - b**2 / c**2) / 4
        phi = (1 + c**2 / b**2 - c**2 / a**2) / 4

        kpoints = {
            "Gamma": [0, 0, 0],
            "C": [1/2, 1/2 - eta, 1 - eta],
            "C1": [1/2, 1/2 + eta, 1 - eta],
            "D": [1/2 - delta, 1/2, 1 - delta],
            "D1": [1/2 + delta, 1/2, delta],
            "L": [1/2, 1/2, 1/2],
            "H": [1 - phi, 1/2 - phi, 1/2],
            "H1": [phi, 1/2 + phi, 1/2],
            "X": [0, 1/2, 1/2],
            "Y": [1/2, 0, 1/2],
            "Z": [1/2, 1/2, 0]}

        path = [["Gamma", "Y", "C", "D", "X", "Gamma", "Z", "D1", "H", "C"],
                ["C1", "Z"],
                ["X", "H1"],
                ["H", "Y"],
                ["L", "Gamma"]]

        return HighSymmetryPoints(kpoints, path)

    def orcf3(self, a, b, c):
        r"""Face-centered orthorhombic (ORCF, oF),

        \dfrac{1}{a^2} = \dfrac{1}{b^2} + \dfrac{1}{c^2}"""

        eta = (1 + a**2 / b**2 + a**2 / c**2) / 4
        zeta = (1 + a**2 / b**2 - a**2 / c**2) / 4

        kpoints = {
            "Gamma": [0, 0, 0],
            "A": [1/2, 1/2 + zeta, zeta],
            "A1": [1/2, 1/2 - zeta, 1 - zeta],
            "L": [1/2, 1/2, 1/2],
            "T": [1, 1/2, 1/2],
            "X": [0, eta, eta],
            "X1": [1, 1 - eta, 1 - eta],
            "Y": [1/2, 0, 1/2],
            "Z": [1/2, 1/2, 0]}

        path = [["Gamma", "Y", "T", "Z", "Gamma", "X", "A1", "Y"],
                ["X", "A", "Z"],
                ["L", "Gamma"]]

        return HighSymmetryPoints(kpoints, path)

    def orci(self, a, b, c):
        r"""Body-centered orthorhombic (ORCI, oI)"""

        eta = (1 + a**2 / c**2) / 4
        zeta = (1 + b**2 / c**2) / 4
        delta = (b**2 - a**2) / (4 * c**2)
        mu = (a**2 + b**2) / (4 * c**2)

        kpoints = {
            "Gamma": [0, 0, 0],
            "L": [-mu, mu, 1/2 - delta],
            "L1": [mu, -mu, 1/2 + delta],
            "L2": [1/2 - delta, 1/2 + delta, -mu],
            "R": [0, 1/2, 0],
            "S": [1/2, 0, 0],
            "T": [0, 0, 1/2],
            "W": [1/4, 1/4, 1/4],
            "X": [-zeta, zeta, zeta],
            "X1": [zeta, 1 - zeta, -zeta],
            "Y": [eta, -eta, eta],
            "Y1": [1 - eta, eta, -eta],
            "Z": [1/2, 1/2, -1/2]}

        path = [["Gamma", "X", "L", "T", "W", "R", "X1", "Z", "Gamma", "Y", "S", "W"],
                ["L1", "Y"],
                ["Y1", "Z"]]

        return HighSymmetryPoints(kpoints, path)

    def orcc(self, a, b):
        r"""C-centered orthorombic (ORCC, oS)"""

        zeta = (1 + a**2 / b**2) / 4

        kpoints = {
            "Gamma": [0, 0, 0],
            "A": [zeta, zeta, 1/2],
            "A1": [-zeta, 1 - zeta, 1/2],
            "R": [0, 1/2, 1/2],
            "S": [0, 1/2, 0],
            "T": [-1/2, 1/2, 1/2],
            "X": [zeta, zeta, 0],
            "X1": [-zeta, 1 - zeta, 0],
            "Y": [-1/2, 1/2, 0],
            "Z": [0, 0, 1/2]}

        path = [["Gamma", "X", "S", "R", "A", "Z", "Gamma", "Y", "X1", "A1", "T", "Y"],
                ["Z", "T"]]

        return HighSymmetryPoints(kpoints, path)

    def hex(self):
        r"""Hexagonal (HEX, hP)"""

        kpoints = {
            "Gamma": [0, 0, 0],
            "A": [0, 0, 1/2],
            "H": [1/3, 1/3, 1/2],
            "K": [1/3, 1/3, 0],
            "L": [1/2, 0, 1/2],
            "M": [1/2, 0, 0]}

        path = [["Gamma", "M", "K", "Gamma", "A", "L", "H", "A"],
                ["L", "M"],
                ["K", "H"]]

        return HighSymmetryPoints(kpoints, path)

    def rhl1(self, alpha):
        r"""Rhombohedral (RHL, hR),

        \alpha < 90^{\circ}"""

        eta = (1 + 4 * cos(alpha)) / (2 + 4 * cos(alpha))
        nu = 3 / 4 - eta / 2

        kpoints = {
            "Gamma": [0, 0, 0],
            "B": [eta, 1/2, 1 - eta],
            "B1": [1/2, 1 - eta, eta - 1],
            "F": [1/2, 1/2, 0],
            "L": [1/2, 0, 0],
            "L1": [0, 0, -1/2],
            "P": [eta, nu, nu],
            "P1": [1 - nu, 1 - nu, 1 - eta],
            "P2": [nu, nu, eta - 1],
            "Q": [1 - nu, nu, 0],
            "X": [nu, 0, -nu],
            "Z": [1/2, 1/2, 1/2]}

        path = [["Gamma", "L", "B1"],
                ["B", "Z", "Gamma", "X"],
                ["Q", "F", "P1", "Z"],
                ["L", "P"]]

        return HighSymmetryPoints(kpoints, path)

    def rhl2(self, alpha):
        r"""Rhombohedral (RHL, hR),

        \alpha > 90^{\circ}"""

        eta = 1 / (2 * tan(alpha / 2)**2)
        nu = 3 / 4 - eta / 2

        kpoints = {
            "Gamma": [0, 0, 0],
            "F": [1/2, -1/2, 0],
            "L": [1/2, 0, 0],
            "P": [1 - nu, -nu, 1 - nu],
            "P1": [nu, nu - 1, nu - 1],
            "Q": [eta, eta, eta],
            "Q1": [1 - eta, -eta, -eta],
            "Z": [1/2, -1/2, 1/2]}

        path = [["Gamma", "P", "Z", "Q", "Gamma", "F", "P1", "Q1", "L", "Z"]]

        return HighSymmetryPoints(kpoints, path)

    def mcl(self, b, c, alpha):
        r"""Monoclinic (MCL, mP)"""

        eta = (1 - b * cos(alpha) / c) / (2 * sin(alpha)**2)
        nu = 1 / 2 - eta * c * cos(alpha) / b

        kpoints = {
            "Gamma": [0, 0, 0],
            "A": [1/2, 1/2, 0],
            "C": [0, 1/2, 1/2],
            "D": [1/2, 0, 1/2],
            "D1": [1/2, 0, -1/2],
            "E": [1/2, 1/2, 1/2],
            "H": [0, eta, 1 - nu],
            "H1": [0, 1 - eta, nu],
            "H2": [0, eta, -nu],
            "M": [1/2, eta, 1 - nu],
            "M1": [1/2, 1 - eta, nu],
            "M2": [1/2, eta, -nu],
            "X": [0, 1/2, 0],
            "Y": [0, 0, 1/2],
            "Y1": [0, 0, -1/2],
            "Z": [1/2, 0, 0]}

        path = [["Gamma", "Y", "H", "C", "E", "M1", "A", "X", "H1"],
                ["M", "D", "Z"],
                ["Y", "D"]]

        return HighSymmetryPoints(kpoints, path)

    def mclc1(self, a, b, c, alpha):
        r"""C-centered monoclinic (MCLC, mS),

        k_{\gamma} > 90^{\circ}"""

        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha)**2)
        eta = 1 / 2 + 2 * eta * c * cos(alpha) / b
        psi = 3 / 4 - a**2 / (4 * b**2 * sin(alpha)**2)
        phi = psi + (3 / 4 - psi) * b * cos(alpha) / c

        kpoints = {
            "Gamma": [0, 0, 0],
            "N": [1/2, 0, 0],
            "N1": [0, -1/2, 0],
            "F": [1 - zeta, 1 - zeta, 1 - eta],
            "F1": [zeta, zeta, eta],
            "F2": [-zeta, -zeta, 1 - eta],
            "F3": [1 - zeta, -zeta, 1 - eta],
            "I": [phi, 1 - phi, 1/2],
            "I1": [1 - phi, phi - 1, 1/2],
            "L": [1/2, 1/2, 1/2],
            "M": [1/2, 0, 1/2],
            "X": [1 - psi, psi - 1, 0],
            "X1": [psi, 1 - psi, 0],
            "X2": [psi - 1, -psi, 0],
            "Y": [1/2, 1/2, 0],
            "Y1": [-1/2, -1/2, 0],
            "Z": [0, 0, 1/2]}

        path = [["Gamma", "Y", "F", "L", "I"],
                ["I1", "Z", "F1"],
                ["Y", "X1"],
                ["X", "Gamma", "N"],
                ["M", "Gamma"]]

        return HighSymmetryPoints(kpoints, path)

    def mclc2(self, a, b, c, alpha):
        r"""C-centered monoclinic (MCLC, mS),

        k_{\gamma} = 90^{\circ}"""

        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha)**2)
        eta = 1 / 2 + 2 * eta * c * cos(alpha) / b
        psi = 3 / 4 - a**2 / (4 * b**2 * sin(alpha)**2)
        phi = psi + (3 / 4 - psi) * b * cos(alpha) / c

        kpoints = {
            "Gamma": [0, 0, 0],
            "N": [1/2, 0, 0],
            "N1": [0, -1/2, 0],
            "F": [1 - zeta, 1 - zeta, 1 - eta],
            "F1": [zeta, zeta, eta],
            "F2": [-zeta, -zeta, 1 - eta],
            "F3": [1 - zeta, -zeta, 1 - eta],
            "I": [phi, 1 - phi, 1/2],
            "I1": [1 - phi, phi - 1, 1/2],
            "L": [1/2, 1/2, 1/2],
            "M": [1/2, 0, 1/2],
            "X": [1 - psi, psi - 1, 0],
            "X1": [psi, 1 - psi, 0],
            "X2": [psi - 1, -psi, 0],
            "Y": [1/2, 1/2, 0],
            "Y1": [-1/2, -1/2, 0],
            "Z": [0, 0, 1/2]}

        path = [["Gamma", "Y", "F", "L", "I"],
                ["I1", "Z", "F1"],
                ["N", "Gamma", "M"]]

        return HighSymmetryPoints(kpoints, path)

    def mclc3(self, a, b, c, alpha):
        r"""C-centered monoclinic (MCLC, mS),

        k_{\gamma} < 90^{\circ}, 
        \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} < 1"""

        mu = (1 + b**2 / a**2) / 4
        delta = b * c * cos(alpha) / (2 * a**2)
        zeta = mu - 1 / 4 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ^ 2)
        eta = 1 / 2 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta

        kpoints = {
            "Gamma": [0, 0, 0],
            "F": [1 - phi, 1 - phi, 1 - psi],
            "F1": [phi, phi - 1, psi],
            "F2": [1 - phi, -phi, 1 - psi],
            "H": [zeta, zeta, eta],
            "H1": [1 - zeta, -zeta, 1 - eta],
            "H2": [-zeta, -zeta, 1 - eta],
            "I": [1/2, -1/2, 1/2],
            "M": [1/2, 0, 1/2],
            "N": [1/2, 0, 0],
            "N1": [0, -1/2, 0],
            "X": [1/2, -1/2, 0],
            "Y": [mu, mu, delta],
            "Y1": [1 - mu, -mu, -delta],
            "Y2": [-mu, -mu, -delta],
            "Y3": [mu, mu - 1, delta],
            "Z": [0, 0, 1/2]}

        path = [["Gamma", "Y", "F", "H", "Z", "I", "F1"],
                ["H1", "Y1", "X", "Gamma", "N"],
                ["M", "Gamma"]]

        return HighSymmetryPoints(kpoints, path)

    def mclc4(self, a, b, c, alpha):
        r"""C-centered monoclinic (MCLC, mS),

        k_{\gamma} < 90^{\circ}, 
        \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} = 1"""

        mu = (1 + b**2 / a**2) / 4
        delta = b * c * cos(alpha) / (2 * a**2)
        zeta = mu - 1 / 4 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ^ 2)
        eta = 1 / 2 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta

        kpoints = {
            "Gamma": [0, 0, 0],
            "F": [1 - phi, 1 - phi, 1 - psi],
            "F1": [phi, phi - 1, psi],
            "F2": [1 - phi, -phi, 1 - psi],
            "H": [zeta, zeta, eta],
            "H1": [1 - zeta, -zeta, 1 - eta],
            "H2": [-zeta, -zeta, 1 - eta],
            "I": [1/2, -1/2, 1/2],
            "M": [1/2, 0, 1/2],
            "N": [1/2, 0, 0],
            "N1": [0, -1/2, 0],
            "X": [1/2, -1/2, 0],
            "Y": [mu, mu, delta],
            "Y1": [1 - mu, -mu, -delta],
            "Y2": [-mu, -mu, -delta],
            "Y3": [mu, mu - 1, delta],
            "Z": [0, 0, 1/2]}

        path = [["Gamma", "Y", "F", "H", "Z", "I"],
                ["H1", "Y1", "X", "Gamma", "N"],
                ["M", "Gamma"]]

        return HighSymmetryPoints(kpoints, path)

    def mclc5(self, a, b, c, alpha):
        r"""C-centered monoclinic (MCLC, mS),

        k_{\gamma} < 90^{\circ}, 
        \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} > 1"""

        zeta = (b**2 / a**2 + (1 - b * cos(alpha) / c) / sin(alpha)**2) / 4
        mu = zeta / 2 + b**2 / (4 * a**2) - b * c * cos(alpha) / (2 * a**2)
        nu = 2 * mu - zeta
        rho = 1 - zeta * a**2 / b**2
        omega = (4 * nu - 1 - b**2 * sin(alpha)**2 / a**2) * \
            c / (2 * b * cos(alpha))
        eta = 1 / 2 + 2 * zeta * c * cos(alpha) / b
        delta = zeta * c * cos(alpha) / b + omega / 2 - 1 / 4

        kpoints = {
            "Gamma": [0, 0, 0],
            "F": [nu, nu, omega],
            "F1": [1 - nu, 1 - nu, 1 - omega],
            "F2": [nu, nu - 1, omega],
            "H": [zeta, zeta, eta],
            "H1": [1 - zeta, -zeta, 1 - eta],
            "H2": [-zeta, -zeta, 1 - eta],
            "I": [rho, 1 - rho, 1/2],
            "I1": [1 - rho, rho - 1, 1/2],
            "L": [1/2, 1/2, 1/2],
            "M": [1/2, 0, 1/2],
            "N": [1/2, 0, 0],
            "N1": [0, -1/2, 0],
            "X": [1/2, -1/2, 0],
            "Y": [mu, mu, delta],
            "Y1": [1 - mu, -mu, -delta],
            "Y2": [-mu, -mu, -delta],
            "Y3": [mu, mu - 1, delta],
            "Z": [0, 0, 1/2]}

        path = [["Gamma", "Y", "F", "L", "I"],
                ["I1", "Z", "H", "F1"],
                ["H1", "Y1", "X", "Gamma", "N"],
                ["M", "Gamma"]]

        return HighSymmetryPoints(kpoints, path)

    def tri1a(self, a):
        r"""Triclinic (TRI, aP),

        k_{'alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} > 90^{circ},
        k_{gamma} = min(k_{\alpha}, k_{\beta}, k_{\gamma))"""

        kpoints = {
            "Gamma": [0, 0, 0],
            "L": [1/2, 1/2, 0],
            "M": [0, 1/2, 1/2],
            "N": [1/2, 0, 1/2],
            "R": [1/2, 1/2, 1/2],
            "X": [1/2, 0, 0],
            "Y": [0, 1/2, 0],
            "Z": [0, 0, 1/2]}

        path = [["X", "Gamma", "Y"],
                ["L", "Gamma", "Z"],
                ["N", "Gamma", "M"],
                "R", "Gamma"]

        return HighSymmetryPoints(kpoints, path)

    def tri1b(self, a):
        r"""Triclinic (TRI, aP),

        k_{'alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} < 90^{circ},
        k_{gamma} = max(k_{\alpha}, k_{\beta}, k_{\gamma))"""

        kpoints = {
            "Gamma": [0, 0, 0],
            "L": [1/2, -1/2, 0],
            "M": [0, 0, 1/2],
            "N": [-1/2, -1/2, 1/2],
            "R": [0, -1/2, 1/2],
            "X": [0, -1/2, 0],
            "Y": [1/2, 0, 0],
            "Z": [-1/2, 0, 1/2]}

        path = [["X", "Gamma", "Y"],
                ["L", "Gamma", "Z"],
                ["N", "Gamma", "M"],
                "R", "Gamma"]

        return HighSymmetryPoints(kpoints, path)

    def tri2a(self, a):
        r"""Triclinic (TRI, aP),

        k_{'alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} = 90^{circ}"""

        kpoints = {
            "Gamma": [0, 0, 0],
            "L": [1/2, 1/2, 0],
            "M": [0, 1/2, 1/2],
            "N": [1/2, 0, 1/2],
            "R": [1/2, 1/2, 1/2],
            "X": [1/2, 0, 0],
            "Y": [0, 1/2, 0],
            "Z": [0, 0, 1/2]}

        path = [["X", "Gamma", "Y"],
                ["L", "Gamma", "Z"],
                ["N", "Gamma", "M"],
                "R", "Gamma"]

        return HighSymmetryPoints(kpoints, path)

    def tri2b(self, a):
        r"""Triclinic (TRI, aP),

        k_{'alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} = 90^{circ}"""

        kpoints = {
            "Gamma": [0, 0, 0],
            "L": [1/2, -1/2, 0],
            "M": [0, 0, 1/2],
            "N": [-1/2, -1/2, 1/2],
            "R": [0, -1/2, 1/2],
            "X": [0, -1/2, 0],
            "Y": [1/2, 0, 0],
            "Z": [-1/2, 0, 1/2]}

        path = [["X", "Gamma", "Y"],
                ["L", "Gamma", "Z"],
                ["N", "Gamma", "M"],
                "R", "Gamma"]

        return HighSymmetryPoints(kpoints, path)
