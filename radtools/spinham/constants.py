__all__ = [
    "PREDEFINED_NOTATIONS",
    "TXT_FLAGS",
]

PREDEFINED_NOTATIONS = {
    "standard": (True, False, -1),
    "tb2j": (True, True, -1),
    "spinw": (True, False, 1),
    "vampire": (True, True, -0.5),
}

TXT_FLAGS = {
    "cell": "Primitive unit cell:",
    "atoms": "Atoms:",
    "magmoms": "Magnetic moments of magnetic atoms:",
    "spins": "Spin of magnetic atoms:",
    "spinham": "Spin Hamiltonian:",
    "iso": "Isotropic exchange:",
    "matrix": "Full exchange matrix:",
    "dmi": "DMI vector:",
    "dmis": "DMI vectors:",
    "dmi_module": "DMI module:",
    "dmi_relative": "|DMI| / |J_iso|:",
    "aniso": "Symmetric anisotropic exchange:",
    "spinmodel": "Spin model:",
    "kpath": "K-path:",
    "kpoints": "Points (in relative coordinates):",
    "klabels": "Labels:",
    "dispersion": "Magnon dispersion:",
    "separate": "Data are in the file:",
    "notation": "Notation of Hamiltonian:",
    "bonds": "Bonds:",
}
