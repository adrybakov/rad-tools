__all__ = [
    "PREDEFINED_NOTATIONS",
    "TXT_FLAGS",
]

PREDEFINED_NOTATIONS = {
    "standard": (True, False, False, False, True),
    "tb2j": (True, True, False, False, True),
    "spinw": (True, False, False, False, False),
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
    "separate": "Data are in the file:"
}
