import radtools as rad

l = rad.lattice_example(f"HEX")
l.plot("wigner-seitz")
# Save an image:
l.savefig(
    "hex_wigner-seitz.png",
    elev=32,
    azim=10,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=32, azim=10)
