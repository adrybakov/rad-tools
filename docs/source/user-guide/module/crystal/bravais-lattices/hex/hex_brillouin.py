import radtools as rad

l = rad.lattice_example(f"HEX")
l.plot("brillouin-kpath")
# Save an image:
l.savefig(
    "hex_brillouin.png",
    elev=19,
    azim=20,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=19, azim=20)
