import radtools as rad

l = rad.lattice_example(f"MCLC1")
l.plot("brillouin-kpath")
# Save an image:
l.savefig(
    "mclc1_brillouin.png",
    elev=21,
    azim=53,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=21, azim=53)
