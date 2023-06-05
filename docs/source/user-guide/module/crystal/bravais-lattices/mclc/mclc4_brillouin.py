import radtools as rad

l = rad.lattice_example(f"MCLC4")
l.plot("brillouin-kpath")
# Save an image:
l.savefig(
    "mclc4_brillouin.png",
    elev=19,
    azim=45,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=19, azim=45)
