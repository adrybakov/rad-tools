import radtools as rad

l = rad.lattice_example(f"BCT1")
l.plot(
    "primitive",
    label="primitive",
)
l.legend()
l.plot(
    "conventional",
    label="conventional",
    colour="black"
)
l.legend()
# Save an image:
l.savefig(
    "bct1_real.png",
    elev=37,
    azim=72,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=37, azim=72)
