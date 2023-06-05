import radtools as rad

l = rad.lattice_example(f"ORCF3")
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
    "orcf3_real.png",
    elev=27,
    azim=36,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=27, azim=36)
