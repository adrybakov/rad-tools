import radtools as rad

l = rad.lattice_example(f"MCLC4")
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
    "mclc4_real.png",
    elev=25,
    azim=70,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=25, azim=70)
