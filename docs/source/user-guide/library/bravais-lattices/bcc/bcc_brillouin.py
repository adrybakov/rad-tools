import radtools as rad

l = rad.lattice_example("BCC")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("bcc_brillouin.png", elev=11, azim=25, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=11, azim=25)
