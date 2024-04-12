# RAD-tools - Sandbox (mainly condense matter plotting).
# Copyright (C) 2022-2024  Andrey Rybakov
#
# e-mail: anry@uv.es, web: rad-tools.org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os

import numpy as np

from radtools.crystal.atom import Atom
from radtools.crystal.crystal import Crystal
from radtools.decorate.array import print_2d_array
from radtools.geometry import absolute_to_relative, volume

__all__ = ["load_poscar", "dump_poscar"]


def load_poscar(file_object=None, return_crystal=True, return_comment=False):
    r"""
    Read the crystal structure from the |POSCAR|_ file.

    Parameters
    ----------
    file_object: str of file-like object, optional
        File to be read. If str, then file is opened with the given name.
        Otherwise it has to have ``.readlines()`` method.
        By default it looks for the "POSCAR" file in the current directory.
        Behaviour for ``str``:

        * Tries to open the file with the name given by the ``file_object``.
        * Tries to open the file with the name "POSCAR" in the directory given by the ``file_object``.
    return_crystal: bool, default True
        If True, returns :py:class:`.Crystal` object. Otherwise returns a tuple of
        ``(cell, atoms)``.
    return_comment: bool, default False
        Whether to return the comment from the first line of the file.

    Returns
    -------
    crystal: :py:class:`.Crystal`
        Crystal structure read from the file. If ``return_crystal`` is ``True``.
    cell: (3, 3) :numpy:`ndarray`
        Cell of the crystal structure. If ``return_crystal`` is ``False``.
    atoms: list of :py:class:`.Atom`
        Atoms of the crystal structure. If ``return_crystal`` is ``False``.
        Positions are always relative to the cell.
    comment: str
        Comment from the first line of the file. If ``return_comment`` is ``True``.
    """

    # Open file if needed
    if file_object is None:
        try:
            file_object = open("POSCAR")
        except FileNotFoundError:
            raise FileNotFoundError("POSCAR file not found")
    elif isinstance(file_object, str):
        try:
            file_object = open(file_object)
        except FileNotFoundError:
            try:
                file_object = open(os.path.join(file_object, "POSCAR"))
            except FileNotFoundError:
                raise FileNotFoundError("POSCAR file not found")

    lines = file_object.readlines()

    comment = lines[0].strip()

    # 1 or 3 numbers
    scale_factor = np.array(list(map(float, lines[1].split())))
    cell = np.array(list(map(lambda x: list(map(float, x.split())), lines[2:5])))
    if len(scale_factor) == 1:
        if scale_factor[0] < 0:
            scale_factor = abs(scale_factor[0] / volume(cell))
        cell *= scale_factor
    elif len(scale_factor) == 3 and np.all(scale_factor > 0):
        cell[0] *= scale_factor[0]
        cell[1] *= scale_factor[1]
        cell[2] *= scale_factor[2]
    else:
        raise ValueError(
            "Scale factor has to be a single positive ot negative number or "
            + f"a list of 3 positive numbers, got: {scale_factor}"
        )
    # Read species name and numbers
    species = lines[5].split()
    GOT_SPECIES_NAMES = False
    index = 6
    for i in species:
        try:
            int(i)
        except ValueError:
            GOT_SPECIES_NAMES = True
    if GOT_SPECIES_NAMES:
        species_names = species
        species = lines[6].split()
        index = 7
    else:
        species_names = None
    ions_per_species = list(map(int, species))

    # Skip selective dynamics
    if lines[index][0] in ["S", "s"]:
        index += 1

    # Get mode
    CARTESIAN = False
    if lines[index][0] in ["C", "c", "K", "k"]:
        CARTESIAN = True
    index += 1
    atoms = []
    for i in range(len(species)):
        for j in range(ions_per_species[i]):
            coordinates = np.array(list(map(float, lines[index].split()[:3])))
            index += 1
            if CARTESIAN:
                # Both cases (1 or 3 numbers) are covered
                coordinates *= scale_factor
                coordinates = absolute_to_relative(cell, coordinates)
            if species_names is None:
                atoms.append(Atom(f"X{i+1}", coordinates))
            else:
                atoms.append(Atom(species_names[i], coordinates))

    if return_crystal:
        if return_comment:
            return Crystal(cell=cell, atoms=atoms), comment
        return Crystal(cell=cell, atoms=atoms)
    if return_comment:
        return cell, atoms, comment
    return cell, atoms


def dump_poscar(
    crystal_like,
    file_object="POSCAR",
    comment: str = None,
    decimals=8,
    mode: str = "Direct",
):
    r"""
    Write :py:class:`.Crystal`-like object  to the |POSCAR|_ file.

    Parameters
    ----------
    crystal_like: any
        Object to be written. It has to have ``.cell`` and ``.atoms`` attributes.
        ``cell`` has to be a 3x3 array-like object, ``atoms`` has to be a list of
        :py:class:`.Atom` objects.
    file_object: str of file-like object, optional
        File to be written. If str, then file is opened with the given name.
        Otherwise it has to have ``.write()`` method.
    comment: str, optional
        Comment to be written in the first line of the file. Has to be a single line.
        All new lines are replaced with spaces.
    decimals: int, default 8
        Number of decimals to be written.
    mode: str, default "Direct"
        Mode of the coordinates to be written. Can be "Direct" or "Cartesian".
    """

    # Prepare comment
    if comment is None:
        try:
            comment = crystal_like.name
        except AttributeError:
            comment = "Data from RAD-tools"
    else:
        comment = comment.replace("\n", " ")

    # Check mode
    if mode not in ("Direct", "Cartesian"):
        raise ValueError(f'mode has to be "Direct" or "Cartesian", given: {mode}')

    # Prepare atoms
    if mode == "Direct":
        atoms = [(atom.type, atom.position) for atom in crystal_like.atoms]
    else:
        atoms = [
            (atom.type, atom.position @ crystal_like.cell)
            for atom in crystal_like.atoms
        ]
    # Sort atoms by type
    atoms = sorted(atoms, key=lambda x: x[0])

    # Prepare atom species and coordinates
    atom_species = {}
    atom_coordinates = []
    for atom in atoms:
        if atom[0] not in atom_species:
            atom_species[atom[0]] = 1
        else:
            atom_species[atom[0]] += 1
        atom_coordinates.append(atom[1])

    # Open file if needed
    if isinstance(file_object, str):
        file_object = open(file_object, "w", encoding="utf-8")

    # Write
    file_object.write(comment + "\n")
    file_object.write("1.0\n")
    file_object.write(
        print_2d_array(
            crystal_like.cell, fmt=f".{decimals}f", print_result=False, borders=False
        )
        + "\n"
    )

    for species in atom_species:
        file_object.write(f"{species} ")
    file_object.write("\n")
    for species in atom_species:
        file_object.write(f"{atom_species[species]} ")
    file_object.write("\n")

    file_object.write(mode + "\n")
    file_object.write(
        print_2d_array(
            atom_coordinates, fmt=f".{decimals}f", print_result=False, borders=False
        )
        + "\n"
    )


if __name__ == "__main__":
    from radtools import HEX, Crystal

    c = Crystal(HEX(1, 2))
    c.add_atom("Cr2", position=[0.5, 0, 0])
    c.add_atom("C2", position=[0, 0, 0])
    c.add_atom("C1", position=[0, 0, 0.5])
    c.add_atom("C3", position=[0, 0.5, 0])
    c.add_atom("Cr1", position=[0, 0.5, 0.5])
    # dump_poscar(c)
    # dump_poscar(c, "test-cart.vasp", mode="Cartesian")
    crystal = load_poscar()
    print_2d_array(crystal.cell)
    for atom in crystal.atoms:
        print(atom.name, atom.position)
    crystal = load_poscar("test-cart.vasp")
    print_2d_array(crystal.cell)
    for atom in crystal.atoms:
        print(atom.name, atom.position)
