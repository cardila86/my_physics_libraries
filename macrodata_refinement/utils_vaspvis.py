'''
This file is based in the vaspvis package.
'''

import os
import numpy as np
import matplotlib.pyplot as plt

from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.io.vasp.outputs import BSVasprun, Eigenval
from pymatgen.io.vasp.inputs import Kpoints, Poscar, Incar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.core.periodic_table import Element
# from vaspvis.unfold import unfold, make_kpath, removeDuplicateKpoints
from pymatgen.core.periodic_table import Element
from .utils_pyprocar import UtilsProcar
from .utils_pyprocar import ProcarParser
# from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
from copy import deepcopy




class data_vaspvis:
    def __init__(self,
        path_read,
        projected=False,
        spin="up",
        kpath=None,
        n=None,
        M=None,
        high_symm_points=None,
        # # bandgap=False,
        # # printbg=True,
        shift_efermi=0,
        interpolate=True,
        new_n=200,
        custom_kpath=None,
        stretch_factor=1.0,
        path_efermi=None
        ):
        self.path_read = path_read
        self.projected = projected
        self.spin = spin
        self.spin_dict = {"up": Spin.up, "down": Spin.down}
        self.kpath = kpath
        self.n = n
        self.M = M
        self.high_symm_points = high_symm_points
        self.interpolate = interpolate
        self.new_n = new_n
        self.custom_kpath = custom_kpath
        self.stretch_factor = stretch_factor
        self.unfold = False

        self.eigenval = Eigenval(os.path.join(self.path_read, "EIGENVAL"))
        self.kpoints_file = Kpoints.from_file(os.path.join(self.path_read, "KPOINTS"))
        self.wavecar = os.path.join(self.path_read, "WAVECAR")
        self.incar = Incar.from_file(os.path.join(self.path_read, "INCAR"))
        self.poscar = Poscar.from_file(
            os.path.join(self.path_read, "POSCAR"),
            check_for_POTCAR=False,
            read_velocities=False,
        )
        self.forbitals = self._check_f_orb()

        self.efermi = shift_efermi
        if path_efermi is not None:
            outcar_path = os.path.join(path_efermi, "OUTCAR")
            efermi_output = os.popen(f'grep E-fermi {outcar_path}').read().strip()
            if not efermi_output:
                raise ValueError(f"No E-fermi value found in {outcar_path}")
            efermi_value = efermi_output.split()[2]
            self.efermi += float(efermi_value)

        if "LSORBIT" in self.incar:
            if self.incar["LSORBIT"]:
                self.lsorbit = True
            else:
                self.lsorbit = False
        else:
            self.lsorbit = False

        if "ISPIN" in self.incar:
            if self.incar["ISPIN"] == 2:
                self.ispin = True
            else:
                self.ispin = False
        else:
            self.ispin = False

        if "LHFCALC" in self.incar:
            if self.incar["LHFCALC"]:
                self.hse = True
            else:
                self.hse = False
        else:
            self.hse = False

        self.pre_loaded_bands = os.path.isfile(
            os.path.join(self.path_read, "eigenvalues.npy")
        )
        self.eigenvalues, self.kpoints = self._load_bands()

        if self.stretch_factor != 1.0:
            self.eigenvalues *= self.stretch_factor

        self.color_dict = {
            0: "#FF0000",
            1: "#0000FF",
            2: "#008000",
            3: "#800080",
            4: "#E09200",
            5: "#FF5C77",
            6: "#778392",
            7: "#07C589",
            8: "#40BAF2",
            9: "#FF0000",
            10: "#0000FF",
            11: "#008000",
            12: "#800080",
            13: "#E09200",
            14: "#FF5C77",
            15: "#778392",
        }
        self.orbital_labels = {
            0: "s",
            1: "p_{y}",
            2: "p_{z}",
            3: "p_{x}",
            4: "d_{xy}",
            5: "d_{yz}",
            6: "d_{z^{2}}",
            7: "d_{xz}",
            8: "d_{x^{2}-y^{2}}",
            9: "f_{y^{3}x^{2}}",
            10: "f_{xyz}",
            11: "f_{yz^{2}}",
            12: "f_{z^{3}}",
            13: "f_{xz^{2}}",
            14: "f_{zx^{3}}",
            15: "f_{x^{3}}",
        }
        self.spd_relations = {
            "s": 0,
            "p": 1,
            "d": 2,
            "f": 3,
        }

        if self.custom_kpath is not None:
            (
                self.custom_kpath_inds,
                self.custom_kpath_flip,
            ) = self._get_custom_kpath()
        #  else:
        #  self.custom_kpath_inds, self.custom_kpath_flip = None, None

        if projected:
            self.pre_loaded_projections = os.path.isfile(
                os.path.join(self.path_read, "projected_eigenvalues.npy")
            )
            self.projected_eigenvalues = self._load_projected_bands()
        # non collinear spin polarization
        if self.spin in ["x", "y", "z"]:
            self.soc_axis = self.spin
            self.spin = "up"
            if self.lsorbit:
                self.pre_loaded_spin_projections = os.path.isfile(
                    os.path.join(self.path_read, "spin_projections.npy")
                )
                self.spin_projections = self._load_soc_spin_projection()
        else:
            self.soc_axis = None

    def _get_custom_kpath(self):
        flip = (-np.sign(self.custom_kpath) + 1).astype(bool)
        inds = (np.abs(self.custom_kpath) - 1).astype(int)

        return inds, flip
    
    def _check_f_orb(self):
        f_elements = [
            "La",
            "Ac",
            "Ce",
            "Tb",
            "Th",
            "Pr",
            "Dy",
            "Pa",
            "Nd",
            "Ho",
            "U",
            "Pm",
            "Er",
            "Np",
            "Sm",
            "Tm",
            "Pu",
            "Eu",
            "Yb",
            "Am",
            "Gd",
            "Lu",
        ]
        f = False
        for element in self.poscar.site_symbols:
            if element in f_elements:
                f = True

        return f

    def _load_bands(self):
        """
        This function is used to load eigenvalues from the vasprun.xml
        file and into a dictionary which is in the form of
        band index --> eigenvalues

        Returns:
            bands_dict (dict[str][np.ndarray]): Dictionary which contains
                the eigenvalues for each band
        """
        if self.spin == "up":
            spin = 0
        if self.spin == "down":
            spin = 1

        if self.pre_loaded_bands:
            with open(
                os.path.join(self.path_read, "eigenvalues.npy"), "rb"
            ) as eigenvals:
                band_data = np.load(eigenvals)

            if self.ispin and not self.lsorbit:
                eigenvalues = band_data[:, :, [0, 2]]
                kpoints = band_data[0, :, 4:]
            else:
                eigenvalues = band_data[:, :, 0]
                kpoints = band_data[0, :, 2:]
        else:
            if len(self.eigenval.eigenvalues.keys()) > 1:
                eigenvalues_up = np.transpose(
                    self.eigenval.eigenvalues[Spin.up], axes=(1, 0, 2)
                )
                eigenvalues_down = np.transpose(
                    self.eigenval.eigenvalues[Spin.down], axes=(1, 0, 2)
                )
                eigenvalues_up[:, :, 0] = eigenvalues_up[:, :, 0] - self.efermi
                eigenvalues_down[:, :, 0] = (
                    eigenvalues_down[:, :, 0] - self.efermi
                )
                eigenvalues = np.concatenate(
                    [eigenvalues_up, eigenvalues_down], axis=2
                )
            else:
                eigenvalues = np.transpose(
                    self.eigenval.eigenvalues[Spin.up], axes=(1, 0, 2)
                )
                eigenvalues[:, :, 0] = eigenvalues[:, :, 0] - self.efermi

            kpoints = np.array(self.eigenval.kpoints)

            if self.hse:
                kpoint_weights = np.array(self.eigenval.kpoints_weights)
                zero_weight = np.where(kpoint_weights == 0)[0]
                eigenvalues = eigenvalues[:, zero_weight]
                kpoints = kpoints[zero_weight]

            band_data = np.append(
                eigenvalues,
                np.tile(kpoints, (eigenvalues.shape[0], 1, 1)),
                axis=2,
            )

            np.save(os.path.join(self.path_read, "eigenvalues.npy"), band_data)

            if len(self.eigenval.eigenvalues.keys()) > 1:
                eigenvalues = eigenvalues[:, :, [0, 2]]
            else:
                eigenvalues = eigenvalues[:, :, 0]

        if len(self.eigenval.eigenvalues.keys()) > 1:
            eigenvalues = eigenvalues[:, :, spin]

        return eigenvalues, kpoints

    def _load_projected_bands(self):
        """
        This function loads the project weights of the orbitals in each band
        from vasprun.xml into a dictionary of the form:
        band index --> atom index --> weights of orbitals

        Returns:
            projected_dict (dict([str][int][pd.DataFrame])): Dictionary containing the projected weights of all orbitals on each atom for each band.
        """

        if self.lsorbit:
            if self.soc_axis is None:
                spin = 0
            elif self.soc_axis == "x":
                spin = 1
            elif self.soc_axis == "y":
                spin = 2
            elif self.soc_axis == "z":
                spin = 3
        else:
            if self.spin == "up":
                spin = 0
            elif self.spin == "down":
                spin = 1

        if not os.path.isfile(os.path.join(self.path_read, "PROCAR_repaired")):
            UtilsProcar().ProcarRepair(
                os.path.join(self.path_read, "PROCAR"),
                os.path.join(self.path_read, "PROCAR_repaired"),
            )

        if self.pre_loaded_projections:
            with open(
                os.path.join(self.path_read, "projected_eigenvalues.npy"), "rb"
            ) as projected_eigenvals:
                projected_eigenvalues = np.load(projected_eigenvals)
        else:
            parser = ProcarParser()
            parser.readFile(os.path.join(self.path_read, "PROCAR_repaired"))
            if (
                self.ispin
                and not self.lsorbit
                and np.sum(self.poscar.natoms) == 1
            ):
                shape = int(parser.spd.shape[1] / 2)
                projected_eigenvalues_up = np.transpose(
                    parser.spd[:, :shape, 0, :, 1:-1], axes=(1, 0, 2, 3)
                )
                projected_eigenvalues_down = np.transpose(
                    parser.spd[:, shape:, 0, :, 1:-1], axes=(1, 0, 2, 3)
                )
                projected_eigenvalues = np.concatenate(
                    [
                        projected_eigenvalues_up[:, :, :, :, np.newaxis],
                        projected_eigenvalues_down[:, :, :, :, np.newaxis],
                    ],
                    axis=4,
                )
                projected_eigenvalues = np.transpose(
                    projected_eigenvalues, axes=(0, 1, 4, 2, 3)
                )
            elif (
                self.ispin
                and not self.lsorbit
                and np.sum(self.poscar.natoms) != 1
            ):
                shape = int(parser.spd.shape[1] / 2)
                projected_eigenvalues_up = np.transpose(
                    parser.spd[:, :shape, 0, :-1, 1:-1], axes=(1, 0, 2, 3)
                )
                projected_eigenvalues_down = np.transpose(
                    parser.spd[:, shape:, 0, :-1, 1:-1], axes=(1, 0, 2, 3)
                )
                projected_eigenvalues = np.concatenate(
                    [
                        projected_eigenvalues_up[:, :, :, :, np.newaxis],
                        projected_eigenvalues_down[:, :, :, :, np.newaxis],
                    ],
                    axis=4,
                )
                projected_eigenvalues = np.transpose(
                    projected_eigenvalues, axes=(0, 1, 4, 2, 3)
                )
            else:
                if np.sum(self.poscar.natoms) == 1:
                    projected_eigenvalues = np.transpose(
                        parser.spd[:, :, :, :, 1:-1], axes=(1, 0, 2, 3, 4)
                    )
                else:
                    projected_eigenvalues = np.transpose(
                        parser.spd[:, :, :, :-1, 1:-1], axes=(1, 0, 2, 3, 4)
                    )

            np.save(
                os.path.join(self.path_read, "projected_eigenvalues.npy"),
                projected_eigenvalues,
            )

        projected_eigenvalues = projected_eigenvalues[:, :, spin, :, :]

        if self.lsorbit and self.soc_axis is not None:
            separated_projections = np.zeros(
                projected_eigenvalues.shape + (2,)
            )
            separated_projections[
                projected_eigenvalues > 0, 0
            ] = projected_eigenvalues[projected_eigenvalues > 0]
            separated_projections[
                projected_eigenvalues < 0, 1
            ] = -projected_eigenvalues[projected_eigenvalues < 0]

            if self.spin == "up":
                soc_spin = 0
            elif self.spin == "down":
                soc_spin = 1

            projected_eigenvalues = separated_projections[..., soc_spin]

        if self.hse:
            kpoint_weights = np.array(self.eigenval.kpoints_weights)
            zero_weight = np.where(kpoint_weights == 0)[0]
            projected_eigenvalues = projected_eigenvalues[:, zero_weight]

        projected_eigenvalues = np.square(projected_eigenvalues)

        return projected_eigenvalues

    def _load_soc_spin_projection(self):
        """
        This function loads the project weights of the orbitals in each band
        from vasprun.xml into a dictionary of the form:
        band index --> atom index --> weights of orbitals

        Returns:
            projected_dict (dict([str][int][pd.DataFrame])): Dictionary containing the projected weights of all orbitals on each atom for each band.
        """

        if not self.lsorbit:
            raise BaseException(
                f"You selected soc_axis='{self.soc_axis}' for a non-soc axis calculation, please set soc_axis=None"
            )
        if self.lsorbit and self.soc_axis == "x":
            spin = 1
        if self.lsorbit and self.soc_axis == "y":
            spin = 2
        if self.lsorbit and self.soc_axis == "z":
            spin = 3

        if not os.path.isfile(os.path.join(self.path_read, "PROCAR_repaired")):
            UtilsProcar().ProcarRepair(
                os.path.join(self.path_read, "PROCAR"),
                os.path.join(self.path_read, "PROCAR_repaired"),
            )

        if self.pre_loaded_spin_projections:
            with open(
                os.path.join(self.path_read, "spin_projections.npy"), "rb"
            ) as spin_projs:
                spin_projections = np.load(spin_projs)
        else:
            parser = ProcarParser()
            parser.readFile(os.path.join(self.path_read, "PROCAR_repaired"))
            spin_projections = np.transpose(
                parser.spd[:, :, :, -1, -1], axes=(1, 0, 2)
            )

            np.save(
                os.path.join(self.path_read, "spin_projections.npy"),
                spin_projections,
            )

        spin_projections = spin_projections[:, :, spin]

        if self.hse:
            kpoint_weights = np.array(self.eigenval.kpoints_weights)
            zero_weight = np.where(kpoint_weights == 0)[0]
            spin_projections = spin_projections[:, zero_weight]

        separated_projections = np.zeros(
            (spin_projections.shape[0], spin_projections.shape[1], 2)
        )
        separated_projections[spin_projections > 0, 0] = spin_projections[
            spin_projections > 0
        ]
        separated_projections[spin_projections < 0, 1] = -spin_projections[
            spin_projections < 0
        ]

        separated_projections = (
            separated_projections / separated_projections.max()
        )

        if self.spin == "up":
            separated_projections = separated_projections[:, :, 0]
        elif self.spin == "down":
            separated_projections = separated_projections[:, :, 1]
        else:
            raise BaseException(
                "The soc_axis feature does not work with spin='both'"
            )

        return separated_projections

    def _sum_spd(self, spd):
        """
        This function sums the weights of the s, p, and d orbitals for each atom
        and creates a dictionary of the form:
        band index --> s,p,d orbital weights

        Returns:
            spd_dict (dict([str][pd.DataFrame])): Dictionary that contains the summed weights for the s, p, and d orbitals for each band
        """

        if not self.forbitals:
            spd_indices = [
                np.array([False for _ in range(9)]) for i in range(3)
            ]
            spd_indices[0][0] = True
            spd_indices[1][1:4] = True
            spd_indices[2][4:] = True
        else:
            spd_indices = [
                np.array([False for _ in range(16)]) for i in range(4)
            ]
            spd_indices[0][0] = True
            spd_indices[1][1:4] = True
            spd_indices[2][4:9] = True
            spd_indices[3][9:] = True

        orbital_contributions = np.sum(self.projected_eigenvalues, axis=2)

        spd_contributions = np.transpose(
            np.array(
                [
                    np.sum(orbital_contributions[:, :, ind], axis=2)
                    for ind in spd_indices
                ]
            ),
            axes=[1, 2, 0],
        )

        #  norm_term = np.sum(spd_contributions, axis=2)[:,:,np.newaxis]
        #  spd_contributions = np.divide(spd_contributions, norm_term, out=np.zeros_like(spd_contributions), where=norm_term!=0)

        spd_contributions = spd_contributions[
            :, :, [self.spd_relations[orb] for orb in spd]
        ]

        return spd_contributions

    def _sum_orbitals(self, orbitals):
        """
        This function finds the weights of desired orbitals for all atoms and
            returns a dictionary of the form:
            band index --> orbital index

        Parameters:
            orbitals (list): List of desired orbitals.
                0 = s
                1 = py
                2 = pz
                3 = px
                4 = dxy
                5 = dyz
                6 = dz2
                7 = dxz
                8 = dx2-y2
                9 = fy3x2
                10 = fxyz
                11 = fyz2
                12 = fz3
                13 = fxz2
                14 = fzx3
                15 = fx3

        Returns:
            orbital_dict (dict[str][pd.DataFrame]): Dictionary that contains the projected weights of the selected orbitals.
        """
        orbital_contributions = self.projected_eigenvalues.sum(axis=2)
        #  norm_term =  np.sum(orbital_contributions, axis=2)[:,:,np.newaxis]
        #  orbital_contributions = np.divide(orbital_contributions, norm_term, out=np.zeros_like(orbital_contributions), where=norm_term!=0)
        orbital_contributions = orbital_contributions[:, :, [orbitals]]

        return orbital_contributions

    def _sum_atoms(self, atoms, spd=False):
        """
        This function finds the weights of desired atoms for all orbitals and
            returns a dictionary of the form:
            band index --> atom index

        Parameters:
            atoms (list): List of desired atoms where atom 0 is the first atom in
                the POSCAR file.

        Returns:
            atom_dict (dict[str][pd.DataFrame]): Dictionary that contains the projected
                weights of the selected atoms.
        """

        if spd:
            if not self.forbitals:
                spd_indices = [
                    np.array([False for _ in range(9)]) for i in range(3)
                ]
                spd_indices[0][0] = True
                spd_indices[1][1:4] = True
                spd_indices[2][4:] = True
            else:
                spd_indices = [
                    np.array([False for _ in range(16)]) for i in range(4)
                ]
                spd_indices[0][0] = True
                spd_indices[1][1:4] = True
                spd_indices[2][4:9] = True
                spd_indices[3][9:] = True

            atoms_spd = np.transpose(
                np.array(
                    [
                        np.sum(
                            self.projected_eigenvalues[:, :, :, ind], axis=3
                        )
                        for ind in spd_indices
                    ]
                ),
                axes=(1, 2, 3, 0),
            )

            #  atoms_spd = atoms_spd[:,:,[atoms], :]

            #  norm_term = np.sum(atoms_spd_to_norm, axis=(2,3))[:,:, np.newaxis]
            #  atoms_spd = np.divide(atoms_spd, norm_term, out=np.zeros_like(atoms_spd), where=norm_term!=0)

            return atoms_spd
        else:
            atoms_array = self.projected_eigenvalues.sum(axis=3)
            #  norm_term = np.sum(atoms_array, axis=2)[:,:,np.newaxis]
            #  atoms_array = np.divide(atoms_array, norm_term, out=np.zeros_like(atoms_array), where=norm_term!=0)
            atoms_array = atoms_array[:, :, [atoms]]

            return atoms_array

    def _sum_elements(
        self, elements, orbitals=False, spd=False, spd_options=None
    ):
        """
        This function sums the weights of the orbitals of specific elements within the
        calculated structure and returns a dictionary of the form:
        band index --> element label --> orbital weights for orbitals = True
        band index --> element label for orbitals = False
        This is useful for structures with many elements because manually entering indicies is
        not practical for large structures.

        Parameters:
            elements (list): List of element symbols to sum the weights of.
            orbitals (bool): Determines whether or not to inclue orbitals or not
                (True = keep orbitals, False = sum orbitals together )
            spd (bool): Determines whether or not to sum the s, p, and d orbitals


        Returns:
            element_dict (dict([str][str][pd.DataFrame])): Dictionary that contains the summed weights for each orbital for a given element in the structure.
        """

        poscar = self.poscar
        natoms = poscar.natoms
        symbols = poscar.site_symbols
        projected_eigenvalues = self.projected_eigenvalues

        element_list = np.hstack(
            [
                [symbols[i] for j in range(natoms[i])]
                for i in range(len(symbols))
            ]
        )

        element_indices = [
            np.where(np.isin(element_list, element))[0] for element in elements
        ]

        element_orbitals = np.transpose(
            np.array(
                [
                    np.sum(projected_eigenvalues[:, :, ind, :], axis=2)
                    for ind in element_indices
                ]
            ),
            axes=(1, 2, 0, 3),
        )

        if orbitals:
            return element_orbitals
        elif spd:
            if not self.forbitals:
                spd_indices = [
                    np.array([False for _ in range(9)]) for i in range(3)
                ]
                spd_indices[0][0] = True
                spd_indices[1][1:4] = True
                spd_indices[2][4:] = True
            else:
                spd_indices = [
                    np.array([False for _ in range(16)]) for i in range(4)
                ]
                spd_indices[0][0] = True
                spd_indices[1][1:4] = True
                spd_indices[2][4:9] = True
                spd_indices[3][9:] = True

            element_spd = np.transpose(
                np.array(
                    [
                        np.sum(element_orbitals[:, :, :, ind], axis=3)
                        for ind in spd_indices
                    ]
                ),
                axes=(1, 2, 3, 0),
            )

            #  norm_term = np.sum(element_spd, axis=(2,3))[:,:,np.newaxis, np.newaxis]
            #  element_spd = np.divide(element_spd, norm_term, out=np.zeros_like(element_spd), where=norm_term!=0)

            return element_spd
        else:
            element_array = np.sum(element_orbitals, axis=3)
            #  norm_term = np.sum(element_array, axis=2)[:,:,np.newaxis]
            #  element_array = np.divide(element_array, norm_term, out=np.zeros_like(element_array), where=norm_term!=0)

            return element_array

    def _get_k_distance(self):
        slices = self._get_slices(unfold=self.unfold, hse=self.hse)
        kdists = []

        if self.custom_kpath is not None:
            #  if self.custom_kpath is None:
            index = self.custom_kpath_inds
        else:
            index = range(len(slices))

        for j, i in enumerate(index):
            inv_cell = deepcopy(self.poscar.structure.lattice.inv_matrix)
            inv_cell_norms = np.linalg.norm(inv_cell, axis=1)
            inv_cell /= inv_cell_norms.min()

            # If you want to be able to compare only identical relative cell lengths
            kpt_c = np.dot(self.kpoints[slices[i]], inv_cell.T)

            # If you want to be able to compare any cell length. Maybe straining an orthorhombic cell or something like that
            # This will mess up relative distances though
            # kpt_c = self.kpoints[slices[i]]
            kdist = np.r_[
                0, np.cumsum(np.linalg.norm(np.diff(kpt_c, axis=0), axis=1))
            ]
            if j == 0:
                kdists.append(kdist)
            else:
                kdists.append(kdist + kdists[-1][-1])

        # kdists = np.array(kdists)

        return kdists

    def _get_kticks(self, ax, wave_vectors, vlinecolor):
        """
        This function extracts the kpoint labels and index locations for a regular
        band structure calculation (non unfolded).

        Parameters:
            ax (matplotlib.pyplot.axis): Axis to append the tick labels
        """

        high_sym_points = self.kpoints_file.kpts

        segements = []
        for i in range(0, len(high_sym_points) - 1):
            if not i % 2:
                segements.append([i, i + 1])

        if self.custom_kpath is not None:
            high_sym_points_inds = []
            for i, b in zip(self.custom_kpath_inds, self.custom_kpath_flip):
                if b:
                    seg = list(reversed(segements[i]))
                else:
                    seg = segements[i]

                high_sym_points_inds.extend(seg)
        else:
            high_sym_points_inds = list(range(len(high_sym_points)))

        num_kpts = self.kpoints_file.num_kpts
        kpts_labels = np.array(
            [
                f"${k}$" if k != "G" else "$\\Gamma$"
                for k in self.kpoints_file.labels
            ]
        )
        all_kpoints = self.kpoints

        group_index = []
        for i, j in enumerate(high_sym_points_inds):
            if i == 0:
                group_index.append([j])
            if i % 2 and not i == len(high_sym_points_inds) - 1:
                group_index.append([j, high_sym_points_inds[i + 1]])
            if i == len(high_sym_points_inds) - 1:
                group_index.append([j])

        labels = []
        index = []

        for i in group_index:
            if len(i) == 1:
                labels.append(kpts_labels[i[0]])
                index.append(i[0])
            else:
                if kpts_labels[i[0]] == kpts_labels[i[1]]:
                    labels.append(kpts_labels[i[0]])
                    index.append(i[0])
                else:
                    merged_label = "|".join(
                        [
                            kpts_labels[i[0]],
                            kpts_labels[i[1]],
                        ]
                    ).replace("$|$", "|")
                    labels.append(merged_label)
                    index.append(i[0])

        kpoints_index = [0] + [
            (i + 1) * num_kpts - 1
            for i in range(int((len(high_sym_points_inds) + 1) / 2))
        ]

        for k in kpoints_index:
            ax.axvline(
                x=wave_vectors[k], color=vlinecolor, alpha=0.7, linewidth=0.5
            )

        ax.set_xticks([wave_vectors[k] for k in kpoints_index])
        ax.set_xticklabels(labels)

    def _get_slices(self, unfold=False, hse=False):
        if not unfold and not hse:
            high_sym_points = self.kpoints_file.kpts
            all_kpoints = self.kpoints
            num_kpts = self.kpoints_file.num_kpts
            num_slices = int(len(high_sym_points) / 2)
            slices = [
                slice(i * num_kpts, (i + 1) * num_kpts, None)
                for i in range(num_slices)
            ]

        if hse and not unfold:
            structure = self.poscar.structure
            kpath_obj = HighSymmKpath(structure)
            kpath_coords = np.array(list(kpath_obj._kpath["kpoints"].values()))
            index = np.where(
                np.isclose(
                    self.kpoints[:, None],
                    kpath_coords,
                )
                .all(-1)
                .any(-1)
                == True
            )[0]

            segements = []
            for i in range(0, len(index) - 1):
                if not i % 2:
                    segements.append([index[i], index[i + 1]])

            # print(segements)

            num_kpts = int(len(self.kpoints) / (len(index) / 2))
            slices = [
                slice(i * num_kpts, (i + 1) * num_kpts, None)
                for i in range(int(len(index) / 2))
            ]
            # print(slices)
            slices = [slice(i[0], i[1] + 1, None) for i in segements]
            # print(slices)

        if unfold and not hse:
            n = int(len(self.kpoints) / len(self.kpath))
            slices = [
                slice(i * n, (i + 1) * n, None)
                for i in range(int(len(self.kpath)))
            ]

        return slices

    def _get_interpolated_data_segment(
        self, wave_vectors, data, crop_zero=False, kind="cubic"
    ):
        data_shape = data.shape

        if len(data_shape) == 1:
            fs = interp1d(wave_vectors, data, kind=kind, axis=0)
        else:
            fs = interp1d(wave_vectors, data, kind=kind, axis=1)

        new_wave_vectors = np.linspace(
            wave_vectors.min(), wave_vectors.max(), self.new_n
        )
        data = fs(new_wave_vectors)

        if crop_zero:
            data[np.where(data < 0)] = 0

        return new_wave_vectors, data

    def _get_interpolated_data(
        self, wave_vectors, data, crop_zero=False, kind="cubic"
    ):
        slices = self._get_slices(unfold=self.unfold, hse=self.hse)
        data_shape = data.shape
        if len(data_shape) == 1:
            data = [data[i] for i in slices]
        else:
            data = [data[:, i] for i in slices]

        wave_vectors = [wave_vectors[i] for i in slices]

        if len(data_shape) == 1:
            fs = [
                interp1d(i, j, kind=kind, axis=0)
                for (i, j) in zip(wave_vectors, data)
            ]
        else:
            fs = [
                interp1d(i, j, kind=kind, axis=1)
                for (i, j) in zip(wave_vectors, data)
            ]

        new_wave_vectors = [
            np.linspace(wv.min(), wv.max(), self.new_n) for wv in wave_vectors
        ]
        data = np.hstack([f(wv) for (f, wv) in zip(fs, new_wave_vectors)])
        wave_vectors = np.hstack(new_wave_vectors)

        if crop_zero:
            data[np.where(data < 0)] = 0

        return wave_vectors, data

    def _filter_bands(self, erange):
        eigenvalues = self.eigenvalues
        where = (eigenvalues >= np.min(erange) - 1) & (
            eigenvalues <= np.max(erange) + 1
        )
        is_true = np.sum(np.isin(where, True), axis=1)
        bands_in_plot = is_true > 0

        return bands_in_plot

    def _add_legend(self, ax, names, colors, fontsize=10, markersize=4):
        legend_lines = []
        legend_labels = []
        for name, color in zip(names, colors):
            legend_lines.append(
                plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    markersize=markersize,
                    linestyle="",
                    color=color,
                )
            )
            legend_labels.append(f"${name}$")

        leg = ax.get_legend()

        if leg is None:
            handles = legend_lines
            labels = legend_labels
        else:
            handles = [l._legmarker for l in leg.legendHandles]
            labels = [text._text for text in leg.texts]
            handles.extend(legend_lines)
            labels.extend(legend_labels)

        ax.legend(
            handles,
            labels,
            ncol=1,
            loc="upper left",
            fontsize=fontsize,
            bbox_to_anchor=(1, 1),
            borderaxespad=0,
            frameon=False,
            handletextpad=0.1,
        )

    def plot_plain(self,
        color="black",
        erange=[-6, 6],
        linewidth=1.25,
        scale_factor=20,
        linestyle="-",
        heatmap=False,
        bins=400,
        sigma=3,
        cmap="hot",
        vlinecolor="black",
        powernorm=False,
        gamma=0.5,
        projection=None,
        
        highlight_band=False,
        highlight_band_color="red",
        band_index=None,
        
        sp_color="red",
        sp_scale_factor=5,
        ax=None):
        # ------------ ax, fig objects -----------
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        # ------  ------
        bands_in_plot = self._filter_bands(erange=erange)
        slices = self._get_slices(unfold=self.unfold, hse=self.hse)
        wave_vector_segments = self._get_k_distance()

        # if self.soc_axis is not None and self.lsorbit:
        #     color = "black"
        #     linestyle = "-"

        if self.soc_axis is not None and self.lsorbit:
            if self.unfold:
                K_indices = np.array(self.K_indices[0], dtype=int)
                spin_projection_full_k = self.spin_projections[:, K_indices]
            else:
                spin_projection_full_k = self.spin_projections

        if self.custom_kpath is not None:
            kpath_inds = self.custom_kpath_inds
            kpath_flip = self.custom_kpath_flip
        else:
            kpath_inds = range(len(slices))
            kpath_flip = [False for _ in range(len(slices))]

        for i, f, wave_vectors in zip(
            kpath_inds, kpath_flip, wave_vector_segments
        ):
            if f:
                eigenvalues = np.flip(
                    self.eigenvalues[bands_in_plot, slices[i]], axis=1
                )
                if self.soc_axis is not None and self.lsorbit:
                    spin_projections = np.flip(
                        spin_projection_full_k[bands_in_plot, slices[i]],
                        axis=1,
                    )
            else:
                eigenvalues = self.eigenvalues[bands_in_plot, slices[i]]
                if self.soc_axis is not None and self.lsorbit:
                    spin_projections = spin_projection_full_k[
                        bands_in_plot, slices[i]
                    ]

            if highlight_band:
                if band_index is not None:
                    if type(band_index) == int:
                        highlight_eigenvalues = self.eigenvalues[
                            int(band_index), slices[i]
                        ]
                    else:
                        highlight_eigenvalues = self.eigenvalues[
                            band_index, slices[i]
                        ]

            wave_vectors_for_kpoints = wave_vectors

            if self.interpolate:
                (
                    wave_vectors,
                    eigenvalues,
                ) = self._get_interpolated_data_segment(
                    wave_vectors_for_kpoints,
                    eigenvalues,
                )
                if self.soc_axis is not None and self.lsorbit:
                    _, spin_projections = self._get_interpolated_data_segment(
                        wave_vectors_for_kpoints,
                        spin_projections,
                        crop_zero=True,
                        kind="linear",
                    )

                if highlight_band:
                    if band_index is not None:
                        (
                            _,
                            highlight_eigenvalues,
                        ) = self._get_interpolated_data_segment(
                            wave_vectors_for_kpoints,
                            highlight_eigenvalues,
                        )

            eigenvalues_ravel = np.ravel(
                np.c_[eigenvalues, np.empty(eigenvalues.shape[0]) * np.nan]
            )
            wave_vectors_tile = np.tile(
                np.append(wave_vectors, np.nan), eigenvalues.shape[0]
            )

            if self.soc_axis is not None and self.lsorbit:
                #  spin_cmap = self._alpha_cmap(color=spin_projection_color, repeats=1)
                spin_projections_ravel = np.ravel(
                    np.c_[
                        spin_projections,
                        np.empty(spin_projections.shape[0]) * np.nan,
                    ]
                )
                #  spin_colors = [spin_cmap(s) for s in spin_projections_ravel]

            if self.unfold:
                spectral_weights = self.spectral_weights[
                    bands_in_plot, slices[i]
                ]
                if f:
                    spectral_weights = np.flip(spectral_weights, axis=1)
                #  spectral_weights = spectral_weights / np.max(spectral_weights)

                if highlight_band:
                    if band_index is not None:
                        highlight_spectral_weights = self.spectral_weights[
                            int(band_index), slices[i]
                        ]

                if self.interpolate:
                    _, spectral_weights = self._get_interpolated_data_segment(
                        wave_vectors_for_kpoints,
                        spectral_weights,
                        crop_zero=True,
                        kind="linear",
                    )

                    if highlight_band:
                        if band_index is not None:
                            (
                                _,
                                highlight_spectral_weights,
                            ) = self._get_interpolated_data_segment(
                                wave_vectors_for_kpoints,
                                highlight_spectral_weights,
                                crop_zero=True,
                                kind="linear",
                            )

                spectral_weights_ravel = np.ravel(
                    np.c_[
                        spectral_weights,
                        np.empty(spectral_weights.shape[0]) * np.nan,
                    ]
                )

                if heatmap:
                    self._heatmap(
                        ax=ax,
                        wave_vectors=wave_vectors,
                        eigenvalues=eigenvalues,
                        weights=spectral_weights,
                        sigma=sigma,
                        cmap=cmap,
                        bins=bins,
                        projection=projection,
                        powernorm=powernorm,
                        gamma=gamma,
                    )
                else:
                    ax.scatter(
                        wave_vectors_tile,
                        eigenvalues_ravel,
                        c=color,
                        ec=[(1, 1, 1, 0)],
                        s=scale_factor * spectral_weights_ravel,
                        zorder=0,
                    )
                    if highlight_band:
                        if band_index is not None:
                            if type(band_index) == int:
                                ax.scatter(
                                    wave_vectors,
                                    highlight_eigenvalues,
                                    c=highlight_band_color,
                                    ec=[(1, 1, 1, 0)],
                                    s=scale_factor
                                    * highlight_spectral_weights,
                                    zorder=100,
                                )
                            else:
                                ax.scatter(
                                    np.tile(
                                        np.append(wave_vectors, np.nan),
                                        highlight_eigenvalues.shape[0],
                                    ),
                                    np.ravel(
                                        np.c_[
                                            highlight_eigenvalues,
                                            np.empty(
                                                highlight_eigenvalues.shape[0]
                                            )
                                            * np.nan,
                                        ]
                                    ),
                                    c=highlight_band_color,
                                    ec=[(1, 1, 1, 0)],
                                    s=scale_factor
                                    * np.ravel(highlight_spectral_weights),
                                    zorder=100,
                                )
                    if self.soc_axis is not None and self.lsorbit:
                        ax.scatter(
                            wave_vectors_tile,
                            eigenvalues_ravel,
                            s=spectral_weights_ravel
                            * sp_scale_factor
                            * spin_projections_ravel,
                            c=sp_color,
                            zorder=100,
                        )
            else:
                if heatmap:
                    self._heatmap(
                        ax=ax,
                        wave_vectors=wave_vectors,
                        eigenvalues=eigenvalues,
                        weights=np.ones(eigenvalues.shape),
                        sigma=sigma,
                        cmap=cmap,
                        bins=bins,
                        projection=projection,
                        powernorm=powernorm,
                        gamma=gamma,
                    )
                else:
                    ax.plot(
                        wave_vectors_tile,
                        eigenvalues_ravel,
                        color=color,
                        linewidth=linewidth,
                        linestyle=linestyle,
                        zorder=0,
                    )
                    if highlight_band:
                        if band_index is not None:
                            if type(band_index) == int:
                                ax.plot(
                                    wave_vectors,
                                    highlight_eigenvalues,
                                    color=highlight_band_color,
                                    linewidth=linewidth,
                                    linestyle=linestyle,
                                    zorder=100,
                                )
                            else:
                                ax.plot(
                                    np.tile(
                                        np.append(wave_vectors, np.nan),
                                        highlight_eigenvalues.shape[0],
                                    ),
                                    np.ravel(
                                        np.c_[
                                            highlight_eigenvalues,
                                            np.empty(
                                                highlight_eigenvalues.shape[0]
                                            )
                                            * np.nan,
                                        ]
                                    ),
                                    color=highlight_band_color,
                                    linewidth=linewidth,
                                    linestyle=linestyle,
                                    zorder=100,
                                )
                    if self.soc_axis is not None and self.lsorbit:
                        ax.scatter(
                            wave_vectors_tile,
                            eigenvalues_ravel,
                            s=sp_scale_factor * spin_projections_ravel,
                            c=sp_color,
                            zorder=100,
                        )

        if self.hse:
            self._get_kticks_hse(
                ax=ax,
                wave_vectors=np.concatenate(self._get_k_distance()),
                kpath=self.kpath,
                vlinecolor=vlinecolor,
            )
        elif self.unfold:
            self._get_kticks_unfold(
                ax=ax,
                wave_vectors=np.concatenate(self._get_k_distance()),
                vlinecolor=vlinecolor,
            )
        else:
            self._get_kticks(
                ax=ax,
                wave_vectors=np.concatenate(self._get_k_distance()),
                vlinecolor=vlinecolor,
            )

        ax.set_xlim(0, np.concatenate(self._get_k_distance()).max())

        return fig, ax    