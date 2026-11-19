import sys, os, io, re
sys.path.append(os.path.join(os.path.split(__file__)[0], '../cell2mol'))

import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import cell2mol
#from cell2mol.refcell import process_refcell
#from cell2mol.unitcell import process_unitcell
#from cell2mol.xyz_molecule import get_molecule
#from cell2mol.read_write import print_molecule

#from cell2mol.read_write import savemolecules
#from cell2mol.readwrite import readinfo, savemolecules

#from cell2mol.final_c2m_driver import handle_cif_file
#from cell2mol.cif2info import cif_2_info
#from cell2mol.c2m_module import save_cell, cell2mol

from pathlib import Path

from cell2mol.utils import config

from cell2mol.write_results import (
    write_unique_species,
    write_plausible_charges,
    get_reference_error_message,
    get_reference_error_messages,
    get_unitcell_error_message,
    get_reference_warning_messages,
)

def get_refcell_interpretation(refcell):
    main = []
    main.append('Total charge comparison result: {}'.format(getattr(refcell, 'total_charge_comparison')))
    main.append('Metal oxidation state comparison result: {}'.format(getattr(refcell, 'metal_os_comparison')))
    reported_os = getattr(refcell, "reported_metal_os", None)
    matched_os = getattr(refcell, "reported_metal_os_matched", None)
    os_confidence = getattr(refcell, "metal_os_match_confidence", None)
    confidence_str = f"{os_confidence:.2f}" if os_confidence is not None else "N/A"
    main.append("  - Reported:   {}".format(reported_os))
    main.append("  - Matched:    {}".format(matched_os))
    main.append("  - Confidence: {}".format(confidence_str))
    main.append(" ")

    cell_mol_info = get_cell_molecules_info(refcell)
    main = main + cell_mol_info
    main.append(" ")

    unique_species = get_unique_species(refcell)
    main = main + unique_species
    main.append(" ")

    #fix to check what to print here
    #plausible_charges = get_plausible_charges(refcell)
    #main = main + plausible_charges
    #main.append(" ")

    return(main)

def get_plausible_charges(object):
    plausible_charges = []
    if object.species_list is not None:
        plausible_charges.append("Possible charges of species in {object.subtype}:")
    for specie in object.species_list:
        info = f"\tunique_index={specie.unique_index} {specie.formula} ({specie.subtype})"

        if specie.subtype == "metal":
            info += f" coord_sphere_formula={getattr(specie, 'coord_sphere_formula', 'N/A')}"

        if specie.subtype == "metal":
            plausible = getattr(specie, "plausible_os", None)
        else:
            plausible = getattr(specie, "plausible_charge_states", None)

        if plausible is not None:
            if specie.subtype == "metal":
                info += f" plausible_os={plausible}"
            else:
                info += f"\n\tplausible_charge_states={plausible}"

        else:
            info += f" {_no_charge_state_reason(specie)}"

        plausible_charges.append(info)

    else:
        plausible_charges.append("No species list found in the cell object.")


    inconsistent = getattr(object, "inconsistent_plausible_charges", None)
    if inconsistent:
        plausible_charges.append("WARNING: plausible charges disagree between copies of the same specie:")

        for unique_index, variants in sorted(inconsistent.items()):
            formula = next(
                    (
                        spec.formula
                        for spec in (object.species_list or [])
                        if getattr(spec, "unique_index", None) == unique_index
                    ),
                    "?",
                )
            plausible_charges.append(f"\tunique_index={unique_index} {formula} enumerated to {variants}")





    return(plausible_charges)


def get_unique_species(object):
    uniqSpe = []
    unique_species = getattr(object, "unique_species", None)
    if not unique_species:
        uniqSpe.append(f"No unique species found in the {object.subtype} object.")
        return(uniqSpe)

    uniqSpe.append(f"Unique Species in {object.subtype}:")

    for specie in unique_species:
        parts = [
                f"unique_index={specie.unique_index}",
                f"{specie.formula}",
                f"({specie.subtype})",
                ]

        if specie.subtype == "metal":
            parts.append(f"coord_sphere_formula={specie.coord_sphere_formula}")
            if getattr(specie, "charge", None) is not None:
                parts.append(f"charge={specie.charge}")

        elif specie.subtype == "ligand":
            parts.append(f"denticity={specie.denticity}")
            if specie.is_haptic:
                parts.append(f"haptic_type={specie.haptic_type}")
            if getattr(specie, "smiles", None) is not None:
                parts.append(f"smiles={specie.smiles}")
            if getattr(specie, "totcharge", None) is not None:
                parts.append(f"totcharge={specie.totcharge}")
            if getattr(specie, "groups", None) is not None:
                parts.append(f"groups={[group.formula for group in specie.groups]}")

        else:
            if getattr(specie, "smiles", None) is not None:
                parts.append(f"smiles={specie.smiles}")
            if getattr(specie, "totcharge", None) is not None:
                parts.append(f"totcharge={specie.totcharge}")
        uniqSpe.append("\t" + " ".join(parts))

    return(uniqSpe)


def get_cell_molecules_info(cell):
    cell_mol_info = []
    if cell.subtype == "reference":
        cell_type = "Reference"
        molecule_list = getattr(cell, "refmoleclist", None)
    elif cell.subtype == "unitcell":
        cell_type = "Unitcell"
        molecule_list = getattr(cell, "moleclist", None)
    else:
        return

    if molecule_list is None:
        cell_mol_info.append('No molecules found in the {} cell.'.format(cell_type))
        return
    
    cell_mol_info.append("Molecules in {}:".format(cell.subtype))

    for i, mol in enumerate(molecule_list):
        get_molecule_info(mol, cell_mol_info, index=i)

    return(cell_mol_info)


def get_molecule_info(mol, cell_mol_info, index=None):
    if mol is None:
        cell_mol_info.append("No molecule object.")
        return

    prefix = "Molecule {}:".format(index)
    cell_mol_info.append("{} {}".format(prefix, mol.formula))

    if mol.contains_metal:
        cell_mol_info.append("(TM Complex)")
    if mol.has_s_block_metal:
        cell_mol_info.append("(Complex with Alkali or Alkaline metals)")
    if mol.is_post_transition_metal_complex:
        cell_mol_info.append("(Complex with Post-Transition metals)")
    if mol.is_metal_free_molecule:
        cell_mol_info.append("(Non-complex)")

    if mol.totcharge is not None:
        cell_mol_info.append(f"totcharge={mol.totcharge}")
    if getattr(mol, "totcharge_cif", None) is not None:
        cell_mol_info.append(f"totcharge_cif={mol.totcharge_cif}")
    if getattr(mol, "totcharge_agree", None) is not None:
        cell_mol_info.append(f"totcharge_agree={mol.totcharge_agree}")
    if getattr(mol, "spin", None) is not None:
        cell_mol_info.append(f"spin_multiplicity={mol.spin}")
    if mol.smiles is not None:
        cell_mol_info.append(f"smiles={mol.smiles}")


    if not mol.is_metal_free_molecule and hasattr(mol, "metals"):
        for met in mol.metals:
            met_info = f"\t{met.formula} ({met.subtype})"

            labels = getattr(met, "atom_site_label", None)
            if labels:
                met_info += f" atom_site_label={labels}"

            if met.charge is not None:
                met_info += f" metal_OS={met.charge}"
            elif getattr(met, "plausible_os", None) is not None:
                met_info += f" metal_plausible_OS={met.plausible_os}"

            if getattr(met, "spin", None) is not None:
                met_info += f" metal_spin={met.spin}"

            cell_mol_info.append(met_info)


            if getattr(met, "coord_sphere_formula", None):
                cell_mol_info.append(f"\t|--coord_sphere_formula={met.coord_sphere_formula}")

            if all(
                    hasattr(met, attr)
                    for attr in ("coord_nr", "coord_geometry", "geom_deviation")
                ):
                cell_mol_info.append(f"\t|--coord_nr={met.coord_nr} coord_geometry={met.coord_geometry} geom_deviation={met.geom_deviation}")

            if getattr(met, "groups", None):
                for group in met.groups:
                    group_info = f"\t|--(group) {group.labels}"
                    group_atom_site_labels = [
                            a.atom_site_label
                            for a in group.atoms
                            if a.atom_site_label is not None
                        ]
                    if group_atom_site_labels:
                        group_info += f" atom_site_labels={group_atom_site_labels}"
                    if hasattr(group, "denticity"):
                        group_info += f" denticity={group.denticity}"
                    if getattr(group, "is_haptic", False):
                        group_info += f" haptic_type={group.haptic_type}"

                    cell_mol_info.append(group_info)


                # Metal-Metal Bonds
                if hasattr(met, "metals") and hasattr(met, "coord_nr_with_metal_bonds"):
                    bonded_metals = [m.label for m in met.metals]
                    if bonded_metals:
                        cell_mol_info.append(
                                f"\t|--bonded metals={bonded_metals} coord_nr_with_metal_bonds={met.coord_nr_with_metal_bonds} "
                                f"coord_geometry_with_metal_bonds={getattr(met, 'coord_geometry_with_metal_bonds', 'N/A')} "
                                f"geom_deviation_with_metal_bonds={getattr(met, 'geom_deviation_with_metal_bonds', 'N/A')}"
                                )

    # --------------------
    # Ligands
    # --------------------

    if not mol.is_metal_free_molecule and getattr(mol, "ligands", None):
        for lig in mol.ligands:
            lig_info = f"\t{lig.formula} ({lig.subtype})"
            for attr in ("smiles", "denticity", "totcharge"):
                val = getattr(lig, attr, None)
                if val is not None:
                    lig_info += f" {attr}={val}"
            if (getattr(lig, "plausible_charge_states", None) is not None
                and getattr(lig, "totcharge", None) is None
                ):
                status = "Exists" if lig.plausible_charge_states else "Does not exist"
                lig_info += f" lig.plausible_charge_states {status}"
            cell_mol_info.append(lig_info)


            if getattr(lig, "groups", None):
                for group in lig.groups:
                    group_info = f"\t|--(group) {group.labels}"
                    group_atom_site_labels = [
                            a.atom_site_label
                            for a in group.atoms
                            if a.atom_site_label is not None
                        ]
                    if group_atom_site_labels:
                        group_info += f" atom_site_labels={group_atom_site_labels}"
                    if hasattr(group, "denticity"):
                        group_info += f" denticity={group.denticity}"
                    if getattr(group, "is_haptic", False):
                        group_info += f" haptic_type={group.haptic_type}"

                    if getattr(group, "metals", None):
                        labels = [
                                getattr(m, "atom_site_label", None) or m.label
                                for m in group.metals
                            ]
                        group_info += f" connected_metals={labels}"
                    cell_mol_info.append(group_info)

    return(cell_mol_info)



def get_refcell_error(refcell):
    refcell_error = []
    if refcell.error_cases:
        all_codes = refcell.error_cases_all or {}
        for err_mode, code in refcell.error_cases.items():
            msg = get_reference_error_messages(all_codes.get(err_mode), code)
            refcell_error.append(f"Reference Error (mode={err_mode}): {msg}")

        return(refcell_error)

def get_unitcell_error(unitcell):
    unitcell_error = []
    for err_mode in unitcell.error_cases.keys():
        code = unitcell.error_cases.get(err_mode, 0)
        msg = get_unitcell_error_message(code)
        unitcell_error.append(f"UnitCell Error (mode={err_mode}): {msg}")

    return(unitcell_error)

def print_molecule(mol, name, ext, folder):
    filename = str(folder) + "/" + str(name) + "." + str(ext)

    if ext == "xyz" or ext == "txt":
        with open(filename, "w") as fil:

            # XYZ
            if ext == "xyz":
                print(mol.natoms, file=fil)
                if hasattr(mol, "totcharge") and hasattr(mol, "spin"):
                    print(mol.totcharge, mol.spin, file=fil)
                if hasattr(mol, "totcharge") and not hasattr(mol, "spin"):
                    print(mol.totcharge, "SPIN", file=fil)
                # if not hasattr(mol, 'totcharge') and not hasattr(mol, 'spin'):
                else:
                    print("", file=fil)

                for a in mol.atoms:
                    print("%s   %.6f   %.6f   %.6f" % (a.label, a.coord[0], a.coord[1], a.coord[2]),file=fil)

            # TXT
            elif ext == "txt":
                print(vars(mol), file=fil)

    elif ext == "gmol" or ext == "mol" or ext == "npy" or ext == "dict":
        with open(filename, "wb") as fil:

            # GMOL
            if ext == "gmol":
                pickle.dump(mol, fil)

            # MOL
            elif ext == "mol":
                pickle.dump(mol.object, fil)

            # NPY
            elif ext == "npy":
                np.save(filename, mol)

            # DICT
            elif ext == "dict":
                mydict = vars(mol)
                pickle.dump(mydict, fil)

    else:
        print(ext, "not found as a valid print extension in print_molecule")


def get_error_code(output_dir, prefix):
    """
    Finds files like:
      unitcell_error_0.out
      reference_error_2.out

    Returns the error number as int, or None if no file exists.
    """
    output_dir = Path(output_dir)

    pattern = re.compile(rf"{prefix}_error_(\d+)\.out$")

    for path in output_dir.glob(f"{prefix}_error_*.out"):
        match = pattern.match(path.name)
        if match:
            return int(match.group(1))

    return None

def save_cell(cell: object, ext: str, output_dir: str, refcode: str):
    #taken from old cell2mol version. 
    cellpath = os.path.join(output_dir, "Cell_{}.gmol".format(refcode))
    with open(cellpath, "wb") as fil:
        if ext == "gmol":
            pickle.dump(cell, fil)
        else:
            print(ext, "not found as a valid print extension in print_molecule")



def savemolecules_tools(moleclist, output_dir, print_types, option_print_repeated=True):
    # DEFAULTS
    print_xyz = True
    print_gmol = True
    print_npy = False
    print_mol = False
    print_txt = False
    print_dict = False

    if "xyz" not in print_types:
        print_xyz = False
    if "gmol" not in print_types:
        print_gmol = False
    if "npy" in print_types:
        print_npy = True
    if "mol" in print_types:
        print_mol = True
    if "txt" in print_types:
        print_txt = True
    if "dict" in print_types:
        print_dict = True

    printedmolecs = []

    for idx, mol in enumerate(moleclist):

        if mol.contains_metal:
            molName = mol.get_parent("reference").name + "_Complex_"+str(idx) 
        else:
            molName = mol.get_parent("reference").name + "_Other_"+ str(idx)
        #return molName

        #shalliprint = False

        #if any((mol.elemcountvec == pmol.elemcountvec).all() for pmol in printedmolecs):
        #    shalliprint = False
        #else:
        #    shalliprint = True
        #if option_print_repeated:  # Overwrites decision if the user decides so
        #    shalliprint = True

        #if shalliprint:
        #    printedmolecs.append(mol)

        if print_xyz:
            print_molecule(mol, molName, "xyz", output_dir)
        if print_gmol:
            print_molecule(mol, molName, "gmol", output_dir)
        if print_npy:
            print_molecule(mol, molName, "npy", output_dir)
        if print_txt:
            print_molecule(mol, molName, "txt", output_dir)
        if print_dict:
            print_molecule(mol, molName, "dict", output_dir)
        if hasattr(mol, "object"):
            if print_mol:
                print_molecule(mol, molName, "mol", output_dir)

    


ELEMENTS = [  # thanks pyscf
    'X',  # Ghost
    'H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' , 'Ca',
    'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I' , 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
]

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = io.StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout


def cell_cmp_lut(cell):
    names = {}
    #for i_mol,mol in enumerate(cell.refmoleclist):
    for i_mol,mol in enumerate(cell.unitcell.moleclist):
        #if mol.type == 'Complex':
        if mol.contains_metal :
            #for i_mtl, mtl in enumerate(mol.metalist):
            for i_mtl, mtl in enumerate(mol.metals):
                mlabel = mtl.label
                #if mtl.totcharge > 0:
                if mtl.charge > 0:
                    mlabel = f'[{mlabel:s}+{mtl.charge:d}]'
                #elif mtl.totcharge < 0:
                elif mtl.charge < 0:
                    mlabel = f'[{mlabel:s}-{-mtl.charge:d}]'
                names.setdefault(mlabel, [])
                names[mlabel].append((i_mol,'m',i_mtl))
            #for i_lig, lig in enumerate(mol.ligandlist):
            for i_lig, lig in enumerate(mol.ligands):
                names.setdefault(lig.smiles, [])
                names[lig.smiles].append((i_mol,'l',i_lig))
        else:
            names.setdefault(mol.smiles, [])
            names[mol.smiles].append((i_mol,'t'))
    return names


def cell_cmp_lut_base(refCell):
    """
    Build compound lookup table for REFERENCE mode.

    In reference mode, ligand/molecule .smiles may still be None because
    no final charge state has been selected. Therefore we use unique_index
    to identify/group species instead of SMILES.

    Returned tuples keep the same structure expected by cell_to_svgs_base:
        (i_mol, 'm', i_mtl)
        (i_mol, 'l', i_lig)
        (i_mol, 't')
    """
    names = {}

    for i_mol, mol in enumerate(refCell.refmoleclist):

        if mol.contains_metal:

            # Metals
            for i_mtl, mtl in enumerate(mol.metals):
                # Include type in the key so a metal label cannot collide
                # with a ligand/molecule identifier.
                key = ('m', mtl.label)

                names.setdefault(key, [])
                names[key].append((i_mol, 'm', i_mtl, mtl.formula))

            # Ligands
            for i_lig, lig in enumerate(mol.ligands):

                # Reference mode should have unique_index assigned by
                # get_unique_species().
                unique_index = getattr(lig, 'unique_index', None)

                if unique_index is None:
                    raise ValueError(
                        "Reference ligand has no unique_index: "
                        f"i_mol={i_mol}, "
                        f"i_lig={i_lig}, "
                        f"formula={getattr(lig, 'formula', None)!r}, "
                        f"origin={getattr(lig, 'origin', None)!r}, "
                        f"ligand={lig!r}"
                    )

                key = ('l', unique_index)

                names.setdefault(key, [])
                names[key].append((i_mol, 'l', i_lig, lig.formula))

        else:
            # Non-complex molecule
            unique_index = getattr(mol, 'unique_index', None)

            if unique_index is None:
                raise ValueError(
                    "Reference molecule has no unique_index: "
                    f"i_mol={i_mol}, "
                    f"formula={getattr(mol, 'formula', None)!r}, "
                    f"origin={getattr(mol, 'origin', None)!r}, "
                    f"molecule={mol!r}"
                )

            key = ('t', unique_index)

            names.setdefault(key, [])
            names[key].append((i_mol, 't', 0, mol.formula))

    return names


def _reference_rdkit_obj(specie):
    """
    Return an RDKit Mol suitable for depiction of a reference-mode specie.

    Reference species do not necessarily have specie.rdkit_obj populated,
    because reference mode enumerates plausible charge states without
    selecting a final one.

    We therefore obtain the RDKit molecule from a valid plausible
    charge state.

    A COPY of the RDKit Mol is returned because Compute2DCoords modifies
    the molecule by adding coordinates.
    """

    states = getattr(specie, 'plausible_charge_states', None)

    if not states:
        raise ValueError(
            "Reference specie has no plausible charge states: "
            f"formula={getattr(specie, 'formula', None)!r}, "
            f"origin={getattr(specie, 'origin', None)!r}, "
            f"specie={specie!r}"
        )

    valid_states = []

    for cs in states:

        # Some Cell2Mol objects expose status; if it does not exist,
        # treat the state as usable unless its RDKit object is missing.
        status = getattr(cs, 'is_valid', True)
        rd = getattr(cs, 'species_rdkit_obj', None)

        if status and rd is not None:
            valid_states.append(cs)

    if not valid_states:
        raise ValueError(
            "Reference specie has plausible charge states, but none "
            "contains a usable RDKit object: "
            f"formula={getattr(specie, 'formula', None)!r}, "
            f"origin={getattr(specie, 'origin', None)!r}, "
            f"n_states={len(states)}, "
            f"specie={specie!r}"
        )

    # IMPORTANT:
    # Reference mode has not selected a definitive charge state.
    # For visualization only, use the first valid plausible state.
    cs = valid_states[0]

    return Chem.Mol(cs.species_rdkit_obj)

#def cell_cmp_lut_base(refCell):
#    names = {}
#
#    for i_mol, mol in enumerate(refCell.refmoleclist):
#
#        if mol.contains_metal:
#
#            for i_mtl, mtl in enumerate(mol.metals):
#                key = ("metal", mtl.label)
#                names.setdefault(key, [])
#                names[key].append((i_mol, 'm', i_mtl))
#
#            for i_lig, lig in enumerate(mol.ligands):
#                key = ("ligand", lig.unique_index)
#                names.setdefault(key, [])
#                names[key].append((i_mol, 'l', i_lig))
#
#        else:
#            key = ("molecule", mol.unique_index)
#            names.setdefault(key, [])
#            names[key].append((i_mol, 't'))
#
#    return names

#def cell_cmp_lut_base(refCell):
#    names = {}
#
#    for i_mol, mol in enumerate(refCell.refmoleclist):
#
#        if mol.contains_metal:
#
#            for i_mtl, mtl in enumerate(mol.metals):
#                mlabel = mtl.label
#
#                names.setdefault(mlabel, [])
#                names[mlabel].append((i_mol, 'm', i_mtl))
#
#            for i_lig, lig in enumerate(mol.ligands):
#
#                if lig.smiles is None:
#                    raise ValueError(
#                        f"Ligand has no SMILES: "
#                        f"i_mol={i_mol}, "
#                        f"i_lig={i_lig}, "
#                        f"rdkit_obj={lig.rdkit_obj!r}, "
#                        f"ligand={lig!r}"
#                    )
#
#                if lig.rdkit_obj is None:
#                    raise ValueError(
#                        f"Ligand has no RDKit object: "
#                        f"i_mol={i_mol}, "
#                        f"i_lig={i_lig}, "
#                        f"smiles={lig.smiles!r}, "
#                        f"ligand={lig!r}"
#                    )
#
#                names.setdefault(lig.smiles, [])
#                names[lig.smiles].append((i_mol, 'l', i_lig))
#
#        else:
#
#            if mol.smiles is None:
#                raise ValueError(
#                    f"Molecule has no SMILES: "
#                    f"i_mol={i_mol}, mol={mol!r}"
#                )
#
#            names.setdefault(mol.smiles, [])
#            names[mol.smiles].append((i_mol, 't'))
#
#    return names

#def cell_cmp_lut_base(refCell):
#    names = {}
#    #for i_mol,mol in enumerate(cell.refmoleclist):
#    for i_mol,mol in enumerate(refCell.refmoleclist):
#        #if mol.type == 'Complex':
#        if mol.contains_metal :
#            #for i_mtl, mtl in enumerate(mol.metalist):
#            for i_mtl, mtl in enumerate(mol.metals):
#                mlabel = mtl.label
#                #charge not predicted in ref mode
#                #if mtl.totcharge > 0:
#                #if mtl.charge > 0:
#                #    mlabel = f'[{mlabel:s}+{mtl.charge:d}]'
#                ##elif mtl.totcharge < 0:
#                #elif mtl.charge < 0:
#                #    mlabel = f'[{mlabel:s}-{-mtl.charge:d}]'
#                names.setdefault(mlabel, [])
#                names[mlabel].append((i_mol,'m',i_mtl))
#            #for i_lig, lig in enumerate(mol.ligandlist):
#            for i_lig, lig in enumerate(mol.ligands):
#                names.setdefault(lig.smiles, [])
#                names[lig.smiles].append((i_mol,'l',i_lig))
#        else:
#            names.setdefault(mol.smiles, [])
#            names[mol.smiles].append((i_mol,'t'))
#    return names


def cell_get_metal_desc(cell, cmplut):
    res = []
    for name, lst in cmplut.items():
        tpl = lst[0]
    #for mol in cell.moleclist:
        mol = cell.unitcell.moleclist[tpl[0]]
        if tpl[1]=='m':
            mtl = mol.metals[tpl[2]]
            res.append('<p>Metal center: {0:s}<br/>predicted charge: {1:+d}</p>'.format(
                mtl.label,
                mtl.charge,
            ))
        else:
            res.append("")
    return res
    

def savemolecules_reference_xyz(refCell, output_dir):
    for idx, mol in enumerate(refCell.refmoleclist):
        if mol.contains_metal:
            molName = mol.get_parent("reference").name + "_Complex_" + str(idx)
        else:
            molName = mol.get_parent("reference").name + "_Other_" + str(idx)

        filename = os.path.join(output_dir, molName + ".xyz")
        with open(filename, "w") as fil:
            print(mol.natoms, file=fil)
            print(molName, file=fil)
            for atom in mol.atoms:
                x, y, z = atom.coord
                print("%s   %.6f   %.6f   %.6f" % (atom.label, x, y, z), file=fil)


def reference_plausible_charges(refCell, cmplut):
    res = []
    for name, lst in cmplut.items():
        tpl = lst[0]
        mol = refCell.refmoleclist[tpl[0]]
        charges = []

        if tpl[1] == 'm':
            specie = mol.metals[tpl[2]]
            states = getattr(specie, 'plausible_os', None)
            if states is not None:
                for charge in states:
                    if isinstance(charge, (int, np.integer)):
                        charges.append(int(charge))
        else:
            if tpl[1] == 'l':
                specie = mol.ligands[tpl[2]]
            else:
                specie = mol

            states = getattr(specie, 'plausible_charge_states', None)
            if states is not None:
                for state in states:
                    if not getattr(state, 'is_valid', True):
                        continue
                    charge = None
                    for attr in ('species_total_charge', 'specie_totcharge', 'specie_charge', 'totcharge', 'charge'):
                        value = getattr(state, attr, None)
                        if isinstance(value, (int, np.integer)):
                            charge = int(value)
                            break
                    if charge is None:
                        rd = getattr(state, 'species_rdkit_obj', None)
                        if rd is not None:
                            charge = int(Chem.GetFormalCharge(rd))
                    if charge is not None:
                        charges.append(charge)

        charges = sorted(set(charges))
        if charges:
            res.append(', '.join(f'{charge:+d}' if charge > 0 else str(charge) for charge in charges))
        else:
            res.append('Not available')

    return res


def cell_get_metal_desc_base(refCell, cmplut):
    res = []
    for name, lst in cmplut.items():
        tpl = lst[0]
    #for mol in cell.moleclist:
        mol = refCell.refmoleclist[tpl[0]]
        if tpl[1]=='m':
            mtl = mol.metals[tpl[2]]
            res.append('<p>Metal center: {0:s}</p>'.format(
                mtl.label,
            ))
        else:
            res.append("")
    return res
    
        
# def cell_to_string_xsf(cell, cmplut=None):
#     xsf = "CRYSTAL\nPRIMVEC\n"
#     for vec in cell.cellvec:
#         xsf += "{:f} {:f} {:f}\n".format(vec[0],vec[1],vec[2])
#     nat = sum(mol.natoms for mol in cell.moleclist)
#     xsf += "PRIMCOORD\n{:d} 1\n".format(nat)
#     for mol in cell.moleclist:
#         for Z,(x,y,z) in zip(mol.atnums, mol.coord):
#             xsf += "{:d} {:f} {:f} {:f}\n".format(Z,x,y,z)
#     return xsf

#def cell_to_string_xyz(cell, cmplut=None):
#    celldesc = "a={:f},b={:f},c={:f},alpha={:f},beta={:f},gamma={:f}".format(*tuple(cell.cellparam))
#    mols = []
#    # for mol in cell.moleclist:
#    #     atms = []
#    #     for Z,(x,y,z) in zip(mol.atnums, mol.coord):
#    #         atms.append(" {:s}    {:8f} {:8f} {:8f}\n".format(ELEMENTS[Z],x,y,z))
#    #     mols.append(
#    #         "{:d}\n{:s}\n".format(len(atms), mol.name) +
#    #         "".join(atms)
#    #     )
#    for name, elems in cmplut.items():
#        atms = []
#        for elem in elems:
#            mol = cell.moleclist[elem[0]]
#            if elem[1] == 't':
#                pass
#            elif elem[1] == 'l':
#                mol = mol.ligandlist[elem[2]]
#            elif elem[1] == 'm':
#                mol = mol.metalist[elem[2]]
#                x,y,z = mol.atom.coord
#                Z = mol.atom.atnum
#                atms.append(" {:s}    {:8f} {:8f} {:8f}\n".format(ELEMENTS[Z],x,y,z))
#                continue
#
#            for Z,(x,y,z) in zip(mol.atnums, mol.coord):
#                atms.append(" {:s}    {:8f} {:8f} {:8f}\n".format(ELEMENTS[Z],x,y,z))
#        mols.append(
#            "{:d}\n{:s}\n".format(len(atms), mol.name) +
#            "".join(atms)
#        )
#    return celldesc, "".join(mols)

def cell_to_string_xyz(cell, cmplut=None):
    celldesc = "{:f},{:f},{:f},{:f},{:f},{:f}".format(*tuple(cell.cell_param))

    allmols = []
    for idx, mols in enumerate(cell.unitcell.moleclist):                                                                                                   

        if mols.contains_metal:
            molName = mols.get_parent("reference").name + "_Complex_"+str(idx) 
        else:
            molName = mols.get_parent("reference").name + "_Other_"+ str(idx)

        atms = []
        for atoms in mols.atoms:                                                                                                  
            x,y,z = atoms.coord
            Z = atoms.label
            atms.append(" {:s}    {:8f} {:8f} {:8f}\\n".format(Z,x,y,z))
        allmols.append( "{:d}\\n {:s}\\n ".format(mols.natoms, molName) + "".join(atms) )

    return celldesc, "".join(allmols)

def refcell_to_string_xyz(cell, cmplut=None):
    celldesc = "{:f},{:f},{:f},{:f},{:f},{:f}".format(*tuple(cell.cell_param))

    allmols = []
    for idx, mols in enumerate(cell.refmoleclist):                                                                                                   

        if mols.contains_metal:
            molName = mols.get_parent("reference").name + "_Complex_"+str(idx) 
        else:
            molName = mols.get_parent("reference").name + "_Other_"+ str(idx)

        atms = []
        for atoms in mols.atoms:                                                                                                  
            x,y,z = atoms.coord
            Z = atoms.label
            atms.append(" {:s}    {:8f} {:8f} {:8f}\\n".format(Z,x,y,z))
        allmols.append( "{:d}\\n {:s}\\n ".format(mols.natoms, molName) + "".join(atms) )

    return celldesc, "".join(allmols)

def interpret_molecule_xyz(input_path, refcode, output_dir, total_charge):
    from cell2mol.process_xyz import interpret_molecule
    previous_bond_info = config.USE_BOND_INFO
    config.USE_BOND_INFO = False
    try:
        return interpret_molecule(input_path, refcode, total_charge, output_dir)
    finally:
        config.USE_BOND_INFO = previous_bond_info


def get_moleculeset_molecules(molset):
    if molset is None:
        raise ValueError('cell2mol returned no molecule. Check the XYZ input and cell2mol logs.')
    if isinstance(molset, (list, tuple)):
        return list(molset)
    for attr in ('moleclist', 'molecules', 'mol_list'):
        molecules = getattr(molset, attr, None)
        if molecules is not None:
            return list(molecules)
    if hasattr(molset, 'atoms'):
        return [molset]
    raise ValueError('The molecule-mode result does not contain a molecule list')


def make_moleculeset_cell(moleclist):
    from types import SimpleNamespace
    return SimpleNamespace(unitcell=SimpleNamespace(moleclist=moleclist), cell_param=(1., 1., 1., 90., 90., 90.))


def moleculeset_molecule_name(mol, idx, refcode):
    return refcode + ('_Complex_' if mol.contains_metal else '_Other_') + str(idx)


def moleculeset_to_string_xyz(moleclist, refcode):
    allmols = []
    for idx, mol in enumerate(moleclist):
        molName = moleculeset_molecule_name(mol, idx, refcode)
        atms = []
        for atom in mol.atoms:
            x, y, z = atom.coord
            atms.append(" {:s}    {:8f} {:8f} {:8f}\n".format(atom.label, x, y, z))
        allmols.append("{:d}\n {:s}\n".format(len(mol.atoms), molName) + "".join(atms))
    return "".join(allmols)


def savemolecules_moleculeset_xyz(moleclist, output_dir, refcode):
    for idx, mol in enumerate(moleclist):
        molName = moleculeset_molecule_name(mol, idx, refcode)
        print_molecule(mol, molName, 'xyz', output_dir)


def molecules_list_moleculeset(moleclist, refcode):
    jmol_list_pos = {}
    for idx, mol in enumerate(moleclist):
        molName = moleculeset_molecule_name(mol, idx, refcode)
        jmol_list_pos[molName] = ' select ' + ' or '.join(
            'within (0.1, {{{} {} {}}})'.format(*atom.coord) for atom in mol.atoms
        )
    return jmol_list_pos


def moleculeset_is_resolved(moleclist):
    for mol in moleclist:
        if getattr(mol, 'totcharge', None) is None:
            return False
        if getattr(mol, 'contains_metal', False):
            for metal in getattr(mol, 'metals', None) or []:
                if getattr(metal, 'charge', None) is None:
                    return False
            for ligand in getattr(mol, 'ligands', None) or []:
                if getattr(ligand, 'totcharge', None) is None or getattr(ligand, 'rdkit_obj', None) is None:
                    return False
        elif getattr(mol, 'rdkit_obj', None) is None:
            return False
    return True


def printing_text_moleculeset(moleclist, output):
    if moleculeset_is_resolved(moleclist):
        cell = make_moleculeset_cell(moleclist)
        return printing_text(cell, output)

    for idx, mol in enumerate(moleclist):
        prefix = '[Complex]' if getattr(mol, 'contains_metal', False) else '[Other]'
        output.append(f'{prefix} Formula : {mol.formula} (molecule : {idx})')
        charge = getattr(mol, 'totcharge', None)
        output.append('   Total charge: ' + (str(charge) if charge is not None else 'Not resolved'))
        spin = getattr(mol, 'spin', None)
        output.append('   Spin: ' + (str(spin) if spin is not None else 'Not available'))
        if getattr(mol, 'smiles', None):
            output.append('   Smiles: ' + mol.smiles)
        for metal in getattr(mol, 'metals', None) or []:
            output.append(f'   >> Metal : {metal.label}')
            oxidation = getattr(metal, 'charge', None)
            output.append('   Oxidation state: ' + (str(oxidation) if oxidation is not None else 'Not resolved'))
            for attr, label in (('coord_nr', 'Coordination number'),
                                ('coord_geometry', 'Coordination geometry'),
                                ('coord_sphere_formula', 'Coordination sphere formula'),
                                ('rel_metal_radius', 'Relative metal radius')):
                value = getattr(metal, attr, None)
                if value is not None:
                    output.append('   ' + label + ': ' + str(value))
        for ligand in getattr(mol, 'ligands', None) or []:
            output.append(f'   >> Ligand Formula : {ligand.formula}')
            charge = getattr(ligand, 'totcharge', None)
            output.append('   Charge: ' + (str(charge) if charge is not None else 'Not resolved'))
            for attr, label in (('denticity', 'Denticity'), ('haptic_type', 'Hapticity'), ('smiles', 'Smiles')):
                value = getattr(ligand, attr, None)
                if value is not None:
                    output.append('   ' + label + ': ' + str(value))
        output.append(' ')
    return output


def get_moleculeset_errors(molset):
    from cell2mol.write_results import get_molecule_error_message

    stages = [
        ('hydrogens', [('has_isolated_H', 1), ('missing_H_in_Water', 2),
                       ('missing_H_on_CoordDonor', 3), ('missing_H_in_Carbon', 4)]),
        ('plausible_charges', [('error_plausible_charges', 5)]),
        ('balance_charges', [('error_multiple_distrib', 6), ('error_empty_distrib', 7)]),
        ('charge_assignment', [('error_assign_charge', 8)]),
        ('bond_creation', [('error_create_bonds', 9)]),
        ('spin_assignment', [('error_get_spin', 10)]),
    ]
    errors = []
    for mode, attributes in stages:
        active = [get_molecule_error_message(code)
                  for attr, code in attributes if getattr(molset, attr, None) is True]
        flags = [getattr(molset, attr, None) for attr, code in attributes]
        if active:
            message = '; '.join(active)
        elif flags and all(flag is False for flag in flags):
            message = 'No errors found'
        elif flags and any(flag is False for flag in flags):
            message = 'No errors detected in evaluated checks; other checks not reported'
        else:
            message = 'Not evaluated or not reported by this cell2mol version'
        errors.append(f'Molecule Error (mode={mode}): {message}')

    overall = getattr(molset, 'error_case', None)
    if overall is None:
        overall_msg = 'Not reported'
    else:
        overall_msg = get_molecule_error_message(overall)
    errors.append(f'Molecule Error (mode=overall): {overall_msg}')
    return errors


def moleculeset_partial_species(moleclist):
    from html import escape
    compound_data = []
    jmol_list_species = {}
    entries = {}

    for idx, mol in enumerate(moleclist):
        items = []
        if getattr(mol, 'contains_metal', False):
            items += [('m', metal) for metal in getattr(mol, 'metals', None) or []]
            items += [('l', lig) for lig in getattr(mol, 'ligands', None) or []]
        else:
            items.append(('t', mol))

        for kind, specie in items:
            formula = getattr(specie, 'formula', None) or getattr(specie, 'label', 'Unknown')
            charge = getattr(specie, 'charge' if kind == 'm' else 'totcharge', None)
            charge_label = 'Oxidation state' if kind == 'm' else 'Charge'
            charge_value = ('+' if charge > 0 else '') + str(int(charge)) if charge is not None else 'Not resolved'
            if kind == 'm':
                title = str(getattr(specie, 'label', formula))
                if charge is not None:
                    title = f'[{title}{charge_value}]'
            else:
                title = getattr(specie, 'smiles', None) or formula

            unique_index = getattr(specie, 'unique_index', None)
            group_id = (kind, unique_index) if unique_index is not None else (kind, formula, title)
            if group_id not in entries:
                name = title
                num = 2
                while name in jmol_list_species:
                    name = title + ' (' + str(num) + ')'
                    num += 1
                entries[group_id] = name
                if kind == 'm':
                    desc = '<p>Metal center: ' + escape(str(getattr(specie, 'label', formula))) + '</p>'
                else:
                    rd = getattr(specie, 'rdkit_obj', None)
                    desc = None
                    if rd is not None:
                        try:
                            copy = Chem.Mol(rd)
                            Chem.rdDepictor.Compute2DCoords(copy)
                            drawer = rdMolDraw2D.MolDraw2DSVG(450, 320)
                            drawer.DrawMolecule(copy)
                            drawer.FinishDrawing()
                            desc = drawer.GetDrawingText().replace('svg:', '')
                            desc = re.sub(re__svghead, '', desc)
                            desc = re.sub(re__svgbackground, '', desc, 1)
                        except (ValueError, RuntimeError):
                            pass
                    if desc is None:
                        desc = '<p>Formula: ' + escape(str(formula)) + '</p>'
                compound_data.append((name, kind == 'm', desc, charge_label, charge_value))
                jmol_list_species[name] = ''

            name = entries[group_id]
            atoms = getattr(specie, 'atoms', None) or []
            if kind == 'm' and not atoms:
                atoms = [specie]
            coords = []
            for atom in atoms:
                coord = getattr(atom, 'coord', None)
                if coord is not None and len(coord) == 3:
                    coords.append('within (0.1, {{{} {} {}}})'.format(*coord))
            if coords:
                selector = ' or '.join(coords)
                if jmol_list_species[name]:
                    jmol_list_species[name] += ' or '
                jmol_list_species[name] += selector
    return compound_data, jmol_list_species


def bond_order_connectivity_moleculeset(moleclist):
    if not moleculeset_is_resolved(moleclist):
        return ' '
    return bond_order_connectivity(make_moleculeset_cell(moleclist))


def molxyz_to_string_xyz(mol, xyzCellVec, cmplut=None):

    celldesc = "{:f}, {:f}, {:f}, ".format(*tuple(xyzCellVec[0]))
    celldesc = celldesc + "{:f}, {:f}, {:f}, ".format(*tuple(xyzCellVec[1]))
    celldesc = celldesc + "{:f}, {:f}, {:f} ".format(*tuple(xyzCellVec[2]))

    allmols = []

    if mol.contains_metal:
        molName = "Complex"
    else:
        molName = "Other"

    atms = []
    for atoms in mol.atoms:                                                                                                  
        x,y,z = atoms.coord
        Z = atoms.label
        atms.append(" {:s}    {:8f} {:8f} {:8f}\\n".format(Z,x,y,z))
    allmols.append( "{:d}\\n {:s}\\n ".format(mol.natoms, molName) + "".join(atms) )

    return celldesc, "".join(allmols)

re__svghead = re.compile(r"<\?xml.*?\?>")
re__svgbackground = re.compile(r"<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='[0-9.]+' height='[0-9.]+' x='[0-9.]+' y='[0-9.]+'> </rect>")
rdMolDraw2D.MolDrawOptions.fixedScale = 1

def cell_to_svgs(cell, cmplut):
    """returns an svg image for every individual compound type in the cell"""
    res = []
    for name, lst in cmplut.items():
        tpl = lst[0]
        #for mol in cell.moleclist:
        mol = cell.unitcell.moleclist[tpl[0]]
        if tpl[1]=='m':
            mol = mol.metals[tpl[2]]
            sm = mol.label
            if mol.charge > 0:
                sm += f'+{mol.charge:d}'
            elif mol.charge < 0:
                sm += f'-{-mol.charge:d}'
            rd = Chem.MolFromSmiles('['+sm+']')
        elif tpl[1]=='l':
            mol = mol.ligands[tpl[2]]
            rd = mol.rdkit_obj
        elif tpl[1]=='t':
            rd = mol.rdkit_obj
        else:
            raise ValueError("internal error juggling compounds")
        try:
            Chem.rdDepictor.Compute2DCoords(rd)
        except Exception:
            raise ValueError(repr(rd) +' '+ sm)

        coords = rd.GetConformer(-1).GetPositions()
        bbox = (coords.max(axis=0) - coords.min(axis=0))[:2]
        size = (30*bbox+30).round()

        try:
            drawer = rdMolDraw2D.MolDraw2DSVG(int(size[0]), int(size[1]))
        except Exception as err:
            raise err
            raise ValueError(repr(size), repr(coords))
        drawer.DrawMolecule(rd)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        #svg = rdMolDraw2D.MolToSVG(rd)
        svg = svg.replace('svg:', '')
        svg = re.sub(re__svghead, '', svg)
        svg = re.sub(re__svgbackground, '', svg, 1)
        #res.append(name+'   '+svg)
        res.append(svg)


    res.append(repr(cell.cell_param))
    return res


def cell_to_svgs_base(refCell, cmplut):
    """
    Return an SVG image for every individual compound type in REFERENCE mode.

    Metals are drawn directly from the element symbol.
    Ligands and standalone molecules are drawn from one of their
    plausible reference-mode charge states.
    """
    res = []

    for name, lst in cmplut.items():

        tpl = lst[0]

        rd = None
        specie = None

        mol = refCell.refmoleclist[tpl[0]]

        # -------------------------------------------------------------
        # METAL
        # -------------------------------------------------------------
        if tpl[1] == 'm':

            specie = mol.metals[tpl[2]]

            label = specie.label
            smiles = '[' + label + ']'

            rd = Chem.MolFromSmiles(smiles)

            if rd is None:
                raise ValueError(
                    "Could not construct RDKit metal molecule: "
                    f"name={name!r}, "
                    f"tpl={tpl!r}, "
                    f"label={label!r}, "
                    f"smiles={smiles!r}, "
                    f"metal={specie!r}"
                )

        # -------------------------------------------------------------
        # LIGAND
        # -------------------------------------------------------------
        elif tpl[1] == 'l':

            specie = mol.ligands[tpl[2]]

            rd = _reference_rdkit_obj(specie)

        # -------------------------------------------------------------
        # STANDALONE MOLECULE
        # -------------------------------------------------------------
        elif tpl[1] == 't':

            specie = mol

            rd = _reference_rdkit_obj(specie)

        else:
            raise ValueError(
                "Internal error juggling reference compounds: "
                f"name={name!r}, tpl={tpl!r}"
            )

        if rd is None:
            raise ValueError(
                "RDKit object unexpectedly None after reference conversion: "
                f"name={name!r}, "
                f"tpl={tpl!r}, "
                f"type={tpl[1]!r}, "
                f"specie={specie!r}"
            )

        #try:
        #    Chem.rdDepictor.Compute2DCoords(rd)

        #except Exception as e:
        #    raise ValueError(
        #        "Compute2DCoords failed in reference mode: "
        #        f"name={name!r}, "
        #        f"tpl={tpl!r}, "
        #        f"type={tpl[1]!r}, "
        #        f"formula={getattr(specie, 'formula', None)!r}, "
        #        f"rd={rd!r}, "
        #        f"specie={specie!r}, "
        #        f"original_error={e!r}"
        #    ) from e
        try:
            Chem.rdDepictor.Compute2DCoords(rd)
        except Exception:
            raise ValueError(repr(rd) +' '+ sm)

        coords = rd.GetConformer(-1).GetPositions()
        bbox = (coords.max(axis=0) - coords.min(axis=0))[:2]
        size = (30*bbox+30).round()

        try:
            drawer = rdMolDraw2D.MolDraw2DSVG(int(size[0]), int(size[1]))
        except Exception as err:
            raise err
            raise ValueError(repr(size), repr(coords))
        drawer.DrawMolecule(rd)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        #svg = rdMolDraw2D.MolToSVG(rd)
        svg = svg.replace('svg:', '')
        svg = re.sub(re__svghead, '', svg)
        svg = re.sub(re__svgbackground, '', svg, 1)
        #res.append(name+'   '+svg)
        res.append(svg)

    res.append(repr(refCell.cell_param))
    return res


#def cell_to_svgs_base(refCell, cmplut):
#    """returns an svg image for every individual compound type in the cell"""
#    res = []
#    for name, lst in cmplut.items():
#        tpl = lst[0]
#        rd = None
#        sm = None
#
#        #for mol in cell.moleclist:
#        mol = refCell.refmoleclist[tpl[0]]
#        if tpl[1]=='m':
#            mol = mol.metals[tpl[2]]
#            sm = mol.label
#            #no charge in ref mode
#            #if mol.charge > 0:
#            #    sm += f'+{mol.charge:d}'
#            #elif mol.charge < 0:
#            #    sm += f'-{-mol.charge:d}'
#            smiles = '[' + sm + ']'
#            rd = Chem.MolFromSmiles(smiles)
#
#            if rd is None:
#                 raise ValueError(
#                        f"MolFromSmiles failed: "
#                        f"name={name!r}, tpl={tpl!r}, "
#                        f"label={sm!r}, smiles={smiles!r}, "
#                        f"mol={mol!r}"
#                        )
#
#        elif tpl[1]=='l':
#            mol = mol.ligands[tpl[2]]
#            rd = mol.rdkit_obj
#        elif tpl[1]=='t':
#            rd = mol.rdkit_obj
#        else:
#            raise ValueError("internal error juggling compounds")
#
#        if rd is None:
#                raise ValueError(
#                f"RDKit object is None: "
#                f"name={name!r}, "
#                f"tpl={tpl!r}, "
#                f"type={tpl[1]!r}, "
#                f"sm={sm!r}, "
#                f"mol={mol!r}"
#            )
#
#        try:
#            Chem.rdDepictor.Compute2DCoords(rd)
#        except Exception:
#            raise ValueError(repr(rd) +' '+ sm)
#
#        coords = rd.GetConformer(-1).GetPositions()
#        bbox = (coords.max(axis=0) - coords.min(axis=0))[:2]
#        size = (30*bbox+30).round()
#
#        try:
#            drawer = rdMolDraw2D.MolDraw2DSVG(int(size[0]), int(size[1]))
#        except Exception as err:
#            raise err
#            raise ValueError(repr(size), repr(coords))
#        drawer.DrawMolecule(rd)
#        drawer.FinishDrawing()
#        svg = drawer.GetDrawingText()
#        #svg = rdMolDraw2D.MolToSVG(rd)
#        svg = svg.replace('svg:', '')
#        svg = re.sub(re__svghead, '', svg)
#        svg = re.sub(re__svgbackground, '', svg, 1)
#        #res.append(name+'   '+svg)
#        res.append(svg)
#
#
#    res.append(repr(refCell.cell_param))
#    return res




def printing_text(cell, output):
    #dicts = {}
    #list_show = []
    #for idx, mol in enumerate(cell.moleclist):
    #    if mol.type == "Complex":        
    #        if mol.formula in dicts.keys():
    #            dicts[mol.formula] +=1 
    #        else:
    #            dicts[mol.formula] = 1
    #            list_show.append(idx)

    #for idx, mol in enumerate(cell.moleclist):
    #    if mol.type == "Other":     
    #        if mol.formula in dicts.keys():
    #            dicts[mol.formula] +=1 
    #        else:
    #            dicts[mol.formula] = 1
    #            list_show.append(idx)

    #for i in list_show:
    #    mol=cell.moleclist[i]
    #    if mol.type == "Complex": 
    #        output.extend([f"[Complex] Formula : {mol.formula}\t(occurrence : {dicts[mol.formula]})"])
    #        output.extend([f"   Total charge : {mol.totcharge}"])
    #        output.extend([f"   Spin : {mol.spin}"])
    #        output.extend([""])

    #        if mol.hapticity == False :
    #            for met in mol.metalist:
    #                output.extend([f"   >> Metal : {met.label}"])
    #                output.extend([f"   Metal oxidation state : {met.totcharge}"])
    #                output.extend([f"   Coordination number: {met.coordination_number}"])
    #                output.extend([f"   Metal-coordinating atoms: {met.coordinating_atoms}"])
    #                output.extend([f"   Coordination geometry: {met.geometry}\t(*deviation value : {met.deviation})"])
    #                output.extend(["   *deviation value: closer to 0, less distortion in a given geometry"])           
    #        else :
    #            for met in mol.metalist:
    #                output.extend([f"   >> Metal : {met.label}"])
    #                output.extend([f"   Metal oxidation state : {met.totcharge}"])
    #                output.extend([f"   Haptic ligand(s) bound to metal. Coordination number and geometry not shown"])
    #        output.extend([""])

    #        for lig in mol.ligandlist:
    #            output.extend([f"   >> Ligand Formula : {lig.formula}"])
    #            output.extend([f"   Charge : {lig.totcharge}"])
    #            if lig.hapticity == True :
    #                output.extend([f"   Hapticity: {lig.hapttype}"])
    #            else : 
    #                output.extend([f"   Denticity: {lig.totmconnec}"])
    #            output.extend([f"   Smiles: {lig.smiles}"])
    #            output.extend([""])

    #    elif mol.type == "Other" :
    #        output.extend([f"[Other] Formula : {mol.formula}\t(occurrence : {dicts[mol.formula]})"])
    #        output.extend([f"   Charge: {mol.totcharge}"])
    #        output.extend([f"   Smiles: {mol.smiles}"])
    #output = "OK"

    dicts = {}
    list_show = []
    for idx, mol in enumerate(cell.unitcell.moleclist):
        if mol.contains_metal:        
            if mol.formula in dicts.keys():
                dicts[mol.formula] +=1 
            else:
                dicts[mol.formula] = 1
                list_show.append(idx)

    for idx, mol in enumerate(cell.unitcell.moleclist):
        if mol.contains_metal == False:     
            if mol.formula in dicts.keys():
                dicts[mol.formula] +=1 
            else:
                dicts[mol.formula] = 1
                list_show.append(idx)

    for i in list_show:
        mol=cell.unitcell.moleclist[i]
        if mol.contains_metal: 
            output.extend([f"[Complex] Formula : {mol.formula}\t(occurrence : {dicts[mol.formula]}) \n"])
            output.extend([f"   Total charge : {mol.totcharge} \n"])
            output.extend([f"   Spin : {mol.spin} \n \n"])
            output.extend(["\n\n"])

            for metal in mol.metals:
                output.extend([f"   >> Metal : {metal.label}\n"])
                output.extend([f"   Oxidation state : {metal.charge}\n"])
                output.extend([f"   Spin : {metal.spin}\n"])
                output.extend([f"   Coordination number: {metal.coord_nr}\n"])
                output.extend([f"   Coordination geometry: {metal.coord_geometry}\n"])
                output.extend([f"   Coordination sphere formula: {metal.coord_sphere_formula}\n"])
                output.extend([f"   Relative metal radius: {metal.rel_metal_radius}\n"])
                output.extend(["\n\n"])
                #output.extend([f"   Valence electrons: {metal.get_valence_electrons(metal.charge)}"])

    
    #        #if mol.is_haptic == False :
    #        #    for met in mol.metalist:
    #        #        output.extend([f"   >> Metal : {met.label}"])
    #        #        output.extend([f"   Metal oxidation state : {met.totcharge}"])
    #        #        output.extend([f"   Coordination number: {met.coordination_number}"])
    #        #        output.extend([f"   Metal-coordinating atoms: {met.coordinating_atoms}"])
    #        #        output.extend([f"   Coordination geometry: {met.geometry}\t(*deviation value : {met.deviation})"])
    #        #        output.extend(["   *deviation value: closer to 0, less distortion in a given geometry"])           
    #        #else :
    #        #    for met in mol.metalist:
    #        #        output.extend([f"   >> Metal : {met.label}"])
    #        #        output.extend([f"   Metal oxidation state : {met.totcharge}"])
    #        #        #output.extend([f"   Haptic type: mol.haptic_type "])

            for lig in mol.ligands:
                output.extend([f"   >> Ligand Formula : {lig.formula} \n"])
                output.extend([f"   Charge : {lig.totcharge} \n"])
                if lig.is_haptic :
                    output.extend([f"   Hapticity: {lig.haptic_type} \n"])
                else : 
                    output.extend([f"   Denticity: {lig.denticity} \n"])
                output.extend([f"   Smiles: {lig.smiles} \n"])
                output.extend(["\n\n"])

        else :
            output.extend([f"[Other] Formula : {mol.formula}\t(occurrence : {dicts[mol.formula]})\n"])
            output.extend([f"   Charge: {mol.totcharge}\n"])
            output.extend([f"   Smiles: {mol.smiles}\n"])
            output.extend(["\n\n"])

    return output


def printing_text_refMol(refcell, output):

    dicts = {}
    list_show = []
    for idx, mol in enumerate(refcell.refmoleclist):
        if mol.contains_metal:        
            if mol.formula in dicts.keys():
                dicts[mol.formula] +=1 
            else:
                dicts[mol.formula] = 1
                list_show.append(idx)

    for idx, mol in enumerate(refcell.refmoleclist):
        if mol.contains_metal == False:     
            if mol.formula in dicts.keys():
                dicts[mol.formula] +=1 
            else:
                dicts[mol.formula] = 1
                list_show.append(idx)

    for i in list_show:
        mol=refcell.refmoleclist[i]
        if mol.contains_metal: 
            output.extend([f"[Complex] Formula : {mol.formula}\t(occurrence : {dicts[mol.formula]}) \n"])
            #output.extend([f"   Total charge : {mol.totcharge} \n"]) #No charge
            #output.extend([f"   Spin : {mol.spin} \n \n"]) #No spin
            output.extend(["\n\n"])

            for metal in mol.metals:
                output.extend([f"   >> Metal : {metal.label}\n"])
                #output.extend([f"   Oxidation state : {metal.charge}\n"])
                #output.extend([f"   Spin : {metal.spin}\n"])
                output.extend([f"   Coordination number: {metal.coord_nr}\n"])
                output.extend([f"   Coordination geometry: {metal.coord_geometry}\n"])
                output.extend([f"   Coordination sphere formula: {metal.coord_sphere_formula}\n"])
                output.extend([f"   Relative metal radius: {metal.rel_metal_radius}\n"])
                output.extend(["\n\n"])

            for lig in mol.ligands:
                output.extend([f"   >> Ligand Formula : {lig.formula} \n"])
                #output.extend([f"   Charge : {lig.totcharge} \n"]) #no charge
                if lig.is_haptic :
                    output.extend([f"   Hapticity: {lig.haptic_type} \n"])
                else : 
                    output.extend([f"   Denticity: {lig.denticity} \n"])
                #output.extend([f"   Smiles: {lig.smiles} \n"]) #no smiles available
                output.extend(["\n\n"])

        else :
            output.extend([f"[Other] Formula : {mol.formula}\t(occurrence : {dicts[mol.formula]})\n"])
            #output.extend([f"   Charge: {mol.totcharge}\n"]) #no charge
            #output.extend([f"   Smiles: {mol.smiles}\n"]) #no smiles
            output.extend(["\n\n"])


    return output


def printing_text_molxyz(mol, output):

    if mol.contains_metal: 
        output.extend([f"[Complex] Formula : {mol.formula}\t \n"])
        #output.extend([f"   Total charge : {mol.totcharge} \n"]) #No charge
        #output.extend([f"   Spin : {mol.spin} \n \n"]) #No spin
        output.extend(["\n\n"])

        for metal in mol.metals:
            output.extend([f"   >> Metal : {metal.label}\n"])
            #output.extend([f"   Oxidation state : {metal.charge}\n"])
            #output.extend([f"   Spin : {metal.spin}\n"])
            output.extend([f"   Coordination number: {metal.coord_nr}\n"])
            output.extend([f"   Coordination geometry: {metal.coord_geometry}\n"])
            output.extend([f"   Coordination sphere formula: {metal.coord_sphere_formula}\n"])
            output.extend([f"   Relative metal radius: {metal.rel_metal_radius}\n"])
            output.extend(["\n\n"])

        for lig in mol.ligands:
            output.extend([f"   >> Ligand Formula : {lig.formula} \n"])
            #output.extend([f"   Charge : {lig.totcharge} \n"]) #no charge
            if lig.is_haptic :
                output.extend([f"   Hapticity: {lig.haptic_type} \n"])
            else : 
                output.extend([f"   Denticity: {lig.denticity} \n"])
            #output.extend([f"   Smiles: {lig.smiles} \n"]) #no smiles available
            output.extend(["\n\n"])

    else :
        output.extend([f"[Other] Formula : {mol.formula}\t \n"])
        #output.extend([f"   Charge: {mol.totcharge}\n"]) #no charge
        #output.extend([f"   Smiles: {mol.smiles}\n"]) #no smiles
        output.extend(["\n\n"])


    return output




def molecules_list(cell):                                                                                      
    ''' Takes the cell2mol cell object and uses its information to make a list of all the molecules with its respectives
    atomic coordinates in a compatible format for jsmol

    Args:
        cell: the output of cell2mol

    Return:
        A list containing all the molecules separeted and its respectives atoms coordinates in a string for jsmol
    '''

    totmol = len(cell.unitcell.moleclist)                                                                                          
    jmol_list_pos = {}                                                                                                    
    for idx, mol in enumerate(cell.unitcell.moleclist):                                                                                            

        if mol.contains_metal:
            molName = mol.get_parent("reference").name + "_Complex_"+str(idx) 
        else:
            molName = mol.get_parent("reference").name + "_Other_"+ str(idx)

        jmol_list_pos[molName] = " select "                                                                              
        cont=0                                                                                                            
        for a in mol.atoms:                                                                                               
            jmol_list_pos[molName] = jmol_list_pos[molName] + " within " +"(0.1, {" + str(a.coord[0]) + " " + str(a.coord[1]) + " " + str(a.coord[2]) + "})"
            cont=cont+1                                                                                                   
            if (cont < mol.natoms):                                                                                       
                jmol_list_pos[molName] = jmol_list_pos[molName] + " or "                                                

    return jmol_list_pos


def molecules_list_reference(cell):                                                                                      
    ''' Takes the cell2mol cell object and uses its information to make a list of all the molecules with its respectives
    atomic coordinates in a compatible format for jsmol

    Args:
        cell: the output of cell2mol

    Return:
        A list containing all the molecules separeted and its respectives atoms coordinates in a string for jsmol
    '''

    totmol = len(cell.refmoleclist)                                                                                          
    jmol_list_pos = {}                                                                                                    
    for idx, mol in enumerate(cell.refmoleclist):                                                                                            

        if mol.contains_metal:
            molName = mol.get_parent("reference").name + "_Complex_"+str(idx) 
        else:
            molName = mol.get_parent("reference").name + "_Other_"+ str(idx)

        jmol_list_pos[molName] = " select "                                                                              
        cont=0                                                                                                            
        for a in mol.atoms:                                                                                               
            jmol_list_pos[molName] = jmol_list_pos[molName] + " within " +"(0.1, {" + str(a.coord[0]) + " " + str(a.coord[1]) + " " + str(a.coord[2]) + "})"
            cont=cont+1                                                                                                   
            if (cont < mol.natoms):                                                                                       
                jmol_list_pos[molName] = jmol_list_pos[molName] + " or "                                                

    return jmol_list_pos


def _molecule_bond_orders_to_jsmol(moleclist):
    jmolCon = " "
    for mol in moleclist or []:
        for bond in getattr(mol, "bonds", None) or []:
            if isinstance(bond, dict):
                order = bond.get("order", 1.0)
                atom_indices = bond.get("atom_indices", [])
                if order <= 1.0 or len(atom_indices) != 2:
                    continue
                atom1 = mol.atoms[atom_indices[0]]
                atom2 = mol.atoms[atom_indices[1]]
            else:
                order = bond.order
                if order <= 1.0:
                    continue
                atom1 = bond.atom1
                atom2 = bond.atom2
            jmolCon += " select within (0.1, {" + " ".join(map(str, atom1.coord)) + "}) or "
            jmolCon += " within (0.1, {" + " ".join(map(str, atom2.coord)) + "}) ; "
            jmolCon += " bondOrder " + str(order) + " ; "
    return jmolCon


def bond_order_connectivity(cell):
    return _molecule_bond_orders_to_jsmol(cell.unitcell.moleclist)


def bond_order_connectivity_reference(refCell):
    return _molecule_bond_orders_to_jsmol(refCell.refmoleclist)


def species_list(cell):
    ''' Takes the cell2mol cell object and uses its information to make a list of all the species(metals, ligands, and others) 
    with its respectives atomic coordinates in a compatible format for jsmol

    Args:
        cell: the output of cell2mol

    Return:
        A dictionary containing all the speciess separeted and its respectives atoms coordinates in a string for jsmol
    '''

    jmol_list_species = {}
    for mol in cell.unitcell.moleclist:
        if mol.contains_metal :
            for ligand in mol.ligands:
                if ligand.smiles not in jmol_list_species:
                    jmol_list_species[ligand.smiles] = " "
                else:
                    jmol_list_species[ligand.smiles] += " or "
                for nat, atms in enumerate(ligand.atoms):
                    jmol_list_species[ligand.smiles] = jmol_list_species[ligand.smiles] + " within (0.1, {" + str(atms.coord[0]) + " "
                    jmol_list_species[ligand.smiles] = jmol_list_species[ligand.smiles] + str(atms.coord[1]) + " "
                    jmol_list_species[ligand.smiles] = jmol_list_species[ligand.smiles] + str(atms.coord[2]) + "})"
                    if nat+1 < len(ligand.atoms):
                        jmol_list_species[ligand.smiles] = jmol_list_species[ligand.smiles] + " or "
            for metal in mol.metals:
                if metal.charge > 0:
                    metalName = f'[{metal.label:s}+{metal.charge:d}]'
                elif metal.charge < 0:  #elif mtl.totcharge < 0: 
                    metalName = f'[{metal.label:s}-{metal.charge:d}]'
                else:
                    metalName = f'{metal.label:s}'
                if metalName not in jmol_list_species:
                    jmol_list_species[metalName] = " "
                else:
                    jmol_list_species[metalName] += " or "
                #for nat in range(metal.natom):
                jmol_list_species[metalName] += " within (0.1, {" + str(metal.coord[0]) + " "
                jmol_list_species[metalName] += str(metal.coord[1]) + " "
                jmol_list_species[metalName] += str(metal.coord[2]) + "}) "
                #    if nat+1 < metal.natom:
                #        jmol_list_species[metalName] += " or "
        else:
            if mol.smiles not in jmol_list_species:
                jmol_list_species[mol.smiles] = " "
            else:
                jmol_list_species[mol.smiles] += " or "
            for nat, atms in enumerate(mol.atoms):
                jmol_list_species[mol.smiles] = jmol_list_species[mol.smiles] + " within (0.1, {" + str(atms.coord[0]) + " "
                jmol_list_species[mol.smiles] = jmol_list_species[mol.smiles] + str(atms.coord[1]) + " "
                jmol_list_species[mol.smiles] = jmol_list_species[mol.smiles] + str(atms.coord[2]) + "})"
                if nat+1 < len(mol.atoms):
                    jmol_list_species[mol.smiles] = jmol_list_species[mol.smiles] + " or "
    
    return jmol_list_species


def species_list_reference(refCell):
    ''' Takes the cell2mol cell object and uses its information to make a list of all the species(metals, ligands, and others) 
    with its respectives atomic coordinates in a compatible format for jsmol

    Args:
        cell: the output of cell2mol

    Return:
        A dictionary containing all the speciess separeted and its respectives atoms coordinates in a string for jsmol
    '''

    jmol_list_species = {}
    for mol in refCell.refmoleclist:
        if mol.contains_metal :
            for ligand in mol.ligands:
                if ligand.plausible_charge_states[0].species_smiles not in jmol_list_species:
                    jmol_list_species[ligand.plausible_charge_states[0].species_smiles] = " "
                else:
                    jmol_list_species[ligand.plausible_charge_states[0].species_smiles] += " or "
                for nat, atms in enumerate(ligand.atoms):
                    jmol_list_species[ligand.plausible_charge_states[0].species_smiles] = jmol_list_species[ligand.plausible_charge_states[0].species_smiles] + " within (0.1, {" + str(atms.coord[0]) + " "
                    jmol_list_species[ligand.plausible_charge_states[0].species_smiles] = jmol_list_species[ligand.plausible_charge_states[0].species_smiles] + str(atms.coord[1]) + " "
                    jmol_list_species[ligand.plausible_charge_states[0].species_smiles] = jmol_list_species[ligand.plausible_charge_states[0].species_smiles] + str(atms.coord[2]) + "})"
                    if nat+1 < len(ligand.atoms):
                        jmol_list_species[ligand.plausible_charge_states[0].species_smiles] = jmol_list_species[ligand.plausible_charge_states[0].species_smiles] + " or "
            for metal in mol.metals:
                #no charge in reference
                #if metal.charge > 0:
                #    metalName = f'[{metal.label:s}+{metal.charge:d}]'
                #elif metal.charge < 0:  #elif mtl.totcharge < 0: 
                #    metalName = f'[{metal.label:s}-{metal.charge:d}]'
                #else:
                metalName = f'{metal.label:s}'
                if metalName not in jmol_list_species:
                    jmol_list_species[metalName] = " "
                else:
                    jmol_list_species[metalName] += " or "
                #for nat in range(metal.natom):
                jmol_list_species[metalName] += " within (0.1, {" + str(metal.coord[0]) + " "
                jmol_list_species[metalName] += str(metal.coord[1]) + " "
                jmol_list_species[metalName] += str(metal.coord[2]) + "}) "
                #    if nat+1 < metal.natom:
                #        jmol_list_species[metalName] += " or "
        else:
            if mol.smiles not in jmol_list_species:
                jmol_list_species[mol.smiles] = " "
            else:
                jmol_list_species[mol.smiles] += " or "
            for nat, atms in enumerate(mol.atoms):
                jmol_list_species[mol.smiles] = jmol_list_species[mol.smiles] + " within (0.1, {" + str(atms.coord[0]) + " "
                jmol_list_species[mol.smiles] = jmol_list_species[mol.smiles] + str(atms.coord[1]) + " "
                jmol_list_species[mol.smiles] = jmol_list_species[mol.smiles] + str(atms.coord[2]) + "})"
                if nat+1 < len(mol.atoms):
                    jmol_list_species[mol.smiles] = jmol_list_species[mol.smiles] + " or "
    
    return jmol_list_species


