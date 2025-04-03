import sys, os, io, re
sys.path.append(os.path.join(os.path.split(__file__)[0], '../cell2mol'))

import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import cell2mol
from cell2mol.refcell import process_refcell
from cell2mol.unitcell import process_unitcell
from cell2mol.xyz_molecule import get_molecule
from cell2mol.read_write import print_molecule

#from cell2mol.read_write import savemolecules
#from cell2mol.readwrite import readinfo, savemolecules

#from cell2mol.final_c2m_driver import handle_cif_file
#from cell2mol.cif2info import cif_2_info
#from cell2mol.c2m_module import save_cell, cell2mol

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

        if mol.iscomplex:
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
    for i_mol,mol in enumerate(cell.moleclist):
        #if mol.type == 'Complex':
        if mol.iscomplex :
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

def cell_get_metal_desc(cell, cmplut):
    res = []
    for name, lst in cmplut.items():
        tpl = lst[0]
    #for mol in cell.moleclist:
        mol = cell.moleclist[tpl[0]]
        if tpl[1]=='m':
            mtl = mol.metals[tpl[2]]
            res.append('<p>Metal center: {0:s}<br/>predicted charge: {1:+d}</p>'.format(
                mtl.label,
                mtl.charge,
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
    celldesc = "a={:f},b={:f},c={:f},alpha={:f},beta={:f},gamma={:f}".format(*tuple(cell.cell_param))

    allmols = []
    for idx, mols in enumerate(cell.moleclist):                                                                                                   

        if mols.iscomplex:
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
    celldesc = "a={:f},b={:f},c={:f},alpha={:f},beta={:f},gamma={:f}".format(*tuple(cell.cell_param))

    allmols = []
    for idx, mols in enumerate(cell.refmoleclist):                                                                                                   

        if mols.iscomplex:
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


re__svghead = re.compile(r"<\?xml.*?\?>")
re__svgbackground = re.compile(r"<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='[0-9.]+' height='[0-9.]+' x='[0-9.]+' y='[0-9.]+'> </rect>")
rdMolDraw2D.MolDrawOptions.fixedScale = 1

def cell_to_svgs(cell, cmplut):
    """returns an svg image for every individual compound type in the cell"""
    res = []
    for name, lst in cmplut.items():
        tpl = lst[0]
        #for mol in cell.moleclist:
        mol = cell.moleclist[tpl[0]]
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
    for idx, mol in enumerate(cell.moleclist):
        if mol.iscomplex:        
            if mol.formula in dicts.keys():
                dicts[mol.formula] +=1 
            else:
                dicts[mol.formula] = 1
                list_show.append(idx)

    for idx, mol in enumerate(cell.moleclist):
        if mol.iscomplex == False:     
            if mol.formula in dicts.keys():
                dicts[mol.formula] +=1 
            else:
                dicts[mol.formula] = 1
                list_show.append(idx)

    for i in list_show:
        mol=cell.moleclist[i]
        if mol.iscomplex: 
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
        if mol.iscomplex:        
            if mol.formula in dicts.keys():
                dicts[mol.formula] +=1 
            else:
                dicts[mol.formula] = 1
                list_show.append(idx)

    for idx, mol in enumerate(refcell.refmoleclist):
        if mol.iscomplex == False:     
            if mol.formula in dicts.keys():
                dicts[mol.formula] +=1 
            else:
                dicts[mol.formula] = 1
                list_show.append(idx)

    for i in list_show:
        mol=refcell.refmoleclist[i]
        if mol.iscomplex: 
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






def molecules_list(cell):                                                                                      
    ''' Takes the cell2mol cell object and uses its information to make a list of all the molecules with its respectives
    atomic coordinates in a compatible format for jsmol

    Args:
        cell: the output of cell2mol

    Return:
        A list containing all the molecules separeted and its respectives atoms coordinates in a string for jsmol
    '''

    totmol = len(cell.moleclist)                                                                                          
    jmol_list_pos = {}                                                                                                    
    for idx, mol in enumerate(cell.moleclist):                                                                                            

        if mol.iscomplex:
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

        if mol.iscomplex:
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


def bond_order_connectivity(cell):
    ''' Takes the cell object and returns a string containig the information needed by jsmol to generate the atomic bonds with 
    the correct bond order

    Args:
        cell: the output of cell2mol
        

    Return:
        A string with the jsmol instructions to select all bonded atoms within the same molecule and specify its bond order
    '''

    jmolCon = " "

    #for mol in cell.moleclist: #loop over all molecules
    #    connectMat = mol.conmat #connectivity matrix
    #    a = mol.atoms #list with all the atoms forming the molecule
    #
    #    for atomi, atomCon in enumerate(connectMat): #loop over atoms 
    #        for atomj, conn in enumerate(atomCon): #loop over atoms
    #            if (conn == 1. and atomj > atomi): #check if there is a bond and avoid double counting
    #                atomiBonds = a[atomi].bond #all bonds of atomi
    #                for bonds in atomiBonds: #loop over all bonds (the first index is always < than the second)
    #                    if (atomi == bonds[0] and atomj == bonds[1]): #if the bonded atoms matches
    #                        #select atomi by its coordinates
    #                        jmolCon = jmolCon + " select within (0.1, {" #+ str(atomi) + " " +str(atomCon)
    #                        jmolCon = jmolCon + str(a[int(atomi)].coord[0]) + " "
    #                        jmolCon = jmolCon + str(a[int(atomi)].coord[1]) + " "
    #                        jmolCon = jmolCon + str(a[int(atomi)].coord[2]) + " "
    #                        jmolCon = jmolCon + "}) or "
    #                        #select atomj by its coordinates
    #                        jmolCon = jmolCon + " within (0.1, {" #+ str(atomi) + " " +str(atomCon)
    #                        jmolCon = jmolCon + str(a[int(atomj)].coord[0]) + " "
    #                        jmolCon = jmolCon + str(a[int(atomj)].coord[1]) + " "
    #                        jmolCon = jmolCon + str(a[int(atomj)].coord[2]) + " "
    #                        jmolCon = jmolCon + "}) ;"
    #                        #bond order
    #                        jmolCon = jmolCon + " bondorder " + str(bonds[2]) + " ; "



    #for mol in cell.moleclist: #loop over all molecules
    #    
    #    #connectMat = cell.moleclist[0].get_adjmatrix()[0] #connectivity matrix
    #    connectMat = mol.get_adjmatrix()[0] #connectivity matrix
    #    atm = mol.atoms
    #        
    #    for atomi, atomCon in enumerate(connectMat): #loop over atoms
    #        for atomj, conn in enumerate(atomCon[atomi:]): #loop over connectivity of each not repeated atom of previous loop
    #            if (conn == 1): #if there is a bond
    #                #select atomi by its coordinates
    #                jmolCon = jmolCon + " select within (0.1, {" #+ str(atomi) + " " +str(atomCon)
    #                jmolCon = jmolCon + str(atm[atomi].coord[0]) + " "
    #                jmolCon = jmolCon + str(atm[atomi].coord[1]) + " "
    #                jmolCon = jmolCon + str(atm[atomi].coord[2]) + " "
    #                jmolCon = jmolCon + "}) or "
    #                #select atomj by its coordinates
    #                jmolCon = jmolCon + " within (0.1, {" #+ str(atomi) + " " +str(atomCon)
    #                jmolCon = jmolCon + str(atm[atomj].coord[0]) + " "
    #                jmolCon = jmolCon + str(atm[atomj].coord[1]) + " "
    #                jmolCon = jmolCon + str(atm[atomj].coord[2]) + " "
    #                jmolCon = jmolCon + "}) ;"
    #                for bond in atm[atomi].bonds: #loop ovre all bonds of atomi
    #                    if (bond.atom2.coord[0] == atm[atomj].coord[0] and
    #                        bond.atom2.coord[1] == atm[atomj].coord[1] and
    #                        bond.atom2.coord[2] == atm[atomj].coord[2]): #if the bond is with atomj

    #                        #jmolCon = str( bond.atom2.coord.all() == atm[atomj].coord.all() ) #if the bond is with atomj
    #                        #jmolCon = str(bond.order)
    #                        #bond order
    #                        jmolCon = jmolCon + " bondorder " + str(bond.order) + " ; "
    #                        break


    #Double looping. atom1-atom2 = atom2-atom1. Can be improved            
    for mol in cell.moleclist:
        for atm in mol.atoms:
            for bond in atm.bonds: #loop over all atoms
                if (bond.order > 1.0):
                    #select atom 1 by its coordinates
                    jmolCon = jmolCon + " select within (0.1, {" #+ str(atomi) + " " +str(atomCon)
                    jmolCon = jmolCon + str(bond.atom1.coord[0]) + " "
                    jmolCon = jmolCon + str(bond.atom1.coord[1]) + " "
                    jmolCon = jmolCon + str(bond.atom1.coord[2]) + " "
                    jmolCon = jmolCon + "}) or "
                    #select atom 2 by its coordinates
                    jmolCon = jmolCon + " within (0.1, {" #+ str(atomi) + " " +str(atomCon)
                    jmolCon = jmolCon + str(bond.atom2.coord[0]) + " "
                    jmolCon = jmolCon + str(bond.atom2.coord[1]) + " "
                    jmolCon = jmolCon + str(bond.atom2.coord[2]) + " "
                    jmolCon = jmolCon + "}) ; "
                    #bond order
                    if (bond.order == 2.0):
                        jmolCon = jmolCon + " bondOrder 2 ; "
                    elif (bond.order == 3.0):
                        jmolCon = jmolCon + " bondOrder 3 ; "
                    else:
                        jmolCon = jmolCon + " bondOrder " + str(bond.order) + " ; "








    
    return jmolCon


def species_list(cell):
    ''' Takes the cell2mol cell object and uses its information to make a list of all the species(metals, ligands, and others) 
    with its respectives atomic coordinates in a compatible format for jsmol

    Args:
        cell: the output of cell2mol

    Return:
        A dictionary containing all the speciess separeted and its respectives atoms coordinates in a string for jsmol
    '''

    jmol_list_species = {}
    for mol in cell.moleclist:
        if mol.iscomplex :
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


