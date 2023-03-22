import sys, os, io, re
sys.path.append(os.path.join(os.path.split(__file__)[0], '../cell2mol'))

import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

import cell2mol
from cell2mol.cif2info import cif_2_info
from cell2mol.c2m_module import save_cell, cell2mol
from cell2mol.readwrite import readinfo, savemolecules



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
        if mol.type == 'Complex':
            for i_mtl, mtl in enumerate(mol.metalist):
                mlabel = mtl.label
                if mtl.totcharge > 0:
                    mlabel = f'[{mlabel:s}+{mtl.totcharge:d}]'
                elif mtl.totcharge < 0:
                    mlabel = f'[{mlabel:s}-{-mtl.totcharge:d}]'
                names.setdefault(mlabel, [])
                names[mlabel].append((i_mol,'m',i_mtl))
            for i_lig, lig in enumerate(mol.ligandlist):
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
            mtl = mol.metalist[tpl[2]]
            res.append('<p>Metal center: {0:s}<br/>predicted charge: {1:+d}</p>'.format(
                mtl.label,
                mtl.totcharge,
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

def cell_to_string_xyz(cell, cmplut=None):
    celldesc = "a={:f},b={:f},c={:f},alpha={:f},beta={:f},gamma={:f}".format(*tuple(cell.cellparam))
    mols = []
    # for mol in cell.moleclist:
    #     atms = []
    #     for Z,(x,y,z) in zip(mol.atnums, mol.coord):
    #         atms.append(" {:s}    {:8f} {:8f} {:8f}\n".format(ELEMENTS[Z],x,y,z))
    #     mols.append(
    #         "{:d}\n{:s}\n".format(len(atms), mol.name) +
    #         "".join(atms)
    #     )
    for name, elems in cmplut.items():
        atms = []
        for elem in elems:
            mol = cell.moleclist[elem[0]]
            if elem[1] == 't':
                pass
            elif elem[1] == 'l':
                mol = mol.ligandlist[elem[2]]
            elif elem[1] == 'm':
                mol = mol.metalist[elem[2]]
                x,y,z = mol.atom.coord
                Z = mol.atom.atnum
                atms.append(" {:s}    {:8f} {:8f} {:8f}\n".format(ELEMENTS[Z],x,y,z))
                continue

            for Z,(x,y,z) in zip(mol.atnums, mol.coord):
                atms.append(" {:s}    {:8f} {:8f} {:8f}\n".format(ELEMENTS[Z],x,y,z))
        mols.append(
            "{:d}\n{:s}\n".format(len(atms), mol.name) +
            "".join(atms)
        )
    return celldesc, "".join(mols)


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
            mol = mol.metalist[tpl[2]]
            sm = mol.label
            if mol.totcharge > 0:
                sm += f'+{mol.totcharge:d}'
            elif mol.totcharge < 0:
                sm += f'-{-mol.totcharge:d}'
            rd = Chem.MolFromSmiles('['+sm+']')
        elif tpl[1]=='l':
            mol = mol.ligandlist[tpl[2]]
            rd = mol.rdkit_mol
        elif tpl[1]=='t':
            rd = mol.rdkit_mol
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


    res.append(repr(cell.cellparam))
    return res




def printing_text(cell, output):
    dicts = {}
    list_show = []
    for idx, mol in enumerate(cell.moleclist):
        if mol.type == "Complex":        
            if mol.formula in dicts.keys():
                dicts[mol.formula] +=1 
            else:
                dicts[mol.formula] = 1
                list_show.append(idx)

    for idx, mol in enumerate(cell.moleclist):
        if mol.type == "Other":     
            if mol.formula in dicts.keys():
                dicts[mol.formula] +=1 
            else:
                dicts[mol.formula] = 1
                list_show.append(idx)

    for i in list_show:
        mol=cell.moleclist[i]
        if mol.type == "Complex": 
            output.extend([f"[Complex] Formula : {mol.formula}\t(occurrence : {dicts[mol.formula]})"])
            output.extend([f"   total charge : {mol.totcharge}"])
            output.extend([f"   spin : {mol.spin}"])
            output.extend([""])

            if mol.hapticity == False :
                for met in mol.metalist:
                    output.extend([f"   >> Metal : {met.label}"])
                    output.extend([f"   Metal oxidation state : {met.totcharge}"])
                    output.extend([f"   coordination number: {met.coordination_number}"])
                    output.extend([f"   metal-coordinating atoms: {met.coordinating_atoms}"])
                    output.extend([f"   coordination geometry: {met.geometry}\t(*deviation value : {met.deviation})"])
                    output.extend(["   *deviation value: closer to 0, less distortion in a given geometry"])           
            else :
                for met in mol.metalist:
                    output.extend([f"   >> Metal : {met.label}"])
                    output.extend([f"   Metal oxidation state : {met.totcharge}"])
                    output.extend([f"   Haptic ligand(s) bound to metal. Coordination number and geometry not shown"])
            output.extend([""])

            for lig in mol.ligandlist:
                output.extend([f"   >> Ligand Formula : {lig.formula}"])
                output.extend([f"   charge : {lig.totcharge}"])
                output.extend([f"   spin : {lig.spin}"])
                if lig.hapticity == True :
                    output.extend([f"   hapticity: {lig.hapttype}"])
                else : 
                    output.extend([f"   denticity: {lig.totmconnec}"])
                output.extend([f"   smiles: {lig.smiles}"])
                output.extend([""])

        elif mol.type == "Other" :
            output.extend([f"[Other] Formula : {mol.formula}\t(occurrence : {dicts[mol.formula]})"])
            output.extend([f"   charge: {mol.totcharge}"])
            output.extend([f"   smiles: {mol.smiles}"])

    return output


