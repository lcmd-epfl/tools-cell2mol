import sys, os, io
sys.path.append(os.path.join(os.path.split(__file__)[0], '../cell2mol'))
import cell2mol
from cell2mol.cif2info import cif_2_info
from cell2mol.c2m_module import save_cell, cell2mol

#import ase
#import nglview as ngl

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = io.StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout


def view_cell(cell):
    xsf = "CRYSTAL\nPRIMVEC\n"
    for vec in cell.cellvec:
        xsf += "{:f} {:f} {:f}\n".format(vec[0],vec[1],vec[2])
    nat = sum(mol.natoms for mol in cell.moleclist)
    xsf += "PRIMCOORD\n{:d} 1\n".format(nat)
    for mol in cell.moleclist:
        for Z,(x,y,z) in zip(mol.atnums, mol.coord):
            xsf += "{:d} {:f} {:f} {:f}\n".format(Z,x,y,z)
    return xsf





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
