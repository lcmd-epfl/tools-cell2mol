import traceback
import pickle
from tools_barebone.structure_importers import get_structure_tuple, UnknownFormatError
import io
import flask
from .interface import *
from .tokens import monitoring, Token, InvalidXYZInput
from .calculation_logging import (
    calculation_logged, start_calculation_log,
    mark_calculation_failure, get_calculation_log_preview,
    calculation_log_path,
)
from .error_reporting import (
    describe_processing_error, molecule_status_message,
    Cell2molNoResultError,
)
#from cell2mol.c2m_module import save_cell
#from cell2mol.read_write import savemolecules, writexyz

#from cell2mol.unitcell import process_unitcell
#from cell2mol.refcell import process_refcell
#from cell2mol.xyz_molecule import get_molecule

from cell2mol.process_unitcell import interpret_unitcell 
from cell2mol.process_reference import interpret_reference 

import os
import glob

import webbrowser

import logging

from flask import send_file
from werkzeug.datastructures import FileStorage

blueprint = flask.Blueprint("compute", __name__, url_prefix="/compute")
    
def _safe_save_molecule_json(molset, path):
    molset.save(path, format="json")


def _cell2mol_reference_interpretation(token, refcell, unitcell=None):
    summary_path = os.path.join(token.get_path(), "reference_summary.txt")
    try:
        with open(summary_path, encoding="utf-8") as summary_file:
            summary_lines = summary_file.read().splitlines()
        if summary_lines:
            for index, line in enumerate(summary_lines):
                if line.startswith((
                    "Possible charges of species in ",
                    "No species list found in the cell object.",
                    "Reference Error (mode=",
                    "Potential issues and data status detected:",
                )):
                    summary_lines = summary_lines[:index]
                    break
            while summary_lines and not summary_lines[-1].strip():
                summary_lines.pop()
            if summary_lines:
                return summary_lines
    except (OSError, UnicodeError):
        logging.getLogger(__name__).warning(
            "Could not read the cell2mol reference summary: %s", summary_path,
            exc_info=True,
        )

    from cell2mol.write_results import write_cell_molecule_details, write_unique_species
    summary = io.StringIO()
    print(token.refcode, file=summary)
    print(
        f"Total charge comparison result: {getattr(refcell, 'total_charge_comparison', None)}",
        file=summary,
    )
    print(
        f"Metal oxidation state comparison result: {getattr(refcell, 'metal_os_comparison', None)}",
        file=summary,
    )
    print(f"  - Reported:   {getattr(refcell, 'reported_metal_os', None)}", file=summary)
    print(f"  - Matched:    {getattr(refcell, 'reported_metal_os_matched', None)}", file=summary)
    confidence = getattr(refcell, "metal_os_match_confidence", None)
    print(f"  - Confidence: {confidence:.2f}" if confidence is not None else "  - Confidence: N/A", file=summary)
    if unitcell is not None:
        method = getattr(unitcell, "os_method", None)
        print(
            f"Metal oxidation states determined by: {' -> '.join(method) if method else 'N/A'}",
            file=summary,
        )
    write_cell_molecule_details(refcell, file=summary)
    write_unique_species(refcell, file=summary)
    return summary.getvalue().splitlines()


def _cell2mol_molecule_interpretation(token, molset):
    summary_path = os.path.join(token.get_path(), "molecule_summary.txt")
    try:
        with open(summary_path, encoding="utf-8") as summary_file:
            summary_lines = summary_file.read().splitlines()
        for index, line in enumerate(summary_lines):
            if line.startswith((
                "Possible charges of species in ",
                "No species list found in the cell object.",
                "Molecule Error (mode=",
                "Potential issues and data status detected:",
            )):
                summary_lines = summary_lines[:index]
                break
        while summary_lines and not summary_lines[-1].strip():
            summary_lines.pop()
        if summary_lines:
            return summary_lines
        logging.getLogger(__name__).warning("Empty cell2mol molecule summary: %s", summary_path)
    except (OSError, UnicodeError):
        logging.getLogger(__name__).warning(
            "Could not read the cell2mol molecule summary: %s", summary_path,
            exc_info=True,
        )

    from cell2mol.write_results import write_molecule_details, write_unique_species
    summary = io.StringIO()
    print(token.refcode, file=summary)
    print(f"{len(molset.moleclist)} molecules, totcharge={molset.totcharge}", file=summary)
    method = getattr(molset, "os_method", None)
    print(
        f"Metal oxidation states determined by: {' -> '.join(method) if method else 'N/A'}",
        file=summary,
    )
    for index, mol in enumerate(molset.moleclist):
        write_molecule_details(mol, file=summary, index=index)
    write_unique_species(molset, file=summary)
    return summary.getvalue().splitlines()


@blueprint.route("/process_structure/", methods=["POST"])
@calculation_logged
def process_structure_init():
    """Example view to process a crystal structure."""

    output = Capturing()
    if flask.request.files["structurefile"].filename == '':
        flask.flash("Please upload a file")
        return flask.redirect(flask.url_for("input_data"))

    
    # Get structure, file format, file content, and form data
    # (needed for additional information, e.g. cell in the case
    # of a XYZ file)
    fileformat = flask.request.form.get("fileformat", "unknown")
    form_data = dict(flask.request.form)
    structurefile = flask.request.files["structurefile"]
    file_ext = os.path.splitext(structurefile.filename)[1]
    system_type = flask.request.form.get("systemtype", "unknown")
    cif_connectivity = flask.request.form.get("cifConnectivity", "no")

    #resp = flask.make_response(flask.render_template(
    #    "user_templates/test.html", prueba = cif_connectivity
    #))
    #return resp

    #Only accept files whose extension is equal to the option extension selected
    #if fileformat in file_ext:

    try:
        token = Token(structurefile, fileformat)
    except InvalidXYZInput as exc:
        flask.flash("XYZ input error: " + str(exc))
        return flask.redirect(flask.url_for("input_data"))
    start_calculation_log(token, system_type, fileformat, form_data.get("xyzTotalCharge"))
    output += [
        "structure_file: "+ repr(dir(structurefile)),
        "form_data: "+ repr(form_data),
        "file_format: "+ repr(fileformat),
        "code: " + token.refcode,
        'req: '+ repr((flask.request, dir(flask.request))),
    ]



    #resp = flask.make_response(flask.render_template(
    #    "user_templates/test.html", struct_name=file_ext
    #))
    #return resp

    # The fileformat and the extension are consistent 
    if fileformat == "cif-pymatgen" and file_ext == ".cif":
            
        #tkn_path = flask.request.cookies.get('token_path', None)
        #if tkn_path is None:
        #    raise ValueError("not token?")
        #
        #token = Token.from_path(tkn_path)
        #if token is None:
        #    raise ValueError("session expired")
        #token.keepalive()

        def _safe_run(func, error_msg):
            """Utility to wrap save operations in try-except."""
            try:
                func()
            except Exception:
                logging.getLogger("compute.calculation").exception(error_msg)

        ##################
        ## Run cell2mol ##
        ##################

        #system_type ==  "unitcell"
        #cell = process_refcell(token.input_path, token.refcode, token.get_path())

        #system_type ==  "unitcell"
        if system_type == "unitcell":

            #flask.flash(token.get_path())
            #return flask.redirect(flask.url_for("input_data"))

            try:

                if cif_connectivity == "no":
                    config.USE_BOND_INFO = False

                #cell = process_unitcell(token.input_path, token.refcode, token.get_path(), cif_bond_info=False)
                cell = interpret_unitcell(token.input_path, token.refcode, token.get_path())
                if cell is None:
                    raise Cell2molNoResultError(token.refcode, "unit-cell")

                #cell.save(token.get_path()+token.refcode+'.json', format="json")
                _safe_run(lambda:cell.save(token.get_path()+'/'+token.refcode+'.json', format="json"), "Failed to save JSON")

                refcell_error = get_refcell_error(cell.reference)
                unitcell_error = get_unitcell_error(cell.unitcell)

                refcell_interpretation = _cell2mol_reference_interpretation(token, cell.reference, cell.unitcell)


                #resp = flask.make_response(flask.render_template(
                #    "user_templates/test.html", prueba=refcell_error
                #))
                #return resp

                save_cell(cell, 'gmol', token.get_path(), token.refcode)
                savemolecules_tools(cell.unitcell.moleclist, token.get_path(), 'xyz')
                savemolecules_tools(cell.unitcell.moleclist, token.get_path(), 'gmol')
                celldata = printing_text(cell, Capturing()) #empty

                cmp_lut = cell_cmp_lut(cell)
                ht_descs = cell_get_metal_desc(cell, cmp_lut)
                svgs = cell_to_svgs(cell, cmp_lut)
                compound_data = []
                for name,desc,svg in zip(cmp_lut.keys(), ht_descs, svgs):
                    # note: this line above uses the assumption that the order of items in a dict is predictable. Only true in recent-ish versions of python3
                    if desc != "":
                        compound_data.append((name, True, desc))
                    else:
                        compound_data.append((name, False, svg))

                ucellparams, xyzdata = cell_to_string_xyz(cell, cmp_lut)

                labels = []
                for mol in cell.unitcell.moleclist:
                    for atm in mol.atoms:
                        labels.append(atm.label)

                jmol_list_pos = molecules_list(cell)
                jmolCon = bond_order_connectivity(cell)
                jmol_list_species =species_list(cell) 

                #resp = flask.make_response(flask.render_template(
                #    "user_templates/test.html", prueba=cmp_lut
                #))
                #return resp

                token.keepalive()
                tkn_path = token.get_path()

                resp = flask.make_response(flask.render_template(
                    "user_templates/c2m-view.html",
                    celldata=celldata,
                    ucellparams=ucellparams,
                    compound_data=compound_data,
                    xyzdata=xyzdata,
                    labels=labels,
                    jmol_list_pos=jmol_list_pos,
                    jmol_list_species = jmol_list_species,
                    jmolCon = jmolCon,
                    totmol = len(cell.unitcell.moleclist),
                    enumerate=enumerate, len=len, zip=zip, # needed
                    struct_name=token.refcode,
                    unitcell_error_reconstruction = str(cell.unitcell.error_reconstruction),
                    unitcell_error_assign_charge = str(cell.unitcell.error_assign_charge),
                    unitcell_error_create_bonds = str(cell.unitcell.error_create_bonds),
                    unitcell_error_get_fragments = str(cell.unitcell.error_get_fragments),
                    unitcell_error_get_spin = str(cell.unitcell.error_get_spin),
                    refcell_error = refcell_error,
                    unitcell_error = unitcell_error,
                    refcell_interpretation = refcell_interpretation,
                    system_type=system_type,
                    calculation_log_preview=get_calculation_log_preview(token),
                ))
                resp.set_cookie("token_path",tkn_path,  secure=False,httponly=True,samesite='Strict') 
                return resp

            except Exception as e:
                mark_calculation_failure(e)
                msg = "Failure…"
                output += traceback.format_tb(e.__traceback__)
                output.append(repr(e))
                resp = flask.make_response(flask.render_template(
                        "user_templates/c2m-debug.html", msg=msg, output_lines=output,
                        calculation_log_preview=get_calculation_log_preview(token),
                        failure=describe_processing_error(e),
                        ))
                resp.set_cookie("token_path", token.get_path(), secure=False, httponly=True, samesite='Strict')
                return resp


        elif system_type == "reference":

            try:

                if cif_connectivity == "no":
                    config.USE_BOND_INFO = False

                #refMol = process_refcell(token.input_path, token.refcode, token.get_path(), cif_bond_info=False)
                #refMol = "test"
                refMol = interpret_reference(token.input_path, token.refcode, token.get_path())
                if refMol is None:
                    raise Cell2molNoResultError(token.refcode, "reference")
                savemolecules_reference_xyz(refMol, token.get_path())

                #Change cell to refMolec to avoid confussions
                _safe_run(lambda:refMol.save(token.get_path()+'/'+token.refcode+'.json', format="json"), "Failed to save JSON")

                refcell_error = get_refcell_error(refMol)
                #unitcell_error = get_unitcell_error(cell.unitcell)

                refcell_interpretation = _cell2mol_reference_interpretation(token, refMol)


                cmp_lut = cell_cmp_lut_base(refMol)
                ht_descs = cell_get_metal_desc_base(refMol, cmp_lut)
                svgs = cell_to_svgs_base(refMol, cmp_lut)
                plausible_charges = reference_plausible_charges(refMol, cmp_lut)
                compound_data = []
                for name,desc,svg,charges in zip(cmp_lut.keys(), ht_descs, svgs, plausible_charges):
                    # note: this line above uses the assumption that the order of items in a dict is predictable. Only true in recent-ish versions of python3
                    if desc != "":
                        compound_data.append((name, True, desc, cmp_lut[name][0][3], charges)) #cmp_lut[name][0][3] = chemical formula
                    else:
                        compound_data.append((name, False, svg, cmp_lut[name][0][3], charges)) #cmp_lut[name][0][3] = chemical formula 

                        #resp = flask.make_response(flask.render_template(
                        #    "user_templates/test.html", prueba=cmp_lut[name][0][3]
                        #))
                        #return resp




                ucellparams, xyzdata = refcell_to_string_xyz(refMol, cmp_lut)

                labels = []
                for mol in refMol.refmoleclist:
                    for atm in mol.atoms:
                        labels.append(atm.label)

                jmol_list_pos = molecules_list_reference(refMol)
                jmolCon = bond_order_connectivity_reference(refMol)
                jmol_list_species =species_list_reference(refMol) 

                token.keepalive()
                tkn_path = token.get_path()

                #resp = flask.make_response(flask.render_template(
                #    "user_templates/test.html", prueba="OK"
                #))
                #return resp

                resp = flask.make_response(flask.render_template(
                    "user_templates/c2m-view-refcell.html",
                    #celldata=celldata,
                    ucellparams=ucellparams,
                    compound_data=compound_data,
                    xyzdata=xyzdata,
                    labels=labels,
                    jmol_list_pos=jmol_list_pos,
                    jmol_list_species = jmol_list_species,
                    jmolCon = jmolCon,
                    totmol = len(refMol.refmoleclist),
                    enumerate=enumerate, len=len, zip=zip, # needed
                    struct_name=token.refcode,
                    unitcell_error_reconstruction = "",#str(refMol.error_reconstruction),
                    unitcell_error_assign_charge = str(refMol.error_assign_charge),
                    unitcell_error_create_bonds = str(refMol.error_create_bonds),
                    unitcell_error_get_fragments = "",#str(refMol.error_get_fragments),
                    unitcell_error_get_spin = str(refMol.error_get_spin),
                    refcell_error = refcell_error,
                    unitcell_error = "",#unitcell_error,
                    refcell_interpretation = refcell_interpretation,
                    system_type = system_type, 
                    calculation_log_preview=get_calculation_log_preview(token),
                ))
                resp.set_cookie("token_path",tkn_path,  secure=False,httponly=True,samesite='Strict') 
                return resp




            except Exception as e:
                mark_calculation_failure(e)
                msg = "Failure…"
                output += traceback.format_tb(e.__traceback__)
                output.append(repr(e))
                resp = flask.make_response(flask.render_template(
                        "user_templates/c2m-debug.html", msg=msg, output_lines=output,
                        calculation_log_preview=get_calculation_log_preview(token),
                        failure=describe_processing_error(e),
                        ))
                resp.set_cookie("token_path", token.get_path(), secure=False, httponly=True, samesite='Strict')
                return resp

            #save_cell(refMol, 'gmol', token.get_path(), token.refcode)
            #savemolecules_tools(refMol.refmoleclist, token.get_path(), 'xyz')
            #savemolecules_tools(refMol.refmoleclist, token.get_path(), 'gmol')
            #celldata = printing_text_refMol(refMol, Capturing()) #empty


            #jmol_list_pos = molecules_list_reference(refMol)

            #ucellparams, xyzdata = refcell_to_string_xyz(refMol)

            ##ref_error_file = glob.glob(os.path.join(token.get_path(), "reference_error_*.out"))

            ##if ref_error_file:

            ##    with open(ref_error_file[0], "r") as f:
            ##        ref_error_content = f.read()

            ##return flask.render_template("user_templates/c2m-debug.html", msg="ok", output_lines=ref_error_file, )
            #token.keepalive()
            #tkn_path = token.get_path()
            #resp = flask.make_response(flask.render_template(
            #    "user_templates/c2m-view-refcell.html",
            #    celldata=celldata,
            #    ucellparams=ucellparams,
            #    xyzdata=xyzdata,
            #    jmol_list_pos=jmol_list_pos,
            #    struct_name=token.refcode,
            #    #ref_error_content=ref_error_content,
            #))
            #resp.set_cookie("token_path",tkn_path,  secure=False,httponly=True,samesite='Strict') 
            #return resp

        else:
            flask.flash("The selected system type is not implemented yet.")
            return flask.redirect(flask.url_for("input_data"))



        #    
        #    #############################
        #    # STEP 3: Run cell2mol on infofile
        #    #############################
        #    #def run_cell2mol(run=False, out=None):
        #    cell = process_unitcell(input_path, name,token.refcode, token.get_path())
        #    #cell = process_unitcell(token.input_path, token.info_path, token.error_path)
        #    #cell = cell2mol(token.info_path, token.refcode, token.get_path(), 3)
        #    #with output as outt:
        #    #    print(cell, dir(cell), type(cell))
        #    #raise RuntimeError("unnamedrte")
        #    save_cell(cell, 'gmol', token.get_path())
        #    savemolecules(cell.moleclist, token.get_path(), 'xyz')
        #    savemolecules(cell.moleclist, token.get_path(), 'gmol')

        #    resp = flask.make_response(flask.render_template(
        #        "user_templates/test.html",
        #    ))
        #    return resp


            #if os.path.exists(token.info_path):
            #    labels, pos, lfracs, fracs, cellvec, cellparam = readinfo(token.info_path)
            #    with open(token.info_path, 'r') as f:
            #        infodata = f.read()
            #        output.append(infodata)
            #    with open(token.analysis_path, 'r') as f:
            #        celldata = f.read()
            #        output.append(celldata)
            #    with open(token.cell_path, 'rb') as f:
            #        cell = pickle.load(f)


            #    cmp_lut = cell_cmp_lut(cell)
            #    ht_descs = cell_get_metal_desc(cell, cmp_lut)
            #    svgs = cell_to_svgs(cell, cmp_lut)
            #    compound_data = []
            #    for name,desc,svg in zip(cmp_lut.keys(), ht_descs, svgs):
            #        # note: this line above uses the assumption that the order of items in a dict is predictable. Only true in recent-ish versions of python3
            #        if desc != "":
            #            compound_data.append((name, True, desc))
            #        else:
            #            compound_data.append((name, False, svg))
            #    

            #    #ucellparams, xyzdata = cell_to_string_xyz(cell, cmp_lut)
            #    ucellparams, xyzdata = cell_to_string_xyz(cell, cmp_lut)

            #    #string used by jsmol to define the molecules/complexes, and the connectivity respectively
            #    jmol_list_pos = molecules_list(cell)
            #    jmolCon = bond_order_connectivity(cell)
            #    jmol_list_species =species_list(cell) 

            #    
            #resp = flask.make_response(flask.render_template(
            #    "user_templates/c2m-view.html",
            #    output_lines=output,
            #    infodata=infodata.strip(),
            #    celldata=celldata,
            #    ucellparams=ucellparams,
            #    compound_data=compound_data,
            #    xyzdata=xyzdata,
            #    labels=labels,
            #    pos=pos,
            #    cellvec=cellvec,
            #    cellparam=cellparam,
            #    jmol_list_pos=jmol_list_pos,
            #    jmol_list_species = jmol_list_species,
            #    jmolCon = jmolCon,
            #    totmol = len(cell.moleclist),
            #    enumerate=enumerate, len=len, zip=zip, # needed
            #    #token_path=tkn_path.replace('/','_'), #blueprint.url_for('process_structure','analysis', token=tkn_path.replace('/','_')),
            #    struct_name=token.refcode,
            #))
            #return resp



            #############################
            # STEP 2: Run cell2info
            #############################
            #with output as _out:
            #    cif_2_info(token.input_path, token.info_path, token.error_path)
            #error = False
            #with open(token.error_path, 'r') as err:
            #    for line in err.readlines():
            #        if "Error" in line:
            #            output.append(line)
            #            error = True
            #if error :
            #    output.append(f"Parsing of .cif file {token.input_path} failed due to the error above.")
            #else :
            #    output.append(f"Infofile {token.info_path} generated from {token.input_path} succesfully.")

        #elif fileformat == "info":
        #    #############################
        #    # STEP 2: Rename input file as info file
        #    #############################
        #    os.rename(token.input_path, token.info_path)

        #elif fileformat =="unknown":
        #    flask.flash("Please select a file fomrat and upload a file ")
        #    return flask.redirect(flask.url_for("input_data"))

        #else :
        #    flask.flash("Please upload a file in the selected format with the proper extension (.'{}') ".format(fileformat))
        #    return flask.redirect(flask.url_for("input_data"))


        #with open(token.info_path, 'r') as f:
        #    infodata = f.read()
        #    output.append(infodata)
        #
        #tkn_path = token.get_path()

        #resp = flask.make_response(flask.render_template(
        #    "user_templates/c2m-infopage.html",
        #    output_lines=output,
        #    infodata=infodata,
        #    enumerate=enumerate, len=len,  # why is this needed?
        #    #token_path=tkn_path.replace('/','_'), #blueprint.url_for('process_structure','analysis', token=tkn_path.replace('/','_')),
        #    struct_name=token.refcode,
        #))
        #resp.set_cookie("token_path",tkn_path,  secure=False,httponly=True,samesite='Strict') #secure must be True for final version
        #return resp
    
    elif fileformat == "xyz-ase" and file_ext == ".xyz":

        system_type = "molecule"
        try:
            total_charge = int(flask.request.form.get("xyzTotalCharge", "0"))
        except (ValueError, TypeError):
            flask.flash("Please enter an integer total charge for the XYZ structure")
            return flask.redirect(flask.url_for("input_data"))

        if system_type == "molecule":
            try:
                molset = interpret_molecule_xyz(token.input_path, token.refcode, token.get_path(), total_charge)
                if molset is None:
                    raise Cell2molNoResultError(token.refcode, "molecule")
                moleclist = get_moleculeset_molecules(molset)

                _safe_save_molecule_json(molset, token.cell_path)
                savemolecules_moleculeset_xyz(moleclist, token.get_path(), token.refcode)

                molecule_interpretation = _cell2mol_molecule_interpretation(token, molset)

                if moleculeset_is_resolved(moleclist):
                    molecule_cell = make_moleculeset_cell(moleclist)
                    cmp_lut = cell_cmp_lut(molecule_cell)
                    ht_descs = cell_get_metal_desc(molecule_cell, cmp_lut)
                    svgs = cell_to_svgs(molecule_cell, cmp_lut)
                    compound_data = []
                    for name, desc, svg in zip(cmp_lut.keys(), ht_descs, svgs):
                        tpl = cmp_lut[name][0]
                        mol = molecule_cell.unitcell.moleclist[tpl[0]]
                        if tpl[1] == 'm':
                            charge = mol.metals[tpl[2]].charge
                            charge_label = 'Oxidation state'
                        elif tpl[1] == 'l':
                            charge = mol.ligands[tpl[2]].totcharge
                            charge_label = 'Charge'
                        else:
                            charge = mol.totcharge
                            charge_label = 'Charge'
                        charge_value = ('+' if charge > 0 else '') + str(int(charge)) if charge is not None else 'Not available'
                        if desc != "":
                            compound_data.append((name, True, desc, charge_label, charge_value))
                        else:
                            compound_data.append((name, False, svg, charge_label, charge_value))
                    jmolCon = bond_order_connectivity(molecule_cell)
                    jmol_list_species = species_list(molecule_cell)
                else:
                    compound_data, jmol_list_species = moleculeset_partial_species(moleclist)
                    jmolCon = bond_order_connectivity_moleculeset(moleclist)

                xyzdata = moleculeset_to_string_xyz(moleclist, token.refcode)
                jmol_list_pos = molecules_list_moleculeset(moleclist, token.refcode)
                molecule_errors = get_moleculeset_errors(molset)
                molecule_errors.insert(0, molecule_status_message(molset, moleculeset_is_resolved(moleclist)))

                token.keepalive()
                tkn_path = token.get_path()
                resp = flask.make_response(flask.render_template(
                    "user_templates/c2m-view-xyzmol.html",
                    molecule_interpretation=molecule_interpretation,
                    compound_data=compound_data,
                    xyzdata=xyzdata,
                    jmol_list_pos=jmol_list_pos,
                    jmol_list_species=jmol_list_species,
                    jmolCon=jmolCon,
                    totmol=len(moleclist),
                    enumerate=enumerate, len=len, zip=zip,
                    struct_name=token.refcode,
                    system_type=system_type,
                    total_charge=total_charge,
                    molecule_errors=molecule_errors,
                    calculation_log_preview=get_calculation_log_preview(token),
                ))
                resp.set_cookie("token_path", tkn_path, secure=False, httponly=True, samesite='Strict')
                return resp

            except Exception as e:
                mark_calculation_failure(e)
                msg = "Failure…"
                output += traceback.format_tb(e.__traceback__)
                output.append(repr(e))
                resp = flask.make_response(flask.render_template(
                    "user_templates/c2m-debug.html", msg=msg, output_lines=output,
                    calculation_log_preview=get_calculation_log_preview(token),
                    failure=describe_processing_error(e),
                ))
                resp.set_cookie("token_path", token.get_path(), secure=False, httponly=True, samesite='Strict')
                return resp

    elif fileformat =="unknown":
        flask.flash("Please select a file format.")
        return flask.redirect(flask.url_for("input_data"))

    else:

        flask.flash("Please upload a file in format '{}' with the proper extension (.'{}') ".format(fileformat, fileformat))
        return flask.redirect(flask.url_for("input_data"))



#@blueprint.route("/analysis", methods=["GET"])
#def process_structure_analysis():
#
#    output = Capturing()
#    try:
#        #tkn_path = flask.request.args['token'].replace('_','/')
#        tkn_path = flask.request.cookies.get('token_path', None)
#        if tkn_path is None:
#            raise ValueError("not token?")
#        
#        token = Token.from_path(tkn_path)
#        if token is None:
#            raise ValueError("session expired")
#        token.keepalive()
#        
#        #############################
#        # STEP 3: Run cell2mol on infofile
#        #############################
#        #def run_cell2mol(run=False, out=None):
#
#        if os.path.exists(token.info_path):
#            with open(token.info_path, 'r') as f:
#                infodata = f.read()
#                output.append(infodata)
#
#            with output as _out:
#                cell = cell2mol(token.info_path, token.refcode, token.get_path(), 3)
#                #with output as outt:
#                #    print(cell, dir(cell), type(cell))
#                #raise RuntimeError("unnamedrte")
#                save_cell(cell, 'gmol', token.get_path())
#                savemolecules(cell.moleclist, token.get_path(), 'xyz')
#                savemolecules(cell.moleclist, token.get_path(), 'gmol')
#            # This line removes cell2mol output from printout, it should be extend not redefine
#            #output += [f"For input {token.input_path}"]
#            celldata = printing_text(cell, Capturing())
#            with open(token.analysis_path, 'w') as f:
#                for line in celldata:
#                    f.write(line+"\n")
#        else:
#            raise ValueError("plz")
#            #output.extend([f"Please, wait until cell2info has finished for this input. Could not find {token.info_path}."])
#
#        #token.remove()
#        resp = flask.make_response(flask.render_template(
#            "user_templates/c2m-analysis.html",
#            output_lines=output,
#            infodata=infodata.strip(),
#            celldata='\n'.join(celldata).strip(),
#            enumerate=enumerate, len=len,  # why TF is this needed?????
#            #token_path=tkn_path.replace('/','_'), #blueprint.url_for('process_structure','analysis', token=tkn_path.replace('/','_')),
#            struct_name=token.refcode,
#        ))
#        return resp
#        #############################
#        # STEP 4: Display structures
#        #############################    
#        #def display_mol(run=False, out=None):
#
#        #with output as _out:
#        #    path = find_gmol(run=run)
#        #    cell = pickle.load(open(path, "rb"))
#        #with out:
#        #    printing_structure_cell(cell)
#    
#    except Exception as err:
#        msg = "Failure…"
#        output += traceback.format_tb(err.__traceback__)
#        output.append(repr(err))
#        return flask.render_template(
#            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
#        )
#    except:
#        msg = "Failure…"
#        output.append("unknown error")
#        return flask.render_template(
#            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
#        )
#    




#@blueprint.route("/view-gmol", methods=["GET"])
#def process_structure_view():
#
#    output = Capturing()
#    try:
#        #tkn_path = flask.request.args['token'].replace('_','/')
#        tkn_path = flask.request.cookies.get('token_path', None)
#        if tkn_path is None:
#            raise ValueError("not token?")
#        
#        token = Token.from_path(tkn_path)
#        if token is None:
#            raise ValueError("session expired")
#        token.keepalive()
#        
#        #############################
#        # STEP 3: Run cell2mol on infofile
#        #############################
#        #def run_cell2mol(run=False, out=None):
#
#        if os.path.exists(token.info_path):
#            labels, pos, lfracs, fracs, cellvec, cellparam = readinfo(token.info_path)
#            with open(token.info_path, 'r') as f:
#                infodata = f.read()
#                output.append(infodata)
#            with open(token.analysis_path, 'r') as f:
#                celldata = f.read()
#                output.append(celldata)
#            with open(token.cell_path, 'rb') as f:
#                cell = pickle.load(f)
#
#
#            cmp_lut = cell_cmp_lut(cell)
#            ht_descs = cell_get_metal_desc(cell, cmp_lut)
#            svgs = cell_to_svgs(cell, cmp_lut)
#            compound_data = []
#            for name,desc,svg in zip(cmp_lut.keys(), ht_descs, svgs):
#                # note: this line above uses the assumption that the order of items in a dict is predictable. Only true in recent-ish versions of python3
#                if desc != "":
#                    compound_data.append((name, True, desc))
#                else:
#                    compound_data.append((name, False, svg))
#            
#
#            #ucellparams, xyzdata = cell_to_string_xyz(cell, cmp_lut)
#            ucellparams, xyzdata = cell_to_string_xyz(cell, cmp_lut)
#
#            #string used by jsmol to define the molecules/complexes, and the connectivity respectively
#            jmol_list_pos = molecules_list(cell)
#            jmolCon = bond_order_connectivity(cell)
#            jmol_list_species =species_list(cell) 
#
#        else:
#            raise ValueError("Please wait until cell2info finishes.")
#            #output.extend([f"Please, wait until cell2info has finished for this input. Could not find {token.info_path}."])
#
#        #token.remove()
#        resp = flask.make_response(flask.render_template(
#            "user_templates/c2m-view_OLD.html",
#            output_lines=output,
#            infodata=infodata.strip(),
#            celldata=celldata,
#            ucellparams=ucellparams,
#            compound_data=compound_data,
#            xyzdata=xyzdata,
#            labels=labels,
#            pos=pos,
#            cellvec=cellvec,
#            cellparam=cellparam,
#            jmol_list_pos=jmol_list_pos,
#            jmol_list_species = jmol_list_species,
#            jmolCon = jmolCon,
#            totmol = len(cell.moleclist),
#            enumerate=enumerate, len=len, zip=zip, # needed
#            #token_path=tkn_path.replace('/','_'), #blueprint.url_for('process_structure','analysis', token=tkn_path.replace('/','_')),
#            struct_name=token.refcode,
#        ))
#        return resp
#        #############################
#        # STEP 4: Display structures
#        #############################    
#        #def display_mol(run=False, out=None):
#
#        #with output as _out:
#        #    path = find_gmol(run=run)
#        #    cell = pickle.load(open(path, "rb"))
#        #with out:
#        #    printing_structure_cell(cell)
#    
#    except Exception as err:
#        msg = "Failure…"
#        output += traceback.format_tb(err.__traceback__)
#        output.append(repr(err))
#        return flask.render_template(
#            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
#        )
#    except:
#        msg = "Failure…"
#        output.append("unknown error")
#        return flask.render_template(
#            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
#        )


#>>> D O W N L O A D   C E L L 2 M O L   O U T P U T <<<
@blueprint.route("/process_structure/download-log", methods=["GET"])
def process_structure_download_log():
    tkn_path = flask.request.cookies.get('token_path', None)
    token = Token.from_path(tkn_path) if tkn_path else None
    if token is None:
        return flask.make_response("Calculation session expired", 404)
    token.keepalive()
    log_path = calculation_log_path(token)
    if not os.path.isfile(log_path):
        return flask.make_response("Calculation log not available", 404)
    return send_file(log_path, as_attachment=True,
                     download_name="calculation.log", mimetype="text/plain")


@blueprint.route("/process_structure/download-c2moutput", methods=["GET"])
def process_structure_download_c2m_output():

    output = Capturing()
    try:
        tkn_path = flask.request.cookies.get('token_path', None)
        if tkn_path is None:
            raise ValueError("no token?")
        token = Token.from_path(tkn_path)
        if token is None:
            raise ValueError("session expired")
        output.append(token.cell_path)
        token.keepalive()
                
        #headers = {"Content-Disposition": f"attachment; filename=cell2mol.out"}
        #with open(token.cell_path, 'rb') as f:
        #    body = f.read()
        #return flask.make_response((body, headers))
        return send_file(
            token.get_path()+'/cell2mol.out',
            as_attachment=True,
            download_name="cell2mol.out",
            mimetype="text/plain"
        )
    
        #if True or res.status_code >= 400:
        #    output.append(repr(res))
        #    raise ValueError("Bad status code: {:d}".format(res.status_code))
        #else:
        #    return res
        
     
    except Exception as err:
        msg = "Failure…"
        output.append(repr(err))
        return flask.render_template(
            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
        )
    except:
        msg = "Failure…"
        output.append("unknown error")
        return flask.render_template(
            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
        )

@blueprint.route("/process_structure/download-referenceSummary", methods=["GET"])
def process_structure_download_reference_summary():

    output = Capturing()
    try:
        tkn_path = flask.request.cookies.get('token_path', None)
        if tkn_path is None:
            raise ValueError("no token?")
        token = Token.from_path(tkn_path)
        if token is None:
            raise ValueError("session expired")
        output.append(token.cell_path)
        token.keepalive()
                
        #headers = {"Content-Disposition": f"attachment; filename=cell2mol.out"}
        #with open(token.cell_path, 'rb') as f:
        #    body = f.read()
        #return flask.make_response((body, headers))
        return send_file(
            token.get_path()+'/reference_summary.out',
            as_attachment=True,
            download_name="reference_summary.out",
            mimetype="text/plain"
        )
    
        #if True or res.status_code >= 400:
        #    output.append(repr(res))
        #    raise ValueError("Bad status code: {:d}".format(res.status_code))
        #else:
        #    return res
        
     
    except Exception as err:
        msg = "Failure…"
        output.append(repr(err))
        return flask.render_template(
            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
        )
    except:
        msg = "Failure…"
        output.append("unknown error")
        return flask.render_template(
            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
        )



#>>> D O W N L O A D   J S O N <<<
@blueprint.route("/process_example_structure/download-json", methods=["GET"])
@blueprint.route("/process_structure/download-json", methods=["GET"])
def process_structure_download_json():

    output = Capturing()
    try:
        tkn_path = flask.request.cookies.get('token_path', None)
        if tkn_path is None:
            raise ValueError("no token?")
        token = Token.from_path(tkn_path)
        if token is None:
            raise ValueError("session expired")
        output.append(token.cell_path)
        token.keepalive()
                
        headers = {"Content-Disposition": f"attachment; filename=Cell_{token.refcode:s}.json"}
        with open(token.cell_path, 'rb') as f:
            body = f.read()
        return flask.make_response((body, headers))
    
        if True or res.status_code >= 400:
            output.append(repr(res))
            raise ValueError("Bad status code: {:d}".format(res.status_code))
        else:
            return res
        
     
    except Exception as err:
        msg = "Failure…"
        output.append(repr(err))
        return flask.render_template(
            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
        )
    except:
        msg = "Failure…"
        output.append("unknown error")
        return flask.render_template(
            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
        )

#>>> D O W N L O A D   C E L L <<<

@blueprint.route("/process_structure/download-gmol", methods=["GET"])
def process_structure_download_gmol():

    output = Capturing()
    try:
        tkn_path = flask.request.cookies.get('token_path', None)
        if tkn_path is None:
            raise ValueError("no token?")
        token = Token.from_path(tkn_path)
        if token is None:
            raise ValueError("session expired")
        output.append(token.cell_path)
        token.keepalive()
                
        headers = {"Content-Disposition": f"attachment; filename=Cell_{token.refcode:s}.gmol"}
        with open(token.cell_path, 'rb') as f:
            body = f.read()
        return flask.make_response((body, headers))
    
        if True or res.status_code >= 400:
            output.append(repr(res))
            raise ValueError("Bad status code: {:d}".format(res.status_code))
        else:
            return res
        
     
    except Exception as err:
        msg = "Failure…"
        output.append(repr(err))
        return flask.render_template(
            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
        )
    except:
        msg = "Failure…"
        output.append("unknown error")
        return flask.render_template(
            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
        )

# >>> D O W N L O A D   M O L <<<
@blueprint.route("/process_example_structure/download-xyzselected", methods=["GET", "POST"])
@blueprint.route("/process_structure/download-xyzselected", methods=["GET", "POST"])
def process_structure_download_xyzselected():
    output = Capturing()
    try:
        tkn_path = flask.request.cookies.get('token_path', None)
        if tkn_path is None:
            raise ValueError("no token?")
        token = Token.from_path(tkn_path)
        if token is None:
            raise ValueError("session expired")
        output.append(token.cell_path)
        token.keepalive()
        #Get name from last part of atmPos
        xyzfile = str(flask.request.form.get("atmPos","unknown")).split('$')[1]
        #add the extension
        xyzfile = xyzfile + ".xyz"
        headers={'Content-disposition': 'attachment; filename='+ xyzfile}
        #read from the tmp directory
        with open(tkn_path+"/"+xyzfile, 'rb') as f:
            body = f.read()
        #return the file
        return flask.make_response(body, headers)

    except Exception as err:
        msg = "Failure…"
        output.append(repr(err))
        return flask.render_template(
            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
        )
    except:
        msg = "Failure…"
        output.append("unknown error")
        return flask.render_template(
            "user_templates/c2m-debug.html", msg=msg, output_lines=output,
        )
    



@blueprint.route("/keepalive", methods=["GET"])
def session_keepalive():

    tkn_path = flask.request.cookies.get('token_path', None)
    if tkn_path is None:
        return flask.make_response("<html>No token path cookie</html>", 404)

    token = Token.from_path(tkn_path)
    if token is None:
        return flask.make_response("<html>session expired</html>", 404)
    else:
        token.keepalive()
        return flask.make_response("", 200)

# >>> YOXKUS example <<<

@blueprint.route("/process_example_structure/", methods=["POST"])
@calculation_logged
def process_structure_example_init():

    example_selected = flask.request.form.get('selected_option')
    if example_selected == "YOXKUS_unitcell":

        output = Capturing()
        example_path = os.path.join(os.path.dirname(__file__), "examples", "cif", "YOXKUS.cif")
        if not os.path.isfile(example_path):
            flask.flash("The preloaded YOXKUS example file is missing. Please contact the administrator.")
            return flask.redirect(flask.url_for("input_data"))

        with open(example_path, 'rb') as structurefile:
            token = Token(FileStorage(structurefile, filename="YOXKUS.cif"), "cif-pymatgen")

        start_calculation_log(token, "unitcell", "cif-pymatgen")
        previous_bond_info = config.USE_BOND_INFO

        try:
            config.USE_BOND_INFO = False
            cell = interpret_unitcell(token.input_path, token.refcode, token.get_path())
            if cell is None:
                raise Cell2molNoResultError(token.refcode, "unit-cell")

            try:
                cell.save(token.get_path()+'/'+token.refcode+'.json', format="json")
            except Exception:
                logging.getLogger(__name__).exception("Failed to save JSON")

            refcell_error = get_refcell_error(cell.reference)
            unitcell_error = get_unitcell_error(cell.unitcell)
            refcell_interpretation = _cell2mol_reference_interpretation(token, cell.reference, cell.unitcell)

            save_cell(cell, 'gmol', token.get_path(), token.refcode)
            savemolecules_tools(cell.unitcell.moleclist, token.get_path(), 'xyz')
            savemolecules_tools(cell.unitcell.moleclist, token.get_path(), 'gmol')
            celldata = printing_text(cell, Capturing())

            cmp_lut = cell_cmp_lut(cell)
            ht_descs = cell_get_metal_desc(cell, cmp_lut)
            svgs = cell_to_svgs(cell, cmp_lut)
            compound_data = []
            for name,desc,svg in zip(cmp_lut.keys(), ht_descs, svgs):
                if desc != "":
                    compound_data.append((name, True, desc))
                else:
                    compound_data.append((name, False, svg))

            ucellparams, xyzdata = cell_to_string_xyz(cell, cmp_lut)

            labels = []
            for mol in cell.unitcell.moleclist:
                for atm in mol.atoms:
                    labels.append(atm.label)

            jmol_list_pos = molecules_list(cell)
            jmolCon = bond_order_connectivity(cell)
            jmol_list_species = species_list(cell)

            token.keepalive()
            resp = flask.make_response(flask.render_template(
                "user_templates/c2m-view.html",
                celldata=celldata,
                ucellparams=ucellparams,
                compound_data=compound_data,
                xyzdata=xyzdata,
                labels=labels,
                jmol_list_pos=jmol_list_pos,
                jmol_list_species=jmol_list_species,
                jmolCon=jmolCon,
                totmol=len(cell.unitcell.moleclist),
                enumerate=enumerate, len=len, zip=zip,
                struct_name=token.refcode,
                unitcell_error_reconstruction=str(cell.unitcell.error_reconstruction),
                unitcell_error_assign_charge=str(cell.unitcell.error_assign_charge),
                unitcell_error_create_bonds=str(cell.unitcell.error_create_bonds),
                unitcell_error_get_fragments=str(cell.unitcell.error_get_fragments),
                unitcell_error_get_spin=str(cell.unitcell.error_get_spin),
                refcell_error=refcell_error,
                unitcell_error=unitcell_error,
                refcell_interpretation=refcell_interpretation,
                system_type="unitcell",
                calculation_log_preview=get_calculation_log_preview(token),
            ))
            resp.set_cookie("token_path", token.get_path(), secure=False, httponly=True, samesite='Strict')
            return resp

        except Exception as e:
            mark_calculation_failure(e)
            msg = "Failure…"
            output += traceback.format_tb(e.__traceback__)
            output.append(repr(e))
            resp = flask.make_response(flask.render_template(
                "user_templates/c2m-debug.html", msg=msg, output_lines=output,
                calculation_log_preview=get_calculation_log_preview(token),
                failure=describe_processing_error(e),
            ))
            resp.set_cookie("token_path", token.get_path(), secure=False, httponly=True, samesite='Strict')
            return resp
        finally:
            config.USE_BOND_INFO = previous_bond_info

    elif example_selected in ("ammonium_chloride_water_molecule", "NH4Cl_H2O_molecule"):

        output = Capturing()
        example_path = os.path.join(os.path.dirname(__file__), "examples", "xyz", "NH4Cl_H2O.xyz")
        if not os.path.isfile(example_path):
            flask.flash("The preloaded NH4Cl_H2O example file is missing. Please contact the administrator.")
            return flask.redirect(flask.url_for("input_data"))

        with open(example_path, 'rb') as structurefile:
            token = Token(FileStorage(structurefile, filename="NH4Cl_H2O.xyz"), "xyz-ase")

        total_charge = 0
        start_calculation_log(token, "molecule", "xyz-ase", total_charge)

        try:
            molset = interpret_molecule_xyz(token.input_path, token.refcode, token.get_path(), total_charge)
            if molset is None:
                raise Cell2molNoResultError(token.refcode, "molecule")
            moleclist = get_moleculeset_molecules(molset)

            _safe_save_molecule_json(molset, token.cell_path)
            savemolecules_moleculeset_xyz(moleclist, token.get_path(), token.refcode)

            molecule_interpretation = _cell2mol_molecule_interpretation(token, molset)

            if moleculeset_is_resolved(moleclist):
                molecule_cell = make_moleculeset_cell(moleclist)
                cmp_lut = cell_cmp_lut(molecule_cell)
                ht_descs = cell_get_metal_desc(molecule_cell, cmp_lut)
                svgs = cell_to_svgs(molecule_cell, cmp_lut)
                compound_data = []
                for name, desc, svg in zip(cmp_lut.keys(), ht_descs, svgs):
                    tpl = cmp_lut[name][0]
                    mol = molecule_cell.unitcell.moleclist[tpl[0]]
                    if tpl[1] == 'm':
                        charge = mol.metals[tpl[2]].charge
                        charge_label = 'Oxidation state'
                    elif tpl[1] == 'l':
                        charge = mol.ligands[tpl[2]].totcharge
                        charge_label = 'Charge'
                    else:
                        charge = mol.totcharge
                        charge_label = 'Charge'
                    charge_value = ('+' if charge > 0 else '') + str(int(charge)) if charge is not None else 'Not available'
                    if desc != "":
                        compound_data.append((name, True, desc, charge_label, charge_value))
                    else:
                        compound_data.append((name, False, svg, charge_label, charge_value))
                jmolCon = bond_order_connectivity(molecule_cell)
                jmol_list_species = species_list(molecule_cell)
            else:
                compound_data, jmol_list_species = moleculeset_partial_species(moleclist)
                jmolCon = bond_order_connectivity_moleculeset(moleclist)

            xyzdata = moleculeset_to_string_xyz(moleclist, token.refcode)
            jmol_list_pos = molecules_list_moleculeset(moleclist, token.refcode)
            molecule_errors = get_moleculeset_errors(molset)
            molecule_errors.insert(0, molecule_status_message(molset, moleculeset_is_resolved(moleclist)))

            token.keepalive()
            tkn_path = token.get_path()
            resp = flask.make_response(flask.render_template(
                "user_templates/c2m-view-xyzmol.html",
                molecule_interpretation=molecule_interpretation,
                compound_data=compound_data,
                xyzdata=xyzdata,
                jmol_list_pos=jmol_list_pos,
                jmol_list_species=jmol_list_species,
                jmolCon=jmolCon,
                totmol=len(moleclist),
                enumerate=enumerate, len=len, zip=zip,
                struct_name=token.refcode,
                system_type="molecule",
                total_charge=total_charge,
                molecule_errors=molecule_errors,
                calculation_log_preview=get_calculation_log_preview(token),
            ))
            resp.set_cookie("token_path", tkn_path, secure=False, httponly=True, samesite='Strict')
            return resp

        except Exception as e:
            mark_calculation_failure(e)
            msg = "Failure…"
            output += traceback.format_tb(e.__traceback__)
            output.append(repr(e))
            resp = flask.make_response(flask.render_template(
                "user_templates/c2m-debug.html", msg=msg, output_lines=output,
                calculation_log_preview=get_calculation_log_preview(token),
                failure=describe_processing_error(e),
            ))
            resp.set_cookie("token_path", token.get_path(), secure=False, httponly=True, samesite='Strict')
            return resp

    else:

        flask.flash("Please select an example.")
        return flask.redirect(flask.url_for("input_data"))


@blueprint.route("/analysis_YOXKUS", methods=["GET"])
def process_structure_analysis_YOXKUS():

    output = Capturing()
    
    #############################
    # STEP 3: Run cell2mol on infofile
    #############################

    info_path="/home/app/code/webservice/compute/examples/cif/info.txt"
    get_path="/home/app/code/webservice/compute/examples/cif/"
    analysis_path="/home/app/code/webservice/compute/examples/cif/cell.txt"
    with open(info_path, 'r') as f:
        infodata = f.read()
        output.append(infodata)


    with output as _out:
        cell = cell2mol(info_path, 'YOXKUS', get_path, 3)
        #with output as outt:
        #    print(cell, dir(cell), type(cell))
        #raise RuntimeError("unnamedrte")
        save_cell(cell, 'gmol', get_path)
        savemolecules(cell.moleclist, get_path, 'xyz')
        savemolecules(cell.moleclist, get_path, 'gmol')
    # This line removes cell2mol output from printout, it should be extend not redefine
    #output += [f"For input {token.input_path}"]
    celldata = printing_text(cell, Capturing())
    with open(analysis_path, 'w') as f:
        for line in celldata:
            f.write(line+"\n")

    resp = flask.make_response(flask.render_template(
        "user_templates/c2m-analysis-YOXKUS.html",
        output_lines=output,
        infodata=infodata.strip(),
        celldata='\n'.join(celldata).strip(),
        enumerate=enumerate, len=len,  # why TF is this needed?????
        #token_path=tkn_path.replace('/','_'), #blueprint.url_for('process_structure','analysis', token=tkn_path.replace('/','_')),
        struct_name="YOXKUS"#token.refcode,
    ))
    return resp





@blueprint.route("compute/view-gmol-YOXKUS/", methods=["GET"])
def process_structure_view_YOXKUS():

    #output = Capturing()

    #input_path="/home/app/code/webservice/compute/examples/cif/YOXKUS.cif"

    #try:
    #    cell = process_unitcell(input_path, "YOXKUS", "/home/app/code/webservice/compute/examples/cif/results")
    #except Exception as e:
    #    exit_with_error_exception(e)

    #celldata = printing_text(cell, Capturing()) #empty

    #cmp_lut = cell_cmp_lut(cell)
    #ht_descs = cell_get_metal_desc(cell, cmp_lut)
    #svgs = cell_to_svgs(cell, cmp_lut)
    #compound_data = []
    #for name,desc,svg in zip(cmp_lut.keys(), ht_descs, svgs):
    #    if desc != "":
    #        compound_data.append((name, True, desc))
    #    else:
    #        compound_data.append((name, False, svg))


    #ucellparams, xyzdata = cell_to_string_xyz(cell, cmp_lut)

    #labels = []
    #for mol in cell.moleclist:
    #    for atm in mol.atoms:
    #        labels.append(atm.label)


    #jmol_list_pos = molecules_list(cell)
    #jmolCon = bond_order_connectivity(cell)
    #jmol_list_species =species_list(cell) 

    resp = flask.make_response(flask.render_template(
        "user_templates/test.html",
        prueba="HOLA",
    #    "user_templates/c2m-view-YOXKUS.html",
    #    celldata=celldata,
    #    ucellparams=ucellparams,
    #    compound_data=compound_data,
    #    xyzdata=xyzdata,
    #    labels=labels,
    #    jmol_list_pos=jmol_list_pos,
    #    jmol_list_species = jmol_list_species,
    #    jmolCon = jmolCon,
    #    totmol = len(cell.moleclist),
    #    enumerate=enumerate, len=len, zip=zip, # needed
    #    struct_name=token.refcode,
    ))
    return resp

    
@blueprint.route("/process_example_structure/download-gmol-YOXKUS", methods=["GET"])
def process_structure_download_gmol_YOXKUS():

    output = Capturing()

    cell_path="/home/app/code/webservice/compute/examples/cif/results/Cell_YOXKUS.gmol"
    headers = {"Content-Disposition": f"attachment; filename=Cell_YOXKUS.gmol"}
    with open(cell_path, 'rb') as f:
        body = f.read()
    return flask.make_response((body, headers))

@blueprint.route("/process_example_structure/download-xyzselected-YOXKUS", methods=["GET", "POST"])
def process_structure_download_xyzselected_YOXKUS():
    output = Capturing()

    #Get name from last part of atmPos
    tkn_path="/home/app/code/webservice/compute/examples/cif/results/"
    xyzfile = str(flask.request.form.get("atmPos","unknown")).split('$')[1]
    #add the extension
    xyzfile = xyzfile + ".xyz"
    headers={'Content-disposition': 'attachment; filename='+ xyzfile}
    #read from the tmp directory
    with open(tkn_path+"/"+xyzfile, 'rb') as f:
        body = f.read()
    #return the file
    return flask.make_response(body, headers)

