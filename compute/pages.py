import traceback
import pickle
from tools_barebone.structure_importers import get_structure_tuple, UnknownFormatError
import flask
from .interface import *
from .tokens import monitoring, Token
from cell2mol.readwrite import savemolecules, writexyz

import os

import webbrowser

blueprint = flask.Blueprint("compute", __name__, url_prefix="/compute")
    
@blueprint.route("/process_structure/", methods=["POST"])
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

    #Only accept files whose extension is equal to the option extension selected
    if fileformat in file_ext:

        token = Token(structurefile)        
        output += [
            "structure_file: "+ repr(dir(structurefile)),
            "form_data: "+ repr(form_data),
            "file_format: "+ repr(fileformat),
            "code: " + token.refcode,
            'req: '+ repr((flask.request, dir(flask.request))),
        ]


        if fileformat == "cif":
            #############################
            # STEP 2: Run cell2info
            #############################
            with output as _out:
                cif_2_info(token.input_path, token.info_path, token.error_path)
            error = False
            with open(token.error_path, 'r') as err:
                for line in err.readlines():
                    if "Error" in line:
                        output.append(line)
                        error = True
            if error :
                output.append(f"Parsing of .cif file {token.input_path} failed due to the error above.")
            else :
                output.append(f"Infofile {token.info_path} generated from {token.input_path} succesfully.")

        elif fileformat == "info":
            #############################
            # STEP 2: Rename input file as info file
            #############################
            os.rename(token.input_path, token.info_path)

        elif fileformat =="unknown":
            flask.flash("Please select a file fomrat and upload a file ")
            return flask.redirect(flask.url_for("input_data"))

        else :
            flask.flash("Please upload a file in the selected format with the proper extension (.'{}') ".format(fileformat))
            return flask.redirect(flask.url_for("input_data"))


        with open(token.info_path, 'r') as f:
            infodata = f.read()
            output.append(infodata)
        
        tkn_path = token.get_path()

        resp = flask.make_response(flask.render_template(
            "user_templates/c2m-infopage.html",
            output_lines=output,
            infodata=infodata,
            enumerate=enumerate, len=len,  # why is this needed?
            #token_path=tkn_path.replace('/','_'), #blueprint.url_for('process_structure','analysis', token=tkn_path.replace('/','_')),
            struct_name=token.refcode,
        ))
        resp.set_cookie("token_path",tkn_path,  secure=False,httponly=True,samesite='Strict') #secure must be True for final version
        return resp

    elif fileformat =="unknown":
        flask.flash("Please select a file format.")
        return flask.redirect(flask.url_for("input_data"))

    else:

        flask.flash("Please upload a file in format '{}' with the proper extension (.'{}') ".format(fileformat, fileformat))
        return flask.redirect(flask.url_for("input_data"))



@blueprint.route("/analysis", methods=["GET"])
def process_structure_analysis():

    output = Capturing()
    try:
        #tkn_path = flask.request.args['token'].replace('_','/')
        tkn_path = flask.request.cookies.get('token_path', None)
        if tkn_path is None:
            raise ValueError("not token?")
        
        token = Token.from_path(tkn_path)
        if token is None:
            raise ValueError("session expired")
        token.keepalive()
        
        #############################
        # STEP 3: Run cell2mol on infofile
        #############################
        #def run_cell2mol(run=False, out=None):

        if os.path.exists(token.info_path):
            with open(token.info_path, 'r') as f:
                infodata = f.read()
                output.append(infodata)

            with output as _out:
                cell = cell2mol(token.info_path, token.refcode, token.get_path(), 3)
                #with output as outt:
                #    print(cell, dir(cell), type(cell))
                #raise RuntimeError("unnamedrte")
                save_cell(cell, 'gmol', token.get_path())
                savemolecules(cell.moleclist, token.get_path(), 'xyz')
                savemolecules(cell.moleclist, token.get_path(), 'gmol')
            # This line removes cell2mol output from printout, it should be extend not redefine
            #output += [f"For input {token.input_path}"]
            celldata = printing_text(cell, Capturing())
            with open(token.analysis_path, 'w') as f:
                for line in celldata:
                    f.write(line+"\n")
        else:
            raise ValueError("plz")
            #output.extend([f"Please, wait until cell2info has finished for this input. Could not find {token.info_path}."])

        #token.remove()
        resp = flask.make_response(flask.render_template(
            "user_templates/c2m-analysis.html",
            output_lines=output,
            infodata=infodata.strip(),
            celldata='\n'.join(celldata).strip(),
            enumerate=enumerate, len=len,  # why TF is this needed?????
            #token_path=tkn_path.replace('/','_'), #blueprint.url_for('process_structure','analysis', token=tkn_path.replace('/','_')),
            struct_name=token.refcode,
        ))
        return resp
        #############################
        # STEP 4: Display structures
        #############################    
        #def display_mol(run=False, out=None):

        #with output as _out:
        #    path = find_gmol(run=run)
        #    cell = pickle.load(open(path, "rb"))
        #with out:
        #    printing_structure_cell(cell)
    
    except Exception as err:
        msg = "Failure…"
        output += traceback.format_tb(err.__traceback__)
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
    




@blueprint.route("/view-gmol", methods=["GET"])
def process_structure_view():

    output = Capturing()
    try:
        #tkn_path = flask.request.args['token'].replace('_','/')
        tkn_path = flask.request.cookies.get('token_path', None)
        if tkn_path is None:
            raise ValueError("not token?")
        
        token = Token.from_path(tkn_path)
        if token is None:
            raise ValueError("session expired")
        token.keepalive()
        
        #############################
        # STEP 3: Run cell2mol on infofile
        #############################
        #def run_cell2mol(run=False, out=None):

        if os.path.exists(token.info_path):
            labels, pos, lfracs, fracs, cellvec, cellparam = readinfo(token.info_path)
            with open(token.info_path, 'r') as f:
                infodata = f.read()
                output.append(infodata)
            with open(token.analysis_path, 'r') as f:
                celldata = f.read()
                output.append(celldata)
            with open(token.cell_path, 'rb') as f:
                cell = pickle.load(f)


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
            

            #ucellparams, xyzdata = cell_to_string_xyz(cell, cmp_lut)
            ucellparams, xyzdata = cell_to_string_xyz(cell, cmp_lut)

            #string used by jsmol to define the molecules/complexes, and the connectivity respectively
            jmol_list_pos = molecules_list(cell)
            jmolCon = bond_order_connectivity(cell)
            jmol_list_species =species_list(cell) 

        else:
            raise ValueError("Please wait until cell2info finishes.")
            #output.extend([f"Please, wait until cell2info has finished for this input. Could not find {token.info_path}."])

        #token.remove()
        resp = flask.make_response(flask.render_template(
            "user_templates/c2m-view.html",
            output_lines=output,
            infodata=infodata.strip(),
            celldata=celldata,
            ucellparams=ucellparams,
            compound_data=compound_data,
            xyzdata=xyzdata,
            labels=labels,
            pos=pos,
            cellvec=cellvec,
            cellparam=cellparam,
            jmol_list_pos=jmol_list_pos,
            jmol_list_species = jmol_list_species,
            jmolCon = jmolCon,
            totmol = len(cell.moleclist),
            enumerate=enumerate, len=len, zip=zip, # needed
            #token_path=tkn_path.replace('/','_'), #blueprint.url_for('process_structure','analysis', token=tkn_path.replace('/','_')),
            struct_name=token.refcode,
        ))
        return resp
        #############################
        # STEP 4: Display structures
        #############################    
        #def display_mol(run=False, out=None):

        #with output as _out:
        #    path = find_gmol(run=run)
        #    cell = pickle.load(open(path, "rb"))
        #with out:
        #    printing_structure_cell(cell)
    
    except Exception as err:
        msg = "Failure…"
        output += traceback.format_tb(err.__traceback__)
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

@blueprint.route("/download-gmol", methods=["GET"])
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
@blueprint.route("/download-xyzselected", methods=["GET", "POST"])
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

