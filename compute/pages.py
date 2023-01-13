import traceback
import pickle
from tools_barebone.structure_importers import get_structure_tuple, UnknownFormatError
import flask
from .interface import *
from .tokens import monitoring, Token

blueprint = flask.Blueprint("compute", __name__, url_prefix="/compute")
    
@blueprint.route("/process_structure/", methods=["POST"])
def process_structure_init():
    """Example view to process a crystal structure."""

    # check if the post request has the file part,
    # otherwise redirect to first page

    output = Capturing()
    try:
        if "structurefile" not in flask.request.files:
            # This will redirect the user to the selection page,
            # that is called `input_data` in tools-barebone
            return flask.redirect(flask.url_for("input_data"))

        
        # Get structure, file format, file content, and form data
        # (needed for additional information, e.g. cell in the case
        # of a XYZ file)
        fileformat = flask.request.form.get("fileformat", "unknown")
        form_data = dict(flask.request.form)
        if fileformat not in ('cif-pymatgen','cif-ase'):
            flask.flash("er… well we will interpret that a a cif file anyway >:)")
        structurefile = flask.request.files["structurefile"]

        token = Token(structurefile)        
        output += [
            "structure_file: "+ repr(dir(structurefile)),
            "form_data: "+ repr(form_data),
            "file_format: "+ repr(fileformat),
            "code: " + token.refcode,
            'req: '+ repr((flask.request, dir(flask.request))),
        ]

        
        #############################
        # STEP 2: Run cell2info
        #############################
        #def run_cell2info(run=False, out=None):
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

        with open(token.info_path, 'r') as f:
            infodata = f.read()
            output.append(infodata)
        
        tkn_path = token.get_path()

        resp = flask.make_response(flask.render_template(
            "user_templates/c2m-infopage.html",
            output_lines=output,
            infodata=infodata,
            enumerate=enumerate, len=len,  # why TF is this needed?????
            #token_path=tkn_path.replace('/','_'), #blueprint.url_for('process_structure','analysis', token=tkn_path.replace('/','_')),
            struct_name=token.refcode,
        ))
        resp.set_cookie("token_path",tkn_path,  secure=True,httponly=True,samesite='Strict')
        return resp
        #return flask.redirect(flask.url_for('process_structure/info', token=tkn_path))
    except Exception as err:
        output += traceback.format_tb(err.__traceback__)
        output.append(repr(err))
        return flask.render_template(
            "user_templates/c2m-debug.html", msg="Failure…", output_lines=output,
        )        
    except:
        return flask.render_template(
            "user_templates/c2m-debug.html", msg="Unknown Failure…", output_lines=output,
        )        


@blueprint.route("/analysis", methods=["GET"])
def process_structure_analysis():

    output = Capturing()
    try:
        #tkn_path = flask.request.args['token'].replace('_','/')
        tkn_path = flask.request.cookies.get('token_path', None)
        if tkn_path is None:
            raise ValueError("not token? :( ")
        
        token = Token.from_path(tkn_path)
        if token is None:
            raise ValueError("session expired :( ")
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
    

@blueprint.route("/download-gmol", methods=["GET"])
def process_structure_download():

    output = Capturing()
    try:
        #tkn_path = flask.request.args['token'].replace('_','/')
        tkn_path = flask.request.cookies.get('token_path', None)
        if tkn_path is None:
            raise ValueError("no token? :( ")
        token = Token.from_path(tkn_path)
        if token is None:
            raise ValueError("session expired :( ")
        output.append(token.cell_path)
        #output.append(str(os.listdir(tkn_path)))
        #output.append(str((TOKENS, TOKEN_LOCK, MONITOR_THR, thr_iters)))
        #output.append(repr(time.monotonic() - token.last_alive ))
        #output.append(os.path.relpath(token.cell_path, '/tmp/cell2mol'))
        #output.append(repr(os.stat(token.cell_path)))
        token.keepalive()
                
#        res= flask.send_from_directory(
#            '/tmp/cell2mol',
#            os.path.relpath(token.cell_path, '/tmp/cell2mol'),
#        res = flask.send_file(
#            token.cell_path,
            #etag=True,
            #as_attachment=True,
            #            conditional=True,
#        )

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


@blueprint.route("/view-gmol", methods=["GET"])
def process_structure_view():

    output = Capturing()
    try:
        #tkn_path = flask.request.args['token'].replace('_','/')
        tkn_path = flask.request.cookies.get('token_path', None)
        if tkn_path is None:
            raise ValueError("not token? :( ")
        
        token = Token.from_path(tkn_path)
        if token is None:
            raise ValueError("session expired :( ")
        token.keepalive()
        
        #############################
        # STEP 3: Run cell2mol on infofile
        #############################
        #def run_cell2mol(run=False, out=None):

        if os.path.exists(token.info_path):
            with open(token.info_path, 'r') as f:
                infodata = f.read()
                output.append(infodata)
            with open(token.analysis_path, 'r') as f:
                celldata = f.read()
                output.append(celldata)
            with open(token.cell_path, 'rb') as f:
                cell = pickle.load(f)

            cmp_lut = cell_cmp_lut(cell)
            #xsf = cell_to_string_xsf(cell, cmp_lut)
            ucellparams, xyzdata = cell_to_string_xyz(cell, cmp_lut)
            svgs = cell_to_svgs(cell, cmp_lut)
            

        else:
            raise ValueError("plz")
            #output.extend([f"Please, wait until cell2info has finished for this input. Could not find {token.info_path}."])

        #token.remove()
        resp = flask.make_response(flask.render_template(
            "user_templates/c2m-view.html",
            output_lines=output,
            infodata=infodata.strip(),
            celldata=celldata,
            #xsfdata=xsf,
            ucellparams=ucellparams,
            compound_names=cmp_lut.keys(),
            xyzdata=xyzdata,
            svg_list=svgs,
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

