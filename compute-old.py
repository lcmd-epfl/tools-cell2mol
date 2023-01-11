#!/usr/bin/env python3
import os, io, sys, logging
sys.path.append(os.path.join(os.path.split(__file__)[0], 'cell2mol'))
#import json
import time
import threading
import tempfile
import flask
from tools_barebone.structure_importers import get_structure_tuple, UnknownFormatError
import cell2mol
from cell2mol.cif2info import cif_2_info
from cell2mol.c2m_module import save_cell, cell2mol

logger = logging.getLogger("tools-app")

os.makedirs('/tmp/cell2mol', exist_ok=True)

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = io.StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout



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




def monitoring():
    """expire all sessions that haven't recieved a keepalive in 5min"""
    global TOKENS, TOKEN_LOCK, MONITOR_THR,thr_iters
    cond = True
    while cond:
        try:
            thr_iters +=1
            TOKEN_LOCK.acquire()
            now = time.monotonic()
            for path, tkn in TOKENS.items():
                if tkn.last_alive + 300 < now:
                    del TOKENS[path]
            cond = len(TOKENS)
        except:
            thr_iters -= 1000
            raise
        finally:
            TOKEN_LOCK.release()
        if cond:
            time.sleep(60)
    MONTIOR_THR = None


TOKEN_LOCK = threading.Lock()
TOKENS = {}
MONITOR_THR = None
thr_iters = 0
class Token:
    def __init__(self, structurefile):
        global TOKENS, TOKEN_LOCK, MONITOR_THR
        self.refcode = structurefile.filename.split('.')[-2]
        #refcode = cif_name.split("/")[-1].split(".")[-2]

        self.localdir = tempfile.TemporaryDirectory(dir='/tmp/cell2mol')
        self.info_path = os.path.join(self.localdir.name, 'info.txt')
        self.error_path = os.path.join(self.localdir.name, 'err.txt')
        self.input_path = os.path.join(self.localdir.name, 'input.cif')
        self.cell_path = os.path.join(self.localdir.name, f'Cell_{self.refcode:s}.gmol')
        self.last_alive = time.monotonic()

        with open(self.input_path,'w') as f:
            f.write(structurefile.read().decode("utf-8"))

        try:
            TOKEN_LOCK.acquire()
            TOKENS[self.localdir.name] = self
        finally:
            TOKEN_LOCK.release()

        if MONITOR_THR is None:
            MONITOR_THR = threading.Thread(target=monitoring)
            MONITOR_THR.start()
            

    def __del__(self):
        self.localdir.cleanup()

    def remove(self):
        global TOKENS, TOKEN_LOCK
        try:
            TOKEN_LOCK.acquire()
            del TOKENS[self.localdir.name]
        finally:
            TOKEN_LOCK.release()
    def keepalive(self):
        self.last_alive = time.monotonic()

    @staticmethod
    def from_path(path):
        global TOKENS, TOKEN_LOCK
        try:
            TOKEN_LOCK.acquire()
            res = TOKENS.get(path, None)
        finally:
            TOKEN_LOCK.release()
        return res
    def get_path(self):
        return self.localdir.name



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
                save_cell(cell, 'gmol', token.get_path())
            # This line removes cell2mol output from printout, it should be extend not redefine
            #output += [f"For input {token.input_path}"]
            celldata = printing_text(cell, Capturing())
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
        output.append(str(os.listdir(tkn_path)))
        output.append(str((TOKENS, TOKEN_LOCK, MONITOR_THR, thr_iters)))
        output.append(repr(time.monotonic() - token.last_alive ))
        output.append(os.path.relpath(token.cell_path, '/tmp/cell2mol'))
        output.append(repr(os.stat(token.cell_path)))
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


@blueprint.route("/keepalive", methods=["GET"])
def session_keepalive():

    tkn_path = flask.request.cookies.get('token_path', None)
    if tkn_path is None:
        return flask.make_response("", 404)

    token = Token.from_path(tkn_path)
    if token is None:
        return flask.make_response("", 404)
    else:
        token.keepalive()
        return flask.make_response("", 200)


@blueprint.route("/view-gmol", methods=["GET"])
def process_structure_view():

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
        output.append(str(os.listdir(tkn_path)))
        output.append(str((TOKENS, TOKEN_LOCK, MONITOR_THR, thr_iters)))
        output.append(repr(time.monotonic() - token.last_alive ))
        token.keepalive()
        raise NotImplementedError()
     
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
