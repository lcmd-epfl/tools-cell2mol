import os, time
import threading
import tempfile
import math
from ase.data import atomic_numbers


class InvalidXYZInput(ValueError):
    pass


def normalize_xyz_content(content):
    lines = content.splitlines()
    while lines and not lines[0].strip():
        lines.pop(0)

    if not lines:
        raise InvalidXYZInput("The XYZ file is empty")

    try:
        natoms = int(lines[0].strip())
    except ValueError as exc:
        raise InvalidXYZInput("The first XYZ line must contain the number of atoms") from exc

    if natoms < 1:
        raise InvalidXYZInput("The XYZ file must contain at least one atom")
    if len(lines) < 2:
        raise InvalidXYZInput("The XYZ file must include a second header line (the comment, which can be empty)")

    comment = lines[1]
    atom_lines = [(i, line) for i, line in enumerate(lines[2:], start=3) if line.strip()]
    if len(atom_lines) != natoms:
        raise InvalidXYZInput(
            f"The XYZ file declares {natoms} atoms, but contains {len(atom_lines)} "
            "coordinate records. Check the atom count and any extra frames."
        )

    for number, line in atom_lines:
        fields = line.split()
        if len(fields) < 4:
            raise InvalidXYZInput(f"XYZ line {number}: expected an element symbol and three coordinates")
        if atomic_numbers.get(fields[0], 0) == 0:
            raise InvalidXYZInput(f"XYZ line {number}: unknown element symbol '{fields[0]}'")
        try:
            xyz = [float(value) for value in fields[1:4]]
        except ValueError as exc:
            raise InvalidXYZInput(f"XYZ line {number}: coordinates must be numeric") from exc
        if not all(math.isfinite(value) for value in xyz):
            raise InvalidXYZInput(f"XYZ line {number}: coordinates must be finite numbers")

    return "\n".join([str(natoms), comment] + [line for _, line in atom_lines]) + "\n"

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
    def __init__(self, structurefile, fileformat):
        global TOKENS, TOKEN_LOCK, MONITOR_THR
        self.refcode = structurefile.filename.split('.')[-2]
        #refcode = cif_name.split("/")[-1].split(".")[-2]

        self.localdir = tempfile.TemporaryDirectory(dir='/tmp/cell2mol')
        self.info_path = os.path.join(self.localdir.name, 'info.txt')
        self.analysis_path = os.path.join(self.localdir.name, 'cell.txt')
        self.error_path = os.path.join(self.localdir.name, 'err.txt')

        if fileformat == "xyz-ase":
            self.input_path = os.path.join(self.localdir.name, 'input.xyz')
        else:
            self.input_path = os.path.join(self.localdir.name, 'input.cif')

        #Do not download gmol, but json 
        #self.cell_path = os.path.join(self.localdir.name, f'Cell_{self.refcode:s}.gmol')
        self.cell_path = os.path.join(self.localdir.name, f'{self.refcode:s}.json')

        self.last_alive = time.monotonic()

        if fileformat == "xyz-ase":
            original_content = structurefile.read()
            try:
                content = original_content.decode("utf-8-sig")
            except UnicodeDecodeError as exc:
                raise InvalidXYZInput("The XYZ file must use UTF-8 text encoding") from exc
            normalized_content = normalize_xyz_content(content)
            with open(os.path.join(self.localdir.name, 'input_original.xyz'), 'wb') as f:
                f.write(original_content)
            with open(self.input_path, 'w') as f:
                f.write(normalized_content)
        else:
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
