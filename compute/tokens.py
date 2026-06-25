import os, time
import threading
import tempfile

def monitoring():
    global MONITOR_THR

    try:
        while True:
            time.sleep(MONITOR_PERIOD)

            now = time.time()
            with TOKEN_LOCK:
                expired_paths = [
                    path for path, tkn in list(TOKENS.items())
                    if (now - tkn.last_alive) > MAX_TOKEN_AGE
                ]
                for path in expired_paths:
                    TOKENS.pop(path, None)

                if not TOKENS:
                    MONITOR_THR = None
                    break

    except Exception:
        # Don’t kill Apache worker because a background cleanup thread crashed
        traceback.print_exc()
        MONITOR_THR = None


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

        self.last_alive = time.time()


        if fileformat == "xyz-ase":
            self.input_path = os.path.join(self.localdir.name, 'input.xyz')
        else:
            self.input_path = os.path.join(self.localdir.name, 'input.cif')

        self.cell_path = os.path.join(self.localdir.name, f'Cell_{self.refcode:s}.cell')
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
