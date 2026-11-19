import contextvars
import functools
import logging
import os
import threading


_ACTIVE_LOG = contextvars.ContextVar("cell2mol_active_calculation", default=None)
_LOGGER_LOCK = threading.Lock()
_LOGGER_COUNT = 0
_PREVIOUS_LEVEL = None


class _CalculationFilter(logging.Filter):
    def __init__(self, session):
        super().__init__()
        self.session = session

    def filter(self, record):
        return _ACTIVE_LOG.get() is self.session


def calculation_log_path(token):
    return os.path.join(token.get_path(), "calculation.log")


def start_calculation_log(token, mode, fileformat, total_charge=None):
    global _LOGGER_COUNT, _PREVIOUS_LEVEL
    if _ACTIVE_LOG.get() is not None:
        raise RuntimeError("A calculation log is already active")
    handler = logging.FileHandler(calculation_log_path(token), encoding="utf-8")
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s - %(name)s - %(message)s"))
    session = {"handler": handler, "failed": False}
    handler.addFilter(_CalculationFilter(session))
    cell2mol_logger = logging.getLogger("cell2mol")
    web_logger = logging.getLogger("compute.calculation")
    web_logger.setLevel(logging.INFO)
    with _LOGGER_LOCK:
        if _LOGGER_COUNT == 0:
            _PREVIOUS_LEVEL = cell2mol_logger.level
            cell2mol_logger.setLevel(logging.DEBUG)
        _LOGGER_COUNT += 1
        cell2mol_logger.addHandler(handler)
        web_logger.addHandler(handler)
    session["context_token"] = _ACTIVE_LOG.set(session)
    web_logger.info("Starting calculation: mode=%s, format=%s, structure=%s", mode, fileformat, token.refcode)
    if fileformat == "xyz-ase":
        web_logger.info("Input total charge: %s", total_charge)


def mark_calculation_failure(error):
    session = _ACTIVE_LOG.get()
    if session is None:
        return
    session["failed"] = True
    logging.getLogger("compute.calculation").error(
        "Calculation failed: %s", error, exc_info=(type(error), error, error.__traceback__)
    )


def finish_calculation_log():
    global _LOGGER_COUNT, _PREVIOUS_LEVEL
    session = _ACTIVE_LOG.get()
    if session is None:
        return
    logging.getLogger("compute.calculation").info(
        "Calculation request finished: %s", "error" if session["failed"] else "completed"
    )
    handler = session["handler"]
    with _LOGGER_LOCK:
        logging.getLogger("cell2mol").removeHandler(handler)
        logging.getLogger("compute.calculation").removeHandler(handler)
        _LOGGER_COUNT -= 1
        if _LOGGER_COUNT == 0:
            logging.getLogger("cell2mol").setLevel(_PREVIOUS_LEVEL)
            _PREVIOUS_LEVEL = None
    _ACTIVE_LOG.reset(session["context_token"])
    handler.close()


def calculation_logged(view):
    @functools.wraps(view)
    def wrapped(*args, **kwargs):
        try:
            return view(*args, **kwargs)
        except Exception as error:
            mark_calculation_failure(error)
            raise
        finally:
            finish_calculation_log()
    return wrapped


def get_calculation_log_preview(token, max_bytes=40000, max_lines=150):
    path = calculation_log_path(token)
    if not os.path.isfile(path):
        return "No calculation log available."
    with open(path, "rb") as handle:
        handle.seek(0, os.SEEK_END)
        size = handle.tell()
        handle.seek(max(0, size - max_bytes))
        content = handle.read().decode("utf-8", errors="replace")
    lines = content.splitlines()
    truncated = size > max_bytes or len(lines) > max_lines
    preview = "\n".join(lines[-max_lines:])
    if truncated:
        return "… Earlier log messages omitted. Download the full log below.\n" + preview
    return preview or "No log messages recorded."
