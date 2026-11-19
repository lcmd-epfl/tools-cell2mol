from cell2mol.utils import config


class Cell2molNoResultError(RuntimeError):
    """The interpreter stopped without returning a displayable result."""

    def __init__(self, structure, mode):
        super().__init__(
            f"cell2mol returned no usable {mode} result for {structure}."
        )



def describe_processing_error(error):
    if isinstance(error, Cell2molNoResultError):
        return {
            "title": "No result returned by cell2mol",
            "summary": "The chemical interpretation did not produce a result that can be displayed or downloaded.",
            "detail": str(error),
            "suggestion": (
                "See Calculation Log for the warnings and errors reported by cell2mol. "
                "If CIF connectivity was enabled, you may also compare the result with "
                "distance-based connectivity. Do not assume a charge or bond assignment succeeded."
            ),
            "partial": None,
            "code": None,
        }

    code = getattr(error, "exit_code", None)
    if code is None:
        if isinstance(error, MemoryError):
            code = config.ERR_MEMORY
        elif isinstance(error, TimeoutError):
            code = config.ERR_TIMEOUT
        elif isinstance(error, (FileNotFoundError, IsADirectoryError, PermissionError)):
            code = config.ERR_INPUT

    categories = {
        config.ERR_INPUT: (
            "Input could not be processed",
            "cell2mol could not access or accept the input structure.",
            "Check the file format, file contents, and uploaded structure.",
        ),
        config.ERR_ASE_PARSE: (
            "Structure could not be read",
            "The structure reader could not interpret the input coordinates.",
            "Check the input format and the coordinate records.",
        ),
        config.ERR_MEMORY: (
            "Calculation stopped: memory limit",
            "The calculation exceeded the available memory.",
            "Try a smaller structure or check the resource limits with the service administrator.",
        ),
        config.ERR_TIMEOUT: (
            "Calculation stopped: time limit",
            "The calculation did not finish within the allowed time.",
            "Try a less demanding structure or check the time limit with the service administrator.",
        ),
        config.ERR_GENERAL: (
            "Calculation interrupted",
            "An error interrupted cell2mol before the calculation could finish.",
            "Open Calculation Log for the detailed error and traceback.",
        ),
    }
    title, summary, suggestion = categories.get(code, (
        "Calculation could not be completed",
        "An unexpected error interrupted the calculation.",
        "Open Calculation Log for the detailed error and traceback.",
    ))

    partial = getattr(error, "result", None)
    if partial is not None:
        molecules = getattr(partial, "moleclist", None)
        if molecules is not None:
            partial_message = (
                f"cell2mol returned a partial result containing {len(molecules)} "
                f"molecule(s). This does not mean charge assignment succeeded."
            )
        else:
            partial_message = "cell2mol returned a partial result before stopping."
    else:
        partial_message = None

    return {
        "title": title,
        "summary": summary,
        "detail": str(error) or type(error).__name__,
        "suggestion": suggestion,
        "partial": partial_message,
        "code": code,
    }


def molecule_status_message(molset, is_resolved):
    from cell2mol.write_results import get_molecule_error_messages

    codes = list(getattr(molset, "error_cases_all", None) or [])
    primary = getattr(molset, "error_case", None)
    if primary not in (None, 0) and primary not in codes:
        codes.insert(0, primary)
    codes = [code for code in codes if code not in (None, 0)]
    prefix = "Molecule Status (mode=interpretation): "

    if codes:
        details = get_molecule_error_messages(codes, primary or 0)
        if 6 in codes:
            return prefix + "Multiple charge distributions satisfy the input charge; charges remain ambiguous. " + details
        if 7 in codes:
            return prefix + "No candidate charge distribution satisfies the input charge. " + details
        if 5 in codes:
            return prefix + "Plausible charge states could not be generated for every species. " + details
        if any(code in codes for code in (1, 2, 3, 4)):
            return prefix + "Hydrogen-related structural issues were detected. " + details
        if 8 in codes or 9 in codes:
            return prefix + "Charge or bond assignment was not completed. " + details
        if 10 in codes:
            return prefix + "Spin assignment could not be completed. " + details
        return prefix + "cell2mol reported a processing problem. " + details

    if primary == 0 and is_resolved:
        return prefix + "Chemical interpretation completed without reported errors."
    if primary == 0:
        return prefix + "No error code was reported, but some chemical states remain unresolved. Check Calculation Log."
    return prefix + "No final processing status was reported. Some results may be incomplete; check Calculation Log."
