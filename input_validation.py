"""CSV parsing and spreadsheet-safe export helpers."""

def parse_csv_bool(value):
    """Parse common CSV boolean representations without truthy-string mistakes."""
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if normalized in {"0", "false", "f", "no", "n", "off", ""}:
        return False
    raise ValueError(f"Invalid boolean in CSV: {value!r}")


def safe_csv_cell(value):
    """Keep untrusted CSV text from becoming a spreadsheet formula on opening."""
    if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")):
        return "'" + value
    return value
