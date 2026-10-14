"""Pure-Python inclusive intervals for the parallel anisotropic search."""
from operator import index


def inclusive_chunks(lower, count, parts):
    """Partition count integers starting at lower into nonempty intervals."""
    lower, count, parts = index(lower), index(count), index(parts)
    if lower < 1 or count < 0 or parts < 1:
        raise ValueError("require lower >= 1, count >= 0 and parts >= 1")
    if count == 0:
        return []
    parts = min(parts, count)
    width, extra = divmod(count, parts)
    result = []
    start = lower
    for i in range(parts):
        stop = start + width + (i < extra) - 1
        result.append((start, stop))
        start = stop + 1
    return result
