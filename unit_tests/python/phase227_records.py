"""Independent validation of authoritative census rows and canonical omega records."""

import math
import struct
from pathlib import Path

CLASSES = {0: "FULL", 1: "BULK", 2: "ITLIM", 3: "STAG", 4: "DIV", 5: "INV"}


def census(path, expected_rays, aggregate):
    result = {}
    for line_number, line in enumerate(Path(path).read_text().splitlines(), 1):
        fields = line.split()
        if len(fields) != 18:
            raise ValueError(f"{path}:{line_number}: expected 18 census columns")
        try:
            refresh, ray, cls, reason, exit_code, bad = map(int, fields[:6])
            residual, impact = map(float, fields[6:8])
            ni = int(fields[8])
            for offset in (9, 13):
                values = [float(fields[offset]), float(fields[offset + 3])]
                if not all(map(math.isfinite, values)):
                    raise ValueError("nonfinite channel values")
                int(fields[offset + 1])
                int(fields[offset + 2])
            if fields[17] not in {"T", "F"}:
                raise ValueError("bad production mask")
            if cls not in CLASSES or ni < 0 or bad < 0:
                raise ValueError("invalid class/geometry/count")
            if not math.isfinite(residual) or not math.isfinite(impact):
                raise ValueError("nonfinite residual/geometry")
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: {error}") from error
        rows = result.setdefault(refresh, {})
        if ray in rows:
            raise ValueError(f"{path}: refresh {refresh}: duplicate ray {ray}")
        rows[ray] = {
            "class": CLASSES[cls],
            "reason": reason,
            "exit": exit_code,
            "bad": bad,
            "residual": residual,
            "impact": impact,
            "ni": ni,
            "production": fields[17] == "T",
        }
    if not result:
        raise ValueError(f"{path}: empty census")
    for refresh, rows in result.items():
        if set(rows) != set(expected_rays):
            raise ValueError(f"{path}: refresh {refresh}: missing/extra ray IDs")
        counts = {name: 0 for name in CLASSES.values()}
        for row in rows.values():
            if row["production"]:
                counts[row["class"]] += 1
        if counts != aggregate[refresh]:
            raise ValueError(f"{path}: refresh {refresh}: aggregate mismatch")
    return result


def normalize_omega(source, destination, ray, refresh, expected_counts):
    rows = []
    seen = [0, 0]
    for line_number, line in enumerate(Path(source).read_text().splitlines(), 1):
        fields = line.split()
        if len(fields) != 8:
            raise ValueError(f"{source}:{line_number}: expected 8 application columns")
        rf, ir, channel, direction, step = map(int, fields[:5])
        if ir != ray:
            raise ValueError(f"{source}:{line_number}: wrong ray {ir}")
        if rf != refresh:
            continue
        if channel not in (1, 2) or direction != 0:
            raise ValueError(f"{source}:{line_number}: channel/direction")
        proposed, applied = map(float, fields[5:7])
        bits = struct.pack(">d", applied).hex().upper()
        if not all(map(math.isfinite, (proposed, applied))):
            raise ValueError(f"{source}:{line_number}: nonfinite omega")
        if bits != fields[7].upper() or not 0.45 <= applied <= 1:
            raise ValueError(f"{source}:{line_number}: binary64/range failure")
        if step != seen[channel - 1] + 1:
            raise ValueError(f"{source}:{line_number}: pass order")
        seen[channel - 1] = step
        rows.append(line)
    if seen != list(expected_counts):
        raise ValueError(
            f"{source}: expected pass counts {expected_counts}, found {seen}"
        )
    with Path(destination).open("x", encoding="ascii") as stream:
        stream.write(f"KNOMG2 {ray} {refresh} {seen[0]} {seen[1]}\n")
        stream.write("\n".join(rows) + "\n")
