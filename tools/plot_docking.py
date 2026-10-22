#!/usr/bin/env python3
"""Validate the saved 1IEP example and render its output-coordinate overlay.

Run from any directory: python tools/plot_docking.py
Requires numpy and matplotlib. No docking engine or research archive is needed.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/fastdock-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
BLUE = "#346b99"
ORANGE = "#e39343"
INK = "#26343e"
MUTED = "#62717b"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_pdbqt(path, model=None):
    atoms, energy = [], None
    active = model is None
    for line in path.read_text().splitlines():
        if line.startswith("MODEL") and model is not None:
            active = int(line.split()[1]) == model
        if line.startswith("ENDMDL") and active and model is not None:
            break
        if not active:
            continue
        if line.startswith("REMARK VINA RESULT:"):
            energy = float(line.split()[3])
        if line[:6].strip() in ("ATOM", "HETATM"):
            atoms.append({
                "serial": int(line[6:11]), "name": line[12:16].strip(),
                "residue": line[17:27], "type": line[77:].strip(),
                "chain": line[21], "resnum": int(line[22:26]),
                "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            })
    return atoms, energy


def load_example():
    record = json.loads((ROOT / "benchmarks/docking-figure.json").read_text())
    for relative, expected in record["inputs"].items():
        assert digest(ROOT / relative) == expected, f"Input changed: {relative}"
    top_path = ROOT / record["topology"]["path"]
    assert digest(top_path) == record["topology"]["sha256"], "Topology changed"
    topology = json.loads(top_path.read_text())
    candidate = ROOT / record["fastdock"]["path"]
    assert digest(candidate) == record["fastdock"]["output_sha256"], "Output changed"
    assert record["original"]["output_sha256"] == digest(candidate)
    atoms, energy = read_pdbqt(candidate, record["model"])
    assert len(atoms) == record["validation"]["fastdock_atom_count"]
    assert len(record["original"]["model1_xyz"]) == record["validation"]["original_atom_count"]
    assert sum(a["type"] not in ("H", "HD", "HS") for a in atoms) == record["validation"]["heavy_atom_count"]
    assert sum(line.startswith("MODEL") for line in candidate.read_text().splitlines()) == record["validation"]["model_count"]
    input_atoms, _ = read_pdbqt(ROOT / "examples/basic/1iep_ligand.pdbqt")
    for actual, original_input, expected in zip(atoms, input_atoms, topology["atoms"], strict=True):
        for key in ("serial", "name", "residue", "type"):
            assert actual[key] == original_input[key] == expected[key], (key, actual)
    original = np.array(record["original"]["model1_xyz"], dtype=float)
    fastdock = np.array([a["xyz"] for a in atoms], dtype=float)
    assert original.shape == fastdock.shape
    rmsd = np.sqrt(np.mean(np.sum((original - fastdock) ** 2, axis=1)))
    maximum = np.max(np.abs(original - fastdock))
    assert rmsd == record["validation"]["all_atom_output_coordinate_rmsd_angstrom"]
    assert maximum == record["validation"]["maximum_output_coordinate_difference_angstrom"]
    assert energy == record["original"]["model1_energy_kcal_mol"]
    assert energy == record["fastdock"]["model1_energy_kcal_mol"]
    for bond in topology["bonds"]:
        assert 0 <= bond["a"] < len(atoms) and 0 <= bond["b"] < len(atoms)
        assert bond["order"] in (1, 2, 3)
    receptor, _ = read_pdbqt(ROOT / "examples/basic/1iep_receptor.pdbqt")
    print(f"Verified MODEL {record['model']}: {len(atoms)} atoms, "
          f"RMSD {rmsd:.3f} Å, max difference {maximum:.3f} Å, "
          f"Vina score {energy:.3f} kcal/mol; no alignment")
    return record, topology, atoms, receptor, original, fastdock, rmsd, energy


def camera(coords):
    # Principal axes expose the long ligand; the tilt keeps its 3D conformation visible.
    _, _, axes = np.linalg.svd(coords - coords.mean(axis=0))
    for row in axes:
        if row[np.argmax(np.abs(row))] < 0:
            row *= -1
    tilt = np.deg2rad(18)
    basis = np.array([axes[0], np.cos(tilt) * axes[1] + np.sin(tilt) * axes[2],
                      -np.sin(tilt) * axes[1] + np.cos(tilt) * axes[2]])
    return basis


def trace_segments(receptor, projected):
    ca = [i for i, a in enumerate(receptor) if a["name"] == "CA"]
    segments = []
    for i, j in zip(ca[:-1], ca[1:]):
        a, b = receptor[i], receptor[j]
        if a["chain"] == b["chain"] and b["resnum"] == a["resnum"] + 1:
            if np.linalg.norm(np.array(a["xyz"]) - b["xyz"]) < 4.5:
                segments.append(projected[[i, j], :2])
    return segments


def draw_ligand(ax, coords, atoms, bonds, color, width, atoms_size, zorder, detail):
    segments = []
    for bond in bonds:
        i, j = bond["a"], bond["b"]
        if atoms[i]["type"] in ("H", "HD", "HS") or atoms[j]["type"] in ("H", "HD", "HS"):
            continue
        start, end = coords[i, :2], coords[j, :2]
        direction = end - start
        perpendicular = np.array([-direction[1], direction[0]])
        perpendicular /= np.linalg.norm(perpendicular)
        order = bond["order"] if detail else 1
        offsets = [0] if order == 1 else ([-0.14, 0.14] if order == 2 else [-0.22, 0, 0.22])
        for offset in offsets:
            segments.append([start + offset * perpendicular, end + offset * perpendicular])
    ax.add_collection(LineCollection(segments, colors=color, linewidths=width,
                                     capstyle="round", zorder=zorder))
    heavy = [i for i, a in enumerate(atoms) if a["type"] not in ("H", "HD", "HS")]
    ax.scatter(coords[heavy, 0], coords[heavy, 1], s=atoms_size,
               color=color, linewidths=0, zorder=zorder + 0.1)


def scale_bar(ax, length, x, y):
    ax.plot([x, x + length], [y, y], color=MUTED, lw=1.5, clip_on=False)
    ax.plot([x, x], [y - .15, y + .15], color=MUTED, lw=1)
    ax.plot([x + length, x + length], [y - .15, y + .15], color=MUTED, lw=1)
    ax.text(x + length / 2, y - .45, f"{length} Å", ha="center", va="top", fontsize=9, color=MUTED)


def plot(output_dir):
    record, topology, atoms, receptor, original, fastdock, rmsd, energy = load_example()
    heavy = np.array([a["type"] not in ("H", "HD", "HS") for a in atoms])
    receptor_heavy = np.array([a["type"] not in ("H", "HD", "HS") for a in receptor])
    center = fastdock[heavy].mean(axis=0)
    basis = camera(fastdock[heavy])
    original_view = (original - center) @ basis.T
    fastdock_view = (fastdock - center) @ basis.T
    receptor_xyz = np.array([a["xyz"] for a in receptor])
    receptor_view = (receptor_xyz - center) @ basis.T
    nearest = np.min(np.linalg.norm(receptor_xyz[:, None, :] - fastdock[heavy][None, :, :], axis=2), axis=1)
    pocket = receptor_heavy & (nearest <= record["rendering"]["pocket_cutoff_angstrom"])
    traces = trace_segments(receptor, receptor_view)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.fonttype": "none", "svg.hashsalt": "fastdock-docking-v1",
                         "axes.unicode_minus": True})
    fig = plt.figure(figsize=(15, 8.5), facecolor="white")
    fig.text(.045, .945, "One docked pose. Identical output coordinates.",
             fontsize=23, fontweight="bold", color=INK)
    fig.text(.045, .904, "1IEP example · independently executed original and FastDock Vina · orthographic projection",
             fontsize=12, color=MUTED)
    left = fig.add_axes([.045, .265, .37, .55])
    right = fig.add_axes([.465, .265, .49, .55])
    for ax in (left, right):
        ax.set_aspect("equal")
        ax.axis("off")

    left.scatter(receptor_view[receptor_heavy, 0], receptor_view[receptor_heavy, 1],
                 s=8, color="#adb8bf", alpha=.10, edgecolors="none", rasterized=True)
    left.add_collection(LineCollection(traces, colors="#9eabb4", linewidths=1.4, alpha=.7))
    draw_ligand(left, original_view, atoms, topology["bonds"], BLUE, 3.2, 12, 4, False)
    draw_ligand(left, fastdock_view, atoms, topology["bonds"], ORANGE, 1.4, 3, 5, False)
    context = receptor_view[receptor_heavy, :2]
    lower, upper = context.min(axis=0) - 3, context.max(axis=0) + 3
    left.set_xlim(lower[0], upper[0]); left.set_ylim(lower[1], upper[1])
    extent = fastdock_view[heavy, :2]
    lo, hi = extent.min(axis=0) - 2.0, extent.max(axis=0) + 2.0
    left.add_patch(Rectangle(lo, *(hi - lo), fill=False, edgecolor=BLUE, lw=1.1,
                             linestyle=(0, (4, 3)), zorder=6))
    left.text(.02, .99, "A   RECEPTOR CONTEXT", transform=left.transAxes,
              fontsize=11, color=INK, fontweight="bold", va="top")
    left.text(.02, -.06, "Gray: receptor heavy atoms + Cα backbone trace",
              transform=left.transAxes, fontsize=9, color=MUTED)

    right.scatter(receptor_view[pocket, 0], receptor_view[pocket, 1],
                  s=115, color="#d5dce0", alpha=.29, edgecolors="none", rasterized=True)
    right.add_collection(LineCollection(traces, colors="#c0c9cf", linewidths=2,
                                        alpha=.30, zorder=1))
    draw_ligand(right, original_view, atoms, topology["bonds"], BLUE, 5.8, 50, 4, True)
    draw_ligand(right, fastdock_view, atoms, topology["bonds"], ORANGE, 2.6, 15, 5, True)
    for i, atom in enumerate(atoms):
        element = atom["type"].rstrip("A") if atom["type"] in ("NA", "OA") else atom["type"]
        if element in ("N", "O"):
            x, y = fastdock_view[i, :2]
            right.text(x, y, element, fontsize=8.5, ha="center", va="center", color=INK,
                       fontweight="bold", zorder=8,
                       bbox={"boxstyle": "circle,pad=.16", "facecolor": "white", "edgecolor": BLUE, "linewidth": .8})
    right.set_xlim(lo[0] - .4, hi[0] + .4)
    right.set_ylim(lo[1] - 1.8, hi[1] + 2.0)
    right.text(0, 1.03, "B   POCKET CLOSE-UP · MODEL 1", transform=right.transAxes,
               fontsize=11, color=INK, fontweight="bold", va="bottom")
    legend = [Line2D([0], [0], color=BLUE, lw=7, solid_capstyle="round"),
              Line2D([0], [0], color=ORANGE, lw=3, solid_capstyle="round")]
    right.legend(legend, ["Original · outer sticks", "FastDock · inner sticks"],
                  loc="upper right", bbox_to_anchor=(1.02, 1.02), frameon=False,
                  labelcolor=INK, fontsize=9, handlelength=1.8)
    scale_bar(right, 5, lo[0] + .1, lo[1] - .6)
    right.text(0, -.06, f"{heavy.sum()} heavy atoms shown · N/O labeled", transform=right.transAxes,
               fontsize=9, color=MUTED)
    right.text(1, -.06, "Identical coordinates · no alignment", transform=right.transAxes,
               ha="right", fontsize=9, color=MUTED)

    fig.add_artist(Line2D([.045, .955], [.226, .226], transform=fig.transFigure,
                          color="#dce3e7", linewidth=1))
    cards = [(f"{rmsd:.3f} Å", f"Output-coordinate RMSD · all {len(atoms)} atoms"),
             (f"{energy:.3f} kcal/mol", "Vina score · both MODEL 1 poses"),
             ("Byte-for-byte identical", f"Both complete PDBQT outputs · {record['validation']['model_count']} models")]
    for x, (number, caption) in zip([.045, .375, .695], cards):
        fig.text(x, .174, number, fontsize=19 if x < .6 else 17, color=BLUE, fontweight="bold")
        fig.text(x, .141, caption, fontsize=9, color=MUTED)
    p = record["parameters"]
    fig.text(.045, .083, f"PDBQT MODEL 1  ·  scoring {p['scoring']}  ·  seed {p['seed']}  ·  "
             f"exhaustiveness {p['exhaustiveness']}  ·  CPU {p['cpu']}", fontsize=10, color=INK)
    fig.text(.045, .045, "Coordinates are rounded to 0.001 Å in PDBQT. "
             "Docking score is not experimental affinity. No interaction lines are inferred.",
             fontsize=9, color=MUTED)
    output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "svg"):
        target = output_dir / f"docking-overlay.{suffix}"
        metadata = {"Creator": "FastDock Vina tools/plot_docking.py"}
        if suffix == "svg":
            metadata["Date"] = None
        fig.savefig(target, dpi=180, facecolor="white", metadata=metadata)
        print(target.relative_to(ROOT) if target.is_relative_to(ROOT) else target)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "docs/figures")
    plot(parser.parse_args().output_dir)
