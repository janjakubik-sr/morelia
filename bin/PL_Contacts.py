#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarize protein-ligand contacts in one or more directories."""

from collections import Counter
import glob
import os
import re
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import wx

import preferences as pref


CONTACT_FILES = (
    "PL-Contacts_HBond.dat",
    "PL-Contacts_Hydrophobic.dat",
    "PL-Contacts_Ionic.dat",
    "PL-Contacts_Pi-Cation.dat",
    "PL-Contacts_Pi-Pi.dat",
    "PL-Contacts_WaterBridge.dat",
)
CONTACT_LABELS = (
    "HBond",
    "Hydrophobic",
    "Ionic",
    "Pi-Cation",
    "Pi-Pi",
    "Water bridge",
)
GRACE_COLORS = (22, 28, 32, 24, 16, 41)


def show_error(message):
    """Display an error message."""
    dialog = wx.MessageDialog(None, message, "Error", wx.OK | wx.ICON_ERROR)
    dialog.ShowModal()
    dialog.Destroy()


def ask_for_value(message, title, default):
    """Ask for a value and return None if the dialog is cancelled."""
    dialog = wx.TextEntryDialog(None, message, title, default)
    value = dialog.GetValue() if dialog.ShowModal() == wx.ID_OK else None
    dialog.Destroy()
    return value


def find_directories(pattern):
    """Return the directories matching a user-supplied glob pattern."""
    expanded_pattern = os.path.expandvars(os.path.expanduser(pattern))
    matches = glob.glob(expanded_pattern, recursive=True)
    return sorted({os.path.abspath(path) for path in matches if os.path.isdir(path)})


def read_contact_counts(path):
    """Count occurrences of each residue number in a contact file."""
    residues = []
    with open(path, encoding="utf-8") as contact_file:
        for line in contact_file:
            fields = line.split()
            if len(fields) > 1 and fields[1].isdigit():
                residues.append(int(fields[1]))
    return Counter(residues)


def read_residue_names(path):
    """Map residue numbers to the identifiers stored in P_RMSF.dat."""
    residue_names = {}
    with open(path, encoding="utf-8") as rmsf_file:
        for line in rmsf_file:
            fields = line.split()
            if len(fields) < 3:
                continue
            match = re.search(r"_(\d+)$", fields[2])
            if match:
                residue_names[int(match.group(1))] = fields[2]
    return residue_names


def write_summary(path, rows):
    """Write the contact-count summary table."""
    with open(path, "w", encoding="utf-8") as summary:
        for row in rows:
            summary.write(" ".join(str(value) for value in row) + "\n")


def write_grace_files(directory, residues, fractions):
    """Write the Grace data files and drawing instructions."""
    indices = np.arange(len(residues))
    data_paths = []
    for index, values in enumerate(fractions):
        data_path = os.path.join(directory, f"s{index}.dat")
        np.savetxt(data_path, np.column_stack((indices, values)), fmt="%s")
        data_paths.append(data_path.replace('"', '\\"'))

    draw_path = os.path.join(directory, "PL_Contacts.draw")
    with open(draw_path, "w", encoding="utf-8") as draw:
        draw.write("page layout portrait \n")
        draw.write("view xmin 0.2 \n")
        draw.write("view xmax 0.8 \n")
        draw.write("view ymin 0.85 \n")
        draw.write("view ymax 1.25 \n")
        draw.write(" G0  type Chart \n")
        draw.write(" G0  stacked true \n")
        draw.write('  xaxis  label "Residue " \n')
        draw.write('  xaxis  label font "Helvetica" \n')
        draw.write("  xaxis  label char size 1.0 \n")
        draw.write('  xaxis  ticklabel font "Helvetica" \n')
        draw.write("  xaxis  ticklabel char size 0.8 \n")
        draw.write("  xaxis  ticklabel angle 45 \n")
        draw.write('  yaxis  label "Interactions fraction" \n')
        draw.write('  yaxis  label font "Helvetica" \n')
        draw.write("  yaxis  label char size 1.0 \n")
        draw.write('  yaxis  ticklabel font "Helvetica" \n')
        draw.write("  yaxis  ticklabel char size 0.8 \n")
        draw.write("  Legend  anchor 0.22 , 0.64 \n")
        draw.write('  Legend  font "Helvetica" \n')
        draw.write("  xaxis  tick off \n")
        draw.write("  xaxis  tick spec type both \n")
        draw.write(f"  xaxis  tick spec {len(residues)} \n")
        for index, residue in enumerate(residues):
            safe_residue = residue.replace('"', '\\"')
            draw.write(f"xaxis  tick major {index} ,{index} \n")
            draw.write(f'xaxis  ticklabel {index} ,"{safe_residue}" \n')

        for index, (data_path, label, color) in enumerate(
            zip(data_paths, CONTACT_LABELS, GRACE_COLORS)
        ):
            draw.write(f' read xy "{data_path}" \n')
            draw.write(f"    S{index} type bar \n")
            draw.write(f"    S{index} symbol size 1.500000 \n")
            draw.write(f"    S{index} symbol color {color} \n")
            draw.write(f"    S{index} symbol fill color {color} \n")
            draw.write(f"    S{index} symbol fill pattern 1 \n")
            draw.write(f"    S{index} symbol linestyle 0 \n")
            draw.write(f"    S{index} line type 0 \n")
            draw.write(f"    S{index} legend on \n")
            draw.write(f'    S{index} legend "{label}" \n')

        draw.write("world ymin 0 \n")
        draw.write("world ymax 2 \n")
        draw.write("world xmin -1 \n")
        draw.write(f"world xmax {len(residues)} \n")
        draw.write("  yaxis  tick major 0.5 \n")
        draw.write("  yaxis  tick minor ticks 1 \n")

    converter = os.path.join(pref.install_dir, "3_to_1_letter.sh")
    subprocess.run([converter], cwd=directory, check=True)


def process_directory(directory, frames, minimum_fraction):
    """Generate contact summaries and plot data in one directory."""
    required_paths = [os.path.join(directory, name) for name in CONTACT_FILES]
    required_paths.append(os.path.join(directory, "P_RMSF.dat"))
    missing = [os.path.basename(path) for path in required_paths if not os.path.isfile(path)]
    if missing:
        raise FileNotFoundError("missing " + ", ".join(missing))

    contact_counts = [read_contact_counts(path) for path in required_paths[:-1]]
    residue_names = read_residue_names(required_paths[-1])
    minimum_count = minimum_fraction * frames
    rows = []
    for residue_number in range(1, 999):
        counts = [counter[residue_number] for counter in contact_counts]
        if sum(counts) > minimum_count:
            if residue_number not in residue_names:
                raise ValueError(f"residue {residue_number} not found in P_RMSF.dat")
            rows.append((residue_names[residue_number], *counts))

    summary_path = os.path.join(directory, "PL_Contacts_Summary.dat")
    write_summary(summary_path, rows)
    residues = [row[0] for row in rows]
    fractions = [
        np.asarray([row[index] / frames for row in rows], dtype=float)
        for index in range(1, 7)
    ]
    write_grace_files(directory, residues, fractions)
    return residues, fractions


def plot_results(results):
    """Show one stacked contact chart for every processed directory."""
    figure, axes = plt.subplots(
        len(results),
        1,
        figsize=(10, max(5, 4 * len(results))),
        squeeze=False,
    )
    for axis, (label, residues, fractions) in zip(axes[:, 0], results):
        indices = np.arange(len(residues))
        bottom = np.zeros(len(residues))
        for values, contact_label in zip(fractions, CONTACT_LABELS):
            axis.bar(indices, values, bottom=bottom, label=contact_label)
            bottom = np.add(bottom, values)
        axis.set_title(label)
        axis.set_xlabel("Residue #")
        axis.set_ylabel("Interactions fraction")
        axis.set_xticks(indices, residues, rotation=0, ha="center")
        axis.legend()
    figure.tight_layout()
    plt.show()


def main():
    pattern = ask_for_value(
        "Enter directories glob\nExample: data?/raw-data/",
        "Input",
        "data?/raw-data/",
    )
    if pattern is None:
        return

    directories = find_directories(pattern.strip())
    if not directories:
        show_error(f"No directories match:\n{pattern}")
        return

    frames_text = ask_for_value("Enter number of frames", "Input", "1000")
    if frames_text is None:
        return
    minimum_text = ask_for_value(
        "Enter minimum frame fraction to be displayed", "Input", "0.20"
    )
    if minimum_text is None:
        return

    try:
        frames = float(frames_text)
        minimum_fraction = float(minimum_text)
        if frames <= 0:
            raise ValueError("number of frames must be greater than zero")
        if minimum_fraction < 0:
            raise ValueError("minimum frame fraction cannot be negative")
    except ValueError as exc:
        show_error(f"Invalid input: {exc}")
        return

    results = []
    errors = []
    for directory in directories:
        try:
            residues, fractions = process_directory(
                directory, frames, minimum_fraction
            )
            results.append((os.path.relpath(directory), residues, fractions))
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            errors.append(f"{directory}: {exc}")

    with open("temp.log", "w", encoding="utf-8") as log:
        log.write("\n".join(errors))

    if errors:
        show_error("Some directories could not be processed:\n\n" + "\n".join(errors))
    if results:
        plot_results(results)


if __name__ == "__main__":
    main()
