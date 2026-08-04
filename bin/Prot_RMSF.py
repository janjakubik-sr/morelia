#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract and plot protein RMSF data from one or more directories."""

import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import wx


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


def read_rmsf(input_path):
    """Read residue, RMSF, and ligand-contact values from P_RMSF.dat."""
    rmsf_rows = []
    contact_rows = []
    with open(input_path, encoding="utf-8") as input_file:
        for line in input_file:
            fields = line.split()
            if len(fields) < 6:
                continue
            try:
                float(fields[0])
                float(fields[5])
            except ValueError:
                # Ignore the column heading and malformed lines.
                continue
            rmsf_rows.append((fields[0], fields[5]))
            contact_rows.append(
                (fields[0], fields[5] if fields[3] == "Yes" else "0")
            )

    if not rmsf_rows:
        raise ValueError("no numeric RMSF rows found")
    return rmsf_rows, contact_rows


def write_grace_file(draw_path, rmsf_path, contacts_path, x_max, y_max):
    """Write the Grace plotting instructions for one directory."""
    grace_rmsf_path = rmsf_path.replace('"', '\\"')
    grace_contacts_path = contacts_path.replace('"', '\\"')
    with open(draw_path, "w", encoding="utf-8") as draw:
        draw.write(f'read xy "{grace_rmsf_path}" \n')
        draw.write("s0 symbol 0 \n")
        draw.write("s0 line type 1 \n")
        draw.write("s0 line color 1 \n")
        draw.write(f'read xy "{grace_contacts_path}" \n')
        draw.write("s1 symbol 0 \n")
        draw.write("s1 line type 0 \n")
        draw.write("s1 line color 2 \n")
        draw.write("s1 dropline on \n")
        draw.write("page layout portrait \n")
        draw.write("view xmin 0.2 \n")
        draw.write("view xmax 0.8 \n")
        draw.write("view ymin 0.85 \n")
        draw.write("view ymax 1.25 \n")
        draw.write("world xmin 0 \n")
        draw.write(f"world xmax {x_max} \n")
        draw.write("world ymin 0 \n")
        draw.write(f"world ymax {y_max} \n")
        draw.write('xaxis label "Residue #"  \n')
        draw.write('xaxis label font "Helvetica"  \n')
        draw.write("xaxis label char size 1.0 \n")
        draw.write('xaxis ticklabel font "Helvetica"  \n')
        draw.write("xaxis ticklabel char size 0.8 \n")
        draw.write("xaxis tick major size 0.8 \n")
        draw.write("xaxis tick minor size 0.4 \n")
        draw.write('yaxis label "RMSF [A]"  \n')
        draw.write('yaxis label font "Helvetica"  \n')
        draw.write("yaxis label char size 1.0 \n")
        draw.write('yaxis ticklabel font "Helvetica"  \n')
        draw.write("yaxis ticklabel char size 0.8 \n")
        draw.write("yaxis tick major size 0.8 \n")
        draw.write("yaxis tick minor size 0.4 \n")
        draw.write("legend 0.2, 0.7 \n")
        draw.write("legend char size 0.8 \n")
        draw.write('legend font "Helvetica" \n')


def process_directory(directory):
    """Create RMSF data and Grace files in one directory."""
    input_path = os.path.join(directory, "P_RMSF.dat")
    if not os.path.isfile(input_path):
        raise FileNotFoundError("P_RMSF.dat not found")

    rmsf_rows, contact_rows = read_rmsf(input_path)
    rmsf_path = os.path.join(directory, "Prot_RMSF.dat")
    contacts_path = os.path.join(directory, "PL_contacts.dat")
    np.savetxt(rmsf_path, np.asarray(rmsf_rows), fmt="%s")
    np.savetxt(contacts_path, np.asarray(contact_rows), fmt="%s")

    data = np.asarray(rmsf_rows, dtype=float)
    x_max = max(float(np.amax(data[:, 0])) * 1.1, 1.0)
    y_max = max(float(np.amax(data[:, 1])) * 1.1, 1.0)
    draw_path = os.path.join(directory, "Prot_RMSF.draw")
    write_grace_file(draw_path, rmsf_path, contacts_path, x_max, y_max)
    return data


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

    results = []
    errors = []
    for directory in directories:
        try:
            data = process_directory(directory)
            results.append((os.path.relpath(directory), data))
        except (OSError, ValueError) as exc:
            errors.append(f"{directory}: {exc}")

    with open("temp.log", "w", encoding="utf-8") as log:
        log.write("\n".join(errors))

    if errors:
        show_error("Some directories could not be processed:\n\n" + "\n".join(errors))
    if not results:
        return

    plt.title("Protein RMSF")
    for label, data in results:
        plt.plot(data[:, 0], data[:, 1], label=label)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel("Residue #")
    plt.ylabel("RMSF [A]")
    if len(results) > 1:
        plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
