#!/usr/bin/env python3
"""
make_usertopo.py -- generate a YASARA usertopo_<NAME>.fof residue topology for a
silicon-containing ligand, so that YASARA can minimize/simulate it despite the
AM1BCC limitation (no Si bond-charge corrections).

How it works (YASARA manual, Recipes > Run an MD simulation > Preparing the force
field, "atoms unknown to (bio)organic chemistry"):
  1. Load the ligand into YASARA and run Clean (fixes valences/protonation,
     assigns canonical atom names) -> canonical mol2 snapshot.
  2. Replace Si by an sp3 carbon ("C-analog"), let AutoSMILES/GAFF2/AM1BCC
     parameterize the all-organic analog (works, AM1BCC covers C,H,N,O,F...).
  3. Swap Si back into the topology (GAFF2 type 'si', formal charge +0.40 by
     default, the charge delta spread over the directly bound atoms), using the
     canonical soup names so the template matches the molecule by residue name
     and atom names.

Usage:
  make_usertopo.py LIGAND.sdf|LIGAND.mol2 [--name NAME] [--resname RES]
                   [--sicharge 0.40] [--ff AMBER03] [--no-upd] [--test]

At runtime, after Clean, rename the molecule's residue to the reported resname:
      NameRes Mol A,<RES>
before starting the minimization/simulation (see the generated *_min.mcr).

Requires: writable YASARA installation (default /home/roshi/yasara),
          force fields rebuilt with 'yasara -txt -upd' (done automatically
          unless --no-upd).
"""

import argparse
import collections
import datetime
import os
import re
import shutil
import subprocess
import sys
import tempfile

import preferences as pref
YASARA_dir=(pref.YASARA_dir)

YASARA = (str(YASARA_dir)+"yasara")
FOFDIR = (str(YASARA_dir)+"fof")
OBABEL = shutil.which("obabel")


def obabel_smi(path):
    """Canonical SMILES for chemistry-fidelity comparison: hydrogens deleted
    (-d) so explicit-vs-implicit H style cannot create false mismatches, and
    no chiral/isotopic markings (-xi) so merely *annotated* stereochemistry
    (unspecified in the input, assigned from 3D later) does not count as a
    change either; true composition/bond-order changes still register."""
    if not OBABEL:
        return None
    p = subprocess.run([OBABEL, path, "-osmi", "-d", "-xi"], capture_output=True, text=True, timeout=120)
    if p.returncode != 0 or not p.stdout.split():
        return None
    return p.stdout.split()[0].strip()


def norm_smi(s):
    """Achiral normalization for chemistry-fidelity comparison (stereo is kept
    in the 3D coordinates; here we check composition and bond orders)."""
    return re.sub(r"@+", "", s.split()[0])


def build_selfnamed_mol2(infile, outfile):
    """Convert ligand to mol2 via obabel (missing hydrogens are added, bond
    orders kept), rewrite atom names deterministically (C, C1, C2, ..., Si ->
    SI, SI1, ...), residue name UNL. Charge column is zeroed for C/Si rows
    (bogus OpenBabel formal charges), hetero-atom charges (e.g. N+) are kept."""
    tmpc = tempfile.mktemp(suffix=".mol2")
    p = subprocess.run([OBABEL, infile, "-O", tmpc], capture_output=True, text=True, timeout=120)
    if p.returncode != 0 or not os.path.exists(tmpc):
        sys.exit("ERROR: obabel conversion failed for " + infile)
    counters = collections.Counter()
    out, sec = [], None
    for line in open(tmpc):
        if line.startswith("@<TRIPOS>"):
            if line.strip().startswith("@<TRIPOS>UNITY"):
                sec = "unity"  # drop UNITY sections: per-atom charge attrs, handled via charge col
                continue
            sec = "a" if line.startswith("@<TRIPOS>ATOM") else \
                  ("b" if line.startswith("@<TRIPOS>BOND") else None)
            out.append(line.rstrip("\n"))
            continue
        if sec == "unity":
            continue
        t = line.split()
        if sec == "a" and len(t) >= 8:
            el = t[5].split(".")[0]
            sym = "SI" if el.lower() == "si" else el.upper()
            n = counters[el]
            counters[el] += 1
            name = sym if n == 0 else "%s%d" % (sym, n)
            charge = t[8] if el.upper() not in ("C", "SI") else "0.0000"
            line = ("%8d %-8s %9.4f %9.4f %9.4f %-6s %4s %-6s %9.4f"
                    % (int(t[0]), name, float(t[2]), float(t[3]), float(t[4]),
                       t[5], int(t[6]), "UNL", float(charge)))
        out.append(line.rstrip("\n"))
    os.remove(tmpc)
    write(outfile, "\n".join(out) + "\n")
    return outfile


# ---------------------------------------------------------------- helpers
def run_yasara(mcr_path, timeout=580):
    """Run a YASARA macro in plain text mode, return console output."""
    out = subprocess.run(
        [YASARA, "-txt", mcr_path], capture_output=True, text=True, timeout=timeout
    )
    return out.stdout + out.stderr


def write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def parse_mol2(fn):
    """Return atoms {seq: (name, element, resname)}, bonds [(a1,a2,type)] (1-based seq)."""
    atoms, bonds, sec = {}, [], None
    for line in open(fn):
        if line.startswith("@<TRIPOS>ATOM"):
            sec = "a"
            continue
        if line.startswith("@<TRIPOS>BOND"):
            sec = "b"
            continue
        if line.startswith("@<TRIPOS>"):
            sec = None
            continue
        t = line.split()
        if not t:
            continue
        if sec == "a" and len(t) >= 8:
            el = t[5].split(".")[0]
            atoms[int(t[0])] = (t[1], el, t[7])
        elif sec == "b" and len(t) >= 4:
            bonds.append((int(t[1]), int(t[2]), t[3]))
    return atoms, bonds


def mol2_lines(fn):
    return open(fn).read().splitlines()


def free_carbon_names(atoms, count):
    """Find carbon names C90, C91, ... not used in the soup."""
    used = {n for n, _, _ in atoms.values()}
    names, i = [], 90
    while len(names) < count:
        cand = "C%d" % i
        if cand not in used:
            names.append(cand)
        i += 1
    return names


def derive_resname(ligname, fofdir):
    """Derive a unique 3-letter residue name: L + last two digits of trailing
    number (PC_02_19 -> L19), else first 3 alphanumerics."""
    m = re.search(r"(\d+)\D*$", ligname)
    if m:
        res = "L%02d" % (int(m.group(1)) % 100)
        res = res[-3:]
    else:
        res = re.sub(r"[^A-Za-z0-9]", "", ligname.upper())[:3].ljust(3, "X")
    # collision check against all PREP topology residue names in the fof dir
    # (the file we are about to (re)write is excluded, so regeneration is OK)
    taken = set()
    dest = "usertopo_%s.fof" % ligname
    for f in os.listdir(fofdir):
        if not f.endswith(".fof") or f == dest:
            continue
        try:
            for line in open(os.path.join(fofdir, f), errors="ignore"):
                mm = re.match(r"^\s*([A-Za-z0-9]{2,4})\s\sINT\s", line)
                if mm:
                    taken.add(mm.group(1).upper())
        except OSError:
            pass
    if res.upper() in taken:
        # automatic fallback: keep the two digits, cycle the leading letter
        # (first letter of the ligand name first, then A..Z)
        digits = res[-2:]
        first = re.sub(r"[^A-Za-z]", "", ligname.upper())
        letters = ([first[0]] if first else []) + [chr(c) for c in range(ord("A"), ord("Z") + 1)]
        for letter in dict.fromkeys(letters):
            cand = (letter + digits)[:3]
            if cand.upper() not in taken:
                log("  residue name '%s' already taken, using '%s' instead" % (res, cand))
                return cand.upper()
        sys.exit(
            "ERROR: no free residue name available for digits '%s'.\n"
            "       Choose one explicitly with --resname XXX." % res
        )
    return res.upper()


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ligand", help="ligand structure file (.sdf or .mol2, single molecule)")
    ap.add_argument("--name", help="ligand name (default: filename stem)")
    ap.add_argument("--resname", help="3-letter residue name for the template (default: derived)")
    ap.add_argument("--sicharge", type=float, default=0.40,
                    help="partial charge assigned to each Si atom (default: +0.40)")
    ap.add_argument("--ff", default="AMBER03", help="force field used for analog parameterization (default: AMBER03)")
    ap.add_argument("--no-upd", action="store_true", help="do not rebuild force fields afterwards")
    ap.add_argument("--test", action="store_true", help="run a test energy minimization at the end")
    ap.add_argument("--keep-tmp", action="store_true", help="keep temporary directory (debug)")
    args = ap.parse_args()

    ligfile = os.path.abspath(args.ligand)
    if not os.path.isfile(ligfile):
        sys.exit("ERROR: file not found: " + ligfile)
    ext = os.path.splitext(ligfile)[1].lower()
    if ext not in (".sdf", ".mol2"):
        sys.exit("ERROR: only .sdf and .mol2 inputs are supported")
    ligname = args.name or os.path.splitext(os.path.basename(ligfile))[0]
    ligname = re.sub(r"[^A-Za-z0-9_-]", "_", ligname)
    resname = (args.resname or derive_resname(ligname, FOFDIR)).upper()
    if not re.match(r"^[A-Z][A-Z0-9]{0,2}$", resname):
        sys.exit("ERROR: residue name must be 1-3 alphanumerics starting with a letter")

    tmp = tempfile.mkdtemp(prefix="usertopo_")
    log("working in " + tmp)

    smi_in = obabel_smi(ligfile)
    if smi_in:
        log("  input SMILES: " + smi_in)
    else:
        log("  WARNING: obabel not available/failed, chemistry validation disabled")

    # ---- step 1a: load + Clean in YASARA, export canonical mol2 -------------
    canon = os.path.join(tmp, "canon.mol2")
    mcr = ["OnError Exit", "ForceField " + args.ff,
           ("LoadSDF " if ext == ".sdf" else "LoadMol2 ") + ligfile,
           "Clean",
           "nsi=CountAtom Element Si",
           "Print (nsi)",
           "SaveMol2 All," + canon,
           "Exit"]
    mcr_file = os.path.join(tmp, "step1.mcr")
    write(mcr_file, "\n".join(mcr) + "\n")
    log("step 1: load + Clean ...")
    out = run_yasara(mcr_file)
    if not os.path.isfile(canon):
        tail = "\n".join(out.splitlines()[-25:])
        sys.exit("ERROR: YASARA Clean/export failed.\n" + tail)
    m = re.search(r"Macro>Print (\d+)", out)
    nsi = int(m.group(1)) if m else -1
    if nsi < 1:
        sys.exit("ERROR: no silicon atom found in the ligand")
    log("  ligand contains %d silicon atom(s)" % nsi)

    # ---- step 1b: chemistry validation; fall back to no-Clean path ---------
    # Clean reassigns bond orders and can change chemistry (e.g. saturate a
    # genuine C=C and add spurious hydrogens). If the cleaned molecule's SMILES
    # differs from the input file's SMILES, bypass Clean entirely and use a
    # self-named mol2 with the file's own bond orders.
    mode = "clean"
    smi_a = obabel_smi(canon)
    if smi_in and smi_a and norm_smi(smi_a) != norm_smi(smi_in):
        mode = "noclean"
        log("  Clean altered the molecule!")
        log("    cleaned : " + smi_a)
        log("    input   : " + smi_in)
        log("    -> using the no-Clean path (file bond orders kept)")
        canon = build_selfnamed_mol2(ligfile, os.path.join(tmp, "canon_nc.mol2"))
        smi_b = obabel_smi(canon)
        if not smi_b or norm_smi(smi_b) != norm_smi(smi_in):
            sys.exit("ERROR: even the no-Clean path does not reproduce the input chemistry.\n"
                     "  input: %s\n  got  : %s" % (smi_in, smi_b))
        # ink the canonical structure next to the ligand so runtime macros use it
        canon_dest = os.path.join(os.path.dirname(ligfile), "%s_canon.mol2" % ligname)
        write(canon_dest, open(canon).read())
        canon = canon_dest
    log("  mode: %s" % mode)

    atoms, bonds = parse_mol2(canon)
    resname_in = next(iter(atoms.values()))[2]
    si_ids = [i for i, (n, e, r) in atoms.items() if e.lower() == "si"]
    adj = collections.defaultdict(list)
    for a1, a2, bt in bonds:
        adj[a1].append(a2)
        adj[a2].append(a1)
    si_neighbors = {i: [j for j in adj[i] if atoms[j][1] != "H"] for i in si_ids}

    # ---- step 2: build carbon analog & parameterize with AutoSMILES ---------
    newnames = free_carbon_names(atoms, len(si_ids))
    standin = {si: nn for si, nn in zip(si_ids, newnames)}  # si seq -> analog C name
    analog = os.path.join(tmp, "analog.mol2")
    replaced = 0
    outlines = []
    in_atoms = False
    for line in mol2_lines(canon):
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
        elif line.startswith("@<TRIPOS>BOND"):
            in_atoms = False
        t = line.split()
        if in_atoms and len(t) >= 8 and t[5].split(".")[0].lower() == "si":
            si = int(t[0])
            line = re.sub(r"(\s)%s(\s)" % re.escape(t[1]), r"\g<1>%s\g<2>" % standin[si], line, count=1)
            line = re.sub(r"(\s)%s(\s)" % re.escape(t[5]), r"\g<1>C.3\g<2>", line, count=1)
            replaced += 1
        if resname_in:
            line = line.replace(" %s " % resname_in, " ANA ")
        outlines.append(line)
    assert replaced == len(si_ids), "Si replacement failed"
    write(analog, "\n".join(outlines) + "\n")

    # AutoSMILES caches successful assignments in fof/<ff>.cache; a cache hit
    # would suppress the topology report we need to parse, so hide the cache
    # file of the selected force field for this run (restored afterwards).
    cache = os.path.join(FOFDIR, args.ff.lower() + ".cache")
    cache_bak = None
    if os.path.exists(cache):
        cache_bak = cache + ".hide"
        os.replace(cache, cache_bak)
        log("  temporarily hid AutoSMILES cache " + cache)
    mcr = ["OnError Exit", "ForceField " + args.ff,
           "LoadMol2 " + analog,
           "Sim Init",
           "Exit"]
    mcr_file = os.path.join(tmp, "step2.mcr")
    write(mcr_file, "\n".join(mcr) + "\n")
    log("step 2: parameterizing carbon analog with AutoSMILES/AM1BCC ...")
    try:
        out = run_yasara(mcr_file, timeout=580)
    finally:
        if cache_bak:
            os.replace(cache_bak, cache)
            log("  restored AutoSMILES cache")
    if "Topology for residue ANA" not in out:
        tail = "\n".join(out.splitlines()[-30:])
        sys.exit("ERROR: analog parameterization failed.\n" + tail)
    start = out.index("Topology for residue ANA")
    end = out.index("DONE", start) + 4
    block = out[start:end].splitlines()

    # ---- step 3: parse analog topology, restore Si --------------------------
    atom_lines, sec_loop, sec_impr, section = [], [], [], None
    for l in block:
        t = l.split()
        if t and t[0].isdigit() and len(t) >= 11 and t[1] != "DUMM":
            atom_lines.append(l)
        elif l.startswith("LOOP"):
            section = "loop"
            continue
        elif l.startswith("IMPROPER"):
            section = "impr"
            continue
        elif l.startswith("DONE"):
            section = None
        elif section and t:
            (sec_loop if section == "loop" else sec_impr).append(t)
    # exclude the 3 DUMM rows when indexing by sequential number
    recs = []
    for l in atom_lines:
        recs.append({"line": l, "seq": int(l[0:4]), "name": l[6:12].strip(),
                     "type": l[12:14].strip(), "na": int(l.split()[4]),
                     "q": float(l[-10:])})
    by_seq = {r["seq"]: r for r in recs}
    standin_names = set(standin.values())
    n_standin = [r for r in recs if r["name"] in standin_names]
    if len(n_standin) != len(si_ids):
        sys.exit("ERROR: stand-in carbon(s) not found in the analog topology")

    inv_standin = {v: k for k, v in standin.items()}  # analog C name -> si seq
    # charge surgery: per Si, spread (q_si - q_standin) over its bound neighbors
    newq = {r["name"]: r["q"] for r in recs}
    for r in n_standin:
        si = inv_standin[r["name"]]
        delta = args.sicharge - r["q"]
        # neighbor prep names of the stand-in row, via canonical mol2 bond list
        nb_prepnames = []
        for j in si_neighbors[si]:
            soup_name = atoms[j][0]
            # soup name == prep name except for hydrogens; neighbors are heavy
            if soup_name in newq:
                nb_prepnames.append(soup_name)
        if not nb_prepnames:
            sys.exit("ERROR: could not map Si neighbors into the topology")
        share = round(delta / len(nb_prepnames), 5)
        for n in nb_prepnames:
            newq[n] = round(newq[n] - share, 5)
        newq[r["name"]] = args.sicharge
        # fix rounding drift on the last neighbor so the sum is exact
        drift = round(delta - share * len(nb_prepnames), 5)
        newq[nb_prepnames[-1]] = round(newq[nb_prepnames[-1]] - drift, 5)
    # force net charge to the integer nearest the analog sum
    total = round(sum(newq.values()), 5)
    target = round(total)
    if abs(total - target) > 1e-6:
        big = max((n for n in newq if n not in standin_names), key=lambda n: abs(newq[n]))
        newq[big] = round(newq[big] - (total - target), 5)
    log("step 3: net charge of the residue = %.5f" % round(sum(newq.values()), 5))

    def soup_h_name(r):
        """YASARA soup hydrogen convention: 'H' + trailing digits of the heavy
        parent's name, duplicated for equivalent H (e.g. H7, H7 on C7).
        Parent without digits (bare 'C' or 'SI') -> plain 'H' for all its H."""
        par = by_seq.get(r["na"])
        if par:
            m = re.search(r"(\d+)$", par["name"])
            if m:
                return "H" + m.group(1)
            if par["name"][0].upper() in ("C", "S"):
                return "H"
        log("  WARNING: keeping topology name for hydrogen '%s' (parent '%s'), verify manually" % (r["name"], par["name"] if par else "?"))
        return r["name"]  # fallback (e.g. H on N/O without digits): keep name

    def final_name(r):
        if r["name"] in standin_names:
            return atoms[inv_standin[r["name"]]][0]  # canonical Si name, e.g. SI
        if mode == "noclean":
            # names come from our own mol2 and are what the runtime macro loads:
            # keep them exactly as the analog report printed them
            return r["name"]
        if r["type"].lower().startswith("h") or r["name"].startswith("H") or re.match(r"^\d+H$", r["name"]):
            return soup_h_name(r)
        return r["name"]

    # ---- step 4: write the usertopo file ------------------------------------
    today = datetime.date.today().isoformat()
    out_l = []
    out_l.append("Topology for residue %s, silicon ligand %s: AutoSMILES/GAFF2/AM1BCC parameters of the C-analog, Si restored (charge %+.2f, delta spread over bound atoms). Net charge %d. Created %s by make_usertopo.py"
                 % (resname, ligname, args.sicharge, int(target), today))
    out_l.append(" Created with YASARA AutoSMILES: www.yasara.org/autosmiles")
    out_l.append(" %s  INT    1" % resname)
    out_l.append(" CORR OMIT DU   BEG")
    out_l.append("   0.00000")
    for l in block[5:8]:
        out_l.append(l)
    for r in recs:
        typ = "si" if r["name"] in standin_names else r["type"]
        name = final_name(r)
        l = r["line"]
        out_l.append(l[0:6] + ("%-6s" % name) + ("%-2s" % typ) + l[14:-10] + ("%10.5f" % newq[r["name"]]))
    out_l.append("")
    out_l.append("LOOP CLOSING EXPLICIT")
    ren = {r["name"]: final_name(r) for r in recs}
    for t in sec_loop:
        out_l.append(" " + "  ".join(ren.get(p, p) for p in t) + " ")
    out_l.append("")
    out_l.append("IMPROPER")
    for t in sec_impr:
        out_l.append(" " + "  ".join(ren.get(p, p) for p in t) + " ")
    out_l.append("")
    out_l.append("DONE")
    dest = os.path.join(FOFDIR, "usertopo_%s.fof" % ligname)
    if os.path.exists(dest):
        os.replace(dest, dest + ".bak")
        log("  existing %s moved to %s.bak" % (dest, dest))
    write(dest, "\n".join(out_l) + "\n")
    log("step 4: wrote " + dest)

    # ---- step 5: companion macro for minimization ---------------------------
    mcr_dest = os.path.join(os.path.dirname(ligfile), "%s_min.mcr" % ligname)
    if mode == "noclean":
        load_lines = ["LoadMol2 " + canon,   # self-named, file bond orders, no Clean
                      "NameRes Mol A," + resname]
    else:
        load_lines = [("LoadSDF " if ext == ".sdf" else "LoadMol2 ") + ligfile,
                      "Clean",
                      "NameRes Mol A," + resname]
    mcr = ["# Energy minimization of %s (silicon ligand, residue %s, %s path)" % (ligname, resname, mode),
           "OnError Exit",
           "ForceField " + args.ff] + load_lines + [
           "Experiment Minimization",
           "Experiment On",
           "Wait ExpEnd",
           "SaveYOB All,%s_minimized.yob" % ligname,
           "SaveMol2 All,%s_minimized.mol2" % ligname,
           "Exit"]
    write(mcr_dest, "\n".join(mcr) + "\n")
    log("        runtime macro: " + mcr_dest)

    # ---- step 6: rebuild force fields ----------------------------------------
    if not args.no_upd:
        log("step 6: rebuilding force fields (yasara -txt -upd) ...")
        out = run_yasara_upd()
        if "ERROR" in out.upper() and "0 errors" not in out.lower():
            tail = "\n".join(out.splitlines()[-30:])
            sys.exit("ERROR: force field rebuild failed.\n" + tail)
        log("  force fields rebuilt")
    else:
        log("step 6: skipped (-upd); run 'yasara -txt -upd' before use")

    log("")
    log("DONE. Residue name to use: %s" % resname)
    log("Runtime recipe: LoadSDF/LoadMol2 -> Clean -> NameRes Mol A,%s -> Simulate." % resname)

    if args.keep_tmp:
        log("temporary files kept in " + tmp)

    if args.test:
        if args.no_upd:
            log("WARNING: --test ignored because force fields were not rebuilt (--no-upd)")
        else:
            log("step 7: test energy minimization ...")
            out = run_yasara(mcr_dest, timeout=580)
            ok = "Energy minimization experiment completed" in out
            for line in out.splitlines():
                if "Energy minimization experiment completed" in line or "ERROR" in line:
                    log("  " + line.strip())
            log("  TEST %s" % ("PASSED" if ok else "FAILED"))
    return dest, resname, mcr_dest


def run_yasara_upd(timeout=580):
    out = subprocess.run([YASARA, "-txt", "-upd"], capture_output=True, text=True, timeout=timeout)
    return out.stdout + out.stderr


def log(msg):
    print(msg, flush=True)


if __name__ == "__main__":
    main()
