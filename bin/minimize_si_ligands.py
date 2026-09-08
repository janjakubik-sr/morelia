#!/usr/bin/env python3
"""
minimize_si_ligands.py -- batch pipeline for silicon-containing ligands.

For every .sdf file in a directory:
  1. scan for silicon atoms (V2000 and V3000 SDF dialects),
  2. if Si found and no force-field template exists yet, generate
     yasara/fof/usertopo_<NAME>.fof via make_usertopo.py,
  3. rebuild the force fields once at the end (if any new template was made),
  4. run the YASARA energy minimization (via the generated <NAME>_min.mcr),
     saving <NAME>_minimized.mol2 and <NAME>_minimized.yob next to the SDF.

Usage:
  minimize_si_ligands.py [directory] [--force] [--dry-run]
                         [--sicharge 0.40] [--ff AMBER03]

Files whose templates already exist are only re-minimized (fast);
use --force to regenerate templates from scratch.
"""

import argparse
import os
import re
import subprocess
import sys

import preferences as pref
install_dir=(pref.install_dir)
YASARA_dir=(pref.YASARA_dir)

HERE = os.path.dirname(os.path.abspath(__file__))
MAKE_TOPO = (str(install_dir)+"make_usertopo.py")
YASARA = (str(YASARA_dir)+"yasara")
FOFDIR = (str(YASARA_dir)+"fof")


def sdf_has_silicon(path):
    """Detect Si atoms in V2000 or V3000 SDF files."""
    try:
        txt = open(path, errors="ignore").read()
    except OSError:
        return False
    if "V3000" in txt:
        return re.search(r"M\s+V30\s+\d+\s+Si\s", txt) is not None
    lines = txt.splitlines()
    for i, l in enumerate(lines):
        if "V2000" in l:
            try:
                na = int(l.split()[0])
            except (ValueError, IndexError):
                return False
            for j in range(i + 1, min(i + 1 + na, len(lines))):
                parts = lines[j].split()
                if len(parts) >= 4 and parts[3] == "Si":
                    return True
            return False
    return False


def run(cmd, timeout=580, cwd=None):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=cwd)
    return p.stdout + p.stderr


def log(msg):
    print(msg, flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("directory", nargs="?", default=".", help="directory with .sdf files (default: current)")
    ap.add_argument("--force", action="store_true", help="regenerate templates even if they exist")
    ap.add_argument("--dry-run", action="store_true", help="only list silicon-containing files")
    ap.add_argument("--sicharge", type=float, default=0.40, help="partial charge for Si atoms (default: +0.40)")
    ap.add_argument("--ff", default="AMBER03", help="force field (default: AMBER03)")
    args = ap.parse_args()

    sdfdir = os.path.abspath(args.directory)
    sdfs = sorted(f for f in os.listdir(sdfdir) if f.lower().endswith(".sdf"))
    if not sdfs:
        sys.exit("No .sdf files found in " + sdfdir)

    si_files, skipped = [], []
    for f in sdfs:
        full = os.path.join(sdfdir, f)
        (si_files if sdf_has_silicon(full) else skipped).append(full)
    log("Found %d .sdf files: %d with silicon, %d without" % (len(sdfs), len(si_files), len(skipped)))
    for f in skipped:
        log("  no Si       : %s" % os.path.basename(f))
    if args.dry_run:
        for f in si_files:
            log("  Si detected : %s" % os.path.basename(f))
        return

    # ---- pass 1: templates --------------------------------------------------
    results = {}
    need_upd = False
    for full in si_files:
        name = re.sub(r"[^A-Za-z0-9_-]", "_", os.path.splitext(os.path.basename(full))[0])
        topo = os.path.join(FOFDIR, "usertopo_%s.fof" % name)
        mcr = os.path.join(sdfdir, "%s_min.mcr" % name)
        results[name] = {"dir": os.path.dirname(full), "mcr": mcr}
        stale = os.path.exists(topo) and os.path.getmtime(topo) < os.path.getmtime(full)
        if stale:
            log("\n>>> SDF is newer than the template for %s - regenerating" % name)
        if args.force or stale or not os.path.exists(topo) or not os.path.exists(mcr):
            log("\n>>> generating template for %s ..." % name)
            cmd = [sys.executable, MAKE_TOPO, full, "--no-upd",
                   "--sicharge", str(args.sicharge), "--ff", args.ff]
            out = run(cmd, timeout=880)
            results[name]["topo"] = "OK" if os.path.exists(topo) else "FAILED"
            if results[name]["topo"] == "OK":
                need_upd = True
            else:
                results[name]["error"] = "\n".join(out.splitlines()[-15:])
                log("  FAILED (see summary)")
        else:
            results[name]["topo"] = "exists"
            log("\n>>> template for %s already exists (skipped, use --force to regenerate)" % name)

    # ---- force field rebuild (once per batch) -------------------------------
    if need_upd:
        log("\n>>> rebuilding force fields (yasara -txt -upd) ...")
        out = run([YASARA, "-txt", "-upd"], timeout=880)
        bad = [n for n, r in results.items() if r.get("topo") == "OK"]
        if "error" in out.lower() and "0 errors" not in out.lower():
            for n in bad:
                if "minimized" not in results[n]:
                    results[n]["topo"] = "UPD-FAILED"
            log("  force field rebuild reported problems - check manually")
        else:
            log("  force fields rebuilt")

    # ---- pass 2: minimization with self-healing retry -----------------------
    def smiles_of(path):
        # -d: hydrogens deleted, -xi: no chiral markings, so that H-style or
        # stereo annotation differences cannot create false mismatches; true
        # composition/bond-order changes still register
        p = subprocess.run(["obabel", path, "-osmi", "-d", "-xi"], capture_output=True, text=True, timeout=120)
        return p.stdout.split()[0].strip() if p.returncode == 0 and p.stdout.split() else None

    def minimize(name):
        return run([YASARA, "-txt", results[name]["mcr"]], timeout=880, cwd=results[name]["dir"])

    def evaluate(name, out):
        """Return (final_energy_str | None, problem | None)."""
        m = re.search(r"Energy minimization experiment completed after (\d+) steps\. Final energy is ([-\d.]+) kJ/mol", out)
        if not m:
            mm = re.search(r"ERROR[^\n]*", out)
            return None, (mm.group(0) if mm else "minimization aborted")
        # chemistry fidelity: minimized structure must match the input file
        smi_in = smiles_of(os.path.join(results[name]["dir"], name + ".sdf"))
        smi_out = smiles_of(os.path.join(results[name]["dir"], name + "_minimized.mol2"))
        if smi_in and smi_out and re.sub(r"@+", "", smi_in) != re.sub(r"@+", "", smi_out):
            return None, "chemistry mismatch: input %s vs minimized %s" % (smi_in, smi_out)
        return "%s steps, %s kJ/mol" % (m.group(1), m.group(2)), None

    need_upd2 = False
    for name, r in results.items():
        if r.get("topo") not in ("OK", "exists"):
            continue
        for attempt in (1, 2):
            log("\n>>> minimizing %s%s ..." % (name, " (retry)" if attempt == 2 else ""))
            out = minimize(name)
            open(os.path.join(r["dir"], "%s_min.log" % name), "w").write(out)
            energy, problem = evaluate(name, out)
            if problem is None:
                r["min"] = "OK (%s)" % energy
                break
            regenerate = ("predefined topology has been found, but" in out
                          or "chemistry mismatch" in problem)
            if not regenerate or attempt == 2:
                r["min"] = "FAILED"
                r["error"] = problem
                break
            log("  %s -> regenerating template" % problem.splitlines()[0][:90])
            out2 = run([sys.executable, MAKE_TOPO, os.path.join(r["dir"], name + ".sdf"),
                        "--no-upd", "--sicharge", str(args.sicharge), "--ff", args.ff], timeout=880)
            if not os.path.exists(os.path.join(FOFDIR, "usertopo_%s.fof" % name)):
                r["min"] = "FAILED"
                r["error"] = "\n".join(out2.splitlines()[-8:])
                break
            r["topo"] = "regenerated"
            need_upd2 = True
            if need_upd2 and attempt == 1:
                log("\n>>> rebuilding force fields (yasara -txt -upd) ...")
                run([YASARA, "-txt", "-upd"], timeout=880)
                log("  force fields rebuilt")
        for ext in ("yob", "mol2"):
            f = os.path.join(r["dir"], "%s_minimized.%s" % (name, ext))
            r[ext] = "written" if os.path.exists(f) else "MISSING"
        if "min" in r:
            log("  " + r["min"])

    # ---- summary -------------------------------------------------------------
    log("\n" + "=" * 72)
    log("SUMMARY")
    log("%-14s %-10s %-28s %-9s %-9s" % ("ligand", "template", "minimization", "mol2", "yob"))
    for name in sorted(results):
        r = results[name]
        log("%-14s %-10s %-28s %-9s %-9s" % (name, r.get("topo", "-"), r.get("min", "-"),
                                            r.get("mol2", "-"), r.get("yob", "-")))
        if "error" in r:
            log("  error: " + r["error"].splitlines()[0])
    n_ok = sum(1 for r in results.values() if str(r.get("min", "")).startswith("OK"))
    log("\n%d/%d silicon ligands minimized successfully" % (n_ok, len(results)))
    sys.exit(0 if n_ok == len(results) else 1)


if __name__ == "__main__":
    main()
