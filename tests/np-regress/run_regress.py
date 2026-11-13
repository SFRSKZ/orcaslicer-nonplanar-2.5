#!/usr/bin/env python3
"""Nonplanar regression test: slice known models headlessly and check the G-code.

Usage: python run_regress.py <orca-slicer.exe> [model-name ...]
Catches the 2.5-port class of bug: a nonplanar slice that stops with an
internal "unwritten extrusion" error, or silently prints flat.
"""
import json, os, re, subprocess, sys, tempfile, glob

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 0.30  # allowed relative drift on npZ / deferred counts
BAD = re.compile(r"unwritten|nonplanar fallback|internal check|exception|fatal", re.I)

def stats(gcode):
    s = open(gcode, encoding="utf-8", errors="ignore").read()
    npz = 0
    for l in s.split("\n"):
        if l.startswith(("G1", "G0")) and " Z" in l:
            m = re.search(r" E([-\d.]+)", l)
            if m and float(m.group(1)) > 0:
                npz += 1
    return dict(npZ=npz, deferred=s.count(";PATH_DEFERRED from"),
                hdr="OrcaSlicer Nonplanar" in s)

def near(got, want):
    return abs(got - want) <= TOL * want

def main():
    exe = sys.argv[1]
    golden = json.load(open(os.path.join(HERE, "golden.json")))
    names = sys.argv[2:] or list(golden)
    fails = 0
    for n in names:
        g = golden[n]
        root = os.environ.get("NP_OUT")  # keep G-code here (for before/after diffs)
        out = os.path.join(root, n) if root else tempfile.mkdtemp(prefix="npreg-")
        os.makedirs(out, exist_ok=True)
        model = os.path.join(HERE, "models", n + ".3mf")
        p = subprocess.run([exe, "--slice", "0", "--outputdir", out, model],
                           capture_output=True, text=True)
        errs = []
        if p.returncode != 0:
            errs.append("exit %d" % p.returncode)
        if BAD.search(p.stderr or ""):
            errs.append("stderr: " + BAD.search(p.stderr).group(0))
        gc = glob.glob(os.path.join(out, "*.gcode"))
        if not gc:
            errs.append("no gcode")
        else:
            st = stats(gc[0])
            if not st["hdr"]:
                errs.append("missing Nonplanar header")
            if g["np"]:
                if st["npZ"] == 0: errs.append("no nonplanar moves (printed flat)")
                elif not near(st["npZ"], g["npZ"]): errs.append("npZ %d vs %d" % (st["npZ"], g["npZ"]))
                if st["deferred"] == 0: errs.append("PATH_DEFERRED count is 0")
                elif not near(st["deferred"], g["deferred"]): errs.append("deferred %d vs %d" % (st["deferred"], g["deferred"]))
            else:
                if st["npZ"] != 0 or st["deferred"] != 0:
                    errs.append("planar model got nonplanar moves: %s" % st)
        print(("FAIL " if errs else "OK   ") + n, "; ".join(errs))
        fails += bool(errs)
    sys.exit(1 if fails else 0)

main()
