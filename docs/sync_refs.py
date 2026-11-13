"""
sync_refs.py -- keep "Theory eq./Sec./Table" references in the code in step
with the compiled PDF.

Every reference in the code is written as   eq. 3.2 [eq:pconv]
i.e. a number followed by the LaTeX label it refers to. After compiling
tlecuil_theory.tex, this script reads the .aux file and:

    python docs/sync_refs.py --check   # exit 1 if any number is out of date
    python docs/sync_refs.py --fix     # rewrite the numbers in place

This is a traceability control: it guarantees that a reader who follows a
code comment to the PDF lands on the right equation. CI runs --check.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
AUX = os.path.join(HERE, "tlecuil_theory.aux")
FILES = ["tlecuil.py", "benchmarks.py", "test_tlecuil.py", "web/tlecuil.html", "docs/make_figures.py"]
REF = re.compile(r"\b(eqs?\.|Sec\.|Table) ([\dA-Z][\d.]*)(?:-([\d.]+))? \[([a-z]+:[a-z]+)(?:, ([a-z]+:[a-z]+))?\]")


def labels():
    out = {}
    with open(AUX) as fh:
        for m in re.finditer(r"\\newlabel\{([a-z]+:[a-z]+)\}\{\{([^}]*)\}", fh.read()):
            out[m.group(1)] = m.group(2)
    return out


def main(fix: bool) -> int:
    lab = labels()
    bad = 0
    for rel in FILES:
        path = os.path.join(ROOT, rel)
        text = open(path).read()

        def repl(m):
            nonlocal bad
            kw, a, b, l1, l2 = m.groups()
            want_a = lab.get(l1)
            want_b = lab.get(l2) if l2 else None
            if want_a is None or (l2 and want_b is None):
                print(f"{rel}: unknown label in '{m.group(0)}'")
                bad += 1
                return m.group(0)
            new = f"{kw} {want_a}" + (f"-{want_b}" if l2 else "") + f" [{l1}" + (f", {l2}" if l2 else "") + "]"
            if new != m.group(0):
                bad += 1
                print(f"{rel}: '{m.group(0)}' -> '{new}'")
            return new

        new_text = REF.sub(repl, text)
        if fix and new_text != text:
            open(path, "w").write(new_text)
    if bad and not fix:
        print(f"{bad} reference(s) out of date; run with --fix")
        return 1
    print("references", "updated" if fix and bad else "OK")
    return 0


if __name__ == "__main__":
    sys.exit(main(fix="--fix" in sys.argv))
