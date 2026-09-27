"""
Random Molecule Generator for MoleditPy (Plugin API v4.0)
=========================================================

Generates random organic molecules and adds them to the 2D canvas.

Three generation strategies:
  * "Drug-like" : a common drug core (ring seed) + 0-4 curated
                  functional-group fragments attached at random positions.
  * "Cyclic"    : a randomly generated ring (3-8 membered, 0-3
                  heteroatoms, optionally aromatic) + random fragments.
  * "Chain"     : Monte-Carlo random walk over the SMILES grammar
                  (branched chains, accidental rings, terminal
                  functional groups and halogens).

Every produced structure is validated by RDKit; invalid attempts are
retried transparently (typically < 200 ms for a whole batch).

Install
-------
Copy this file to your plugin directory and restart MoleditPy, or use
    Plugins > Reload All Plugins
  * Windows : C:\\Users\\<You>\\.moleditpy\\plugins\\
  * Linux   : ~/.moleditpy/plugins/
  * macOS   : ~/.moleditpy/plugins/

Menu entries
------------
  * Tools > Random Molecule > Generator...     (options dialog)
  * Tools > Random Molecule > Quick: Drug-like
  * Tools > Random Molecule > Quick: Cyclic
  * Tools > Random Molecule > Quick: Chain
  * Plugin Toolbar > "Random Mol" button
"""

import logging
import random

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
except Exception:  # pragma: no cover - host always ships RDKit
    Chem = None

# Qt bindings: MoleditPy ships PyQt6 (see the v4.0 plugin manual).
try:  # pragma: no cover - depends on host environment
    from PyQt6.QtWidgets import (
        QApplication, QComboBox, QDialog, QFormLayout, QHBoxLayout,
        QLabel, QPushButton, QPlainTextEdit, QSpinBox, QVBoxLayout,
    )
    _HAS_QT = True
except Exception as exc:  # pragma: no cover
    _HAS_QT = False
    logger.debug("PyQt6 import failed: %s", exc, exc_info=True)

# ---------------------------------------------------------------------------
# Metadata (read by the MoleditPy plugin manager)
# ---------------------------------------------------------------------------
PLUGIN_NAME = "Random Molecule Generator"
PLUGIN_VERSION = "1.0.0"
PLUGIN_AUTHOR = "Author"
PLUGIN_DESCRIPTION = (
    "Generate random organic molecules (drug-like, cyclic or random-walk "
    "chains) and add them to the 2D canvas. All structures are validated "
    "by RDKit."
)
PLUGIN_CATEGORY = "Tools"
PLUGIN_TAGS = ["Generator", "Utility", "Random", "RDKit"]
PLUGIN_DEPENDENCIES = ["rdkit"]
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0"

MODES = ("Drug-like", "Cyclic", "Chain")

# ---------------------------------------------------------------------------
# Core chemistry data
# ---------------------------------------------------------------------------

# Drug cores: SMILES of common rings (digits are part of the structure).
_SEED_RINGS = [
    "C1CCCCC1",               # cyclopentane
    "C1CCCCC1",
    "O1CCCCC1",               # tetrahydrofuran
    "N1CCCCC1",               # pyrrolidine
    "C1NCCCN1",               # piperidine
    "C1NCCNCC1",              # piperazine
    "C1CCOCC1",               # tetrahydropyran
    "O1CCOCC1",               # 1,4-dioxane
    "C1CCC(=O)CC1",           # cyclopentanone
    "c1ccccc1",               # benzene
    "c1ccncc1",               # pyridine
    "c1ccoc1",                # furan
    "c1ccsc1",                # thiophene
    "c1ccc2[nH]ccc2c1",       # indole
    "c1ccc2ccccc2c1",         # naphthalene
    "c1ccc2ncccc2c1",         # quinoline
]

# (fragment SMILES, weight). Weights bias the sampler toward common,
# small groups. Fragments are attached at the molecular-graph level
# (RDKit RWMol), so their internal ring numbering never clashes with
# the core.
_FRAGMENTS = [
    ("C", 6), ("CC", 5), ("CCC", 3), ("C(C)C", 3), ("CC(C)C", 1),
    ("O", 2), ("N", 4), ("F", 2), ("Cl", 2), ("Br", 1),
    ("OC", 4), ("OCC", 3), ("OCCC", 2), ("OC(=O)C", 3),
    ("OC(=O)N", 2), ("OC(=O)CC", 2),
    ("NC", 4), ("N(C)C", 4), ("CCN", 3), ("NC(=O)C", 4),
    ("NC(=O)N(C)C", 2),
    ("C(=O)C", 4), ("C(=O)CC", 2), ("C(=O)O", 3),
    ("C(=O)N", 3), ("C(=O)N(C)C", 2),
    ("SC", 2), ("SCC", 1), ("S(=O)(=O)N", 1),
    ("CC=O", 2), ("CC#N", 1), ("CC(=O)OC", 1),
    ("C2=CC=CC=C2", 3),   # phenyl
    ("c3ccccc3", 2),      # phenyl (aromatic form)
    ("c3ccncc3", 2),      # pyridyl
    ("C3CCCC3", 2),       # cyclopentyl
    ("C3CCNCC3", 2),      # piperidyl
]

# Random-walk (Chain mode) atom pools. (symbol, weight)
_WALK_MID = [("C", 55), ("N", 12), ("O", 10), ("S", 5)]
_WALK_LAST = [("C", 28), ("N", 12), ("O", 15), ("S", 5),
              ("F", 10), ("Cl", 16), ("Br", 6)]
_BRANCH_HALOGENS = [("F", 4), ("Cl", 4), ("Br", 2), ("C", 2), ("N", 1)]

# ---------------------------------------------------------------------------
# Pure generation core (RDKit only -- no Qt needed, so it can be unit-tested)
# ---------------------------------------------------------------------------

def _weighted_choice(pairs, rng):
    """Pick a token from [(token, weight), ...] using the given rng."""
    total = sum(w for _, w in pairs)
    r = rng.random() * total
    acc = 0.0
    for token, w in pairs:
        acc += w
        if r <= acc:
            return token
    return pairs[-1][0]


def _attach_fragment(mol_smiles, fragment, rng):
    """Attach `fragment` (SMILES) to a random atom of `mol_smiles`.

    Works at the molecular-graph level (RWMol), which avoids all SMILES
    token splicing issues. Returns the new structure as a canonical
    SMILES string, or None if the attachment is chemically invalid.
    """
    if Chem is None:
        return None
    core = Chem.MolFromSmiles(mol_smiles)
    if core is None:
        return None
    frag = Chem.MolFromSmiles(fragment)
    if frag is None:
        return None

    # Prefer neutral heavy atoms; fall back to any heavy atom.
    heavy = [i for i in range(core.GetNumAtoms())
             if core.GetAtomWithIdx(i).GetAtomicNum() > 1]
    if not heavy:
        return None
    neutral = [i for i in heavy if core.GetAtomWithIdx(i).GetFormalCharge() == 0]
    attach_atom = rng.choice(neutral if neutral else heavy)

    try:
        rw = Chem.RWMol(core)
        idx_map = {}
        root_idx = None
        for i, f_atom in enumerate(frag.GetAtoms()):
            new_idx = rw.AddAtom(f_atom)
            idx_map[i] = new_idx
            if i == 0:
                root_idx = new_idx
        for f_atom in frag.GetAtoms():
            for bond in f_atom.GetBonds():
                a0, a1 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                other = a1 if a0 == f_atom.GetIdx() else a0
                if other > f_atom.GetIdx():  # each internal bond only once
                    rw.AddBond(idx_map[a0], idx_map[a1], bond.GetBondType())
        rw.AddBond(attach_atom, root_idx, Chem.BondType.SINGLE)
        out = rw.GetMol()
        Chem.SanitizeMol(out)
        return Chem.MolToSmiles(out)
    except Exception:
        return None


def _grow_with_fragments(seed, n_target, rng, max_attach=6):
    """Attach random fragments to `seed` until ~n_target heavy atoms."""
    s = seed
    attached = 0
    m = Chem.MolFromSmiles(s) if Chem else None
    if m is None:
        return None
    while (m.GetNumHeavyAtoms() < n_target and attached < max_attach):
        frag = _weighted_choice(_FRAGMENTS, rng)
        s2 = _attach_fragment(s, frag, rng)
        if s2 is None:
            break  # no clean attachment left; keep what we have
        s = s2
        attached += 1
        m = Chem.MolFromSmiles(s)
        if m is None:
            break
    if m is not None:
        return Chem.MolToSmiles(m)
    return None


def _random_ring(rng):
    """Build a random SMILES ring (optionally aromatic, 0-3 heteroatoms)."""
    size = rng.choice((5, 6, 6, 5, 4, 7))
    aromatic = size in (5, 6) and rng.random() < 0.4

    nhet = 0
    if rng.random() < 0.55:
        nhet = min([rng.choice((1, 1, 2)), size - 2, 3])
        if size == 3:
            nhet = min(nhet, 1)

    atoms = ["C"] * size
    for i in rng.sample(range(size), nhet):
        atoms[i] = rng.choice(["N", "O", "S"])
    if aromatic:
        atoms = [a.lower() if a in "COS" else a for a in atoms]

    digit = rng.choice("123456789")
    toks = []
    for i, a in enumerate(atoms):
        toks.append(a + digit if i in (0, size - 1) else a)
    return "".join(toks)


def _do_walk(n, rng):
    """One attempt at a random-but-grammatically-valid SMILES of n atoms.

    Grammar constraints kept simple so RDKit can reject the rest:
      * a ring digit is opened at most once and closed at most once,
        and only after >= 2 intervening atoms (ring size >= 3),
      * branches "(" ")" are opened <= 3 deep and always closed.
    """
    toks = []
    depth = 0
    ring = None          # ring digit opened but not closed
    ring_pos = -1        # atom count where the ring was opened
    count = 0

    while count < n:
        remaining_after = n - count - 1
        is_last = remaining_after == 0

        # Where does this atom sit?
        if count == 0:
            choice = "cont"
        elif is_last and depth > 0:
            choice = "close"          # never leave a branch dangling
        elif depth > 0 and rng.random() < 0.4:
            choice = "close"
        elif (depth < 3 and remaining_after >= 1 and not is_last
              and rng.random() < 0.35):
            choice = "branch"
        else:
            choice = "cont"

        if choice == "branch":
            # Sometimes the whole branch is a single halogen/methyl.
            if rng.random() < 0.4:
                hal = _weighted_choice(_BRANCH_HALOGENS, rng)
                toks.append("(" + hal + ")")
                count += 1
                if count == n:
                    break
                continue
            toks.append("(")
            depth += 1

        sym = _weighted_choice(_WALK_LAST if is_last else _WALK_MID, rng)

        # Terminal C=O looks like a carbonyl: worth the extra double bond.
        bond = ""
        if is_last and sym == "O":
            bond = "=" if rng.random() < 0.5 else ""
        if bond:
            toks.append(bond)

        # Ring operations (attached to this atom).
        dec = ""
        if (ring is None and count >= 2 and remaining_after >= 2
                and rng.random() < 0.2):
            ring = rng.choice("23456789")
            ring_pos = count
            dec = ring
        elif (ring is not None and (is_last or rng.random() < 0.4)
                and (count - ring_pos) >= 2):
            dec = ring
            ring = None

        toks.append(sym + dec)
        count += 1

        if choice == "close":
            toks.append(")")
            depth -= 1
        if count == n:
            break

    return "".join(toks)


def _one_druglike(nmin, nmax, rng):
    n_target = rng.randint(nmin, nmax)
    for _ in range(50):
        seed = rng.choice(_SEED_RINGS)
        s = _grow_with_fragments(seed, n_target, rng, max_attach=6)
        if s is None:
            continue
        m = Chem.MolFromSmiles(s)
        if m is not None and nmin <= m.GetNumHeavyAtoms() <= n_target + 3:
            return s
    return None


def _one_cyclic(nmin, nmax, rng):
    n_target = rng.randint(nmin, nmax)
    for _ in range(50):
        seed = _random_ring(rng)
        s = _grow_with_fragments(seed, n_target, rng, max_attach=8)
        if s is None:
            continue
        m = Chem.MolFromSmiles(s)
        if m is not None and nmin <= m.GetNumHeavyAtoms() <= n_target + 3:
            return s
    return None


def _one_chain(nmin, nmax, rng):
    n = rng.randint(nmin, nmax)
    for _ in range(300):
        smi = _do_walk(n, rng)
        try:
            m = Chem.MolFromSmiles(smi)
        except Exception:
            m = None
        if m is not None and m.GetNumHeavyAtoms() >= max(2, nmin):
            return Chem.MolToSmiles(m)
    return None


def _mode_key(mode):
    key = str(mode).strip().lower().replace("-", "").replace(" ", "")
    aliases = {
        "drug": "Drug-like", "druglike": "Drug-like", "medicinal": "Drug-like",
        "ring": "Cyclic", "cyclic": "Cyclic", "heterocycle": "Cyclic",
        "chain": "Chain", "walk": "Chain", "aliphatic": "Chain",
        "randomwalk": "Chain", "random": "Drug-like",
    }
    return aliases.get(key, "Drug-like")


def generate_random_molecule(mode="Drug-like", nmin=6, nmax=16, rng=None):
    """Generate one valid random molecule (SMILES string) or None."""
    if Chem is None:
        raise RuntimeError("RDKit is not installed; cannot generate molecules.")
    # Intentionally non-cryptographic: sampling molecules, not secrets.
    rng = rng or random.Random()  # nosec B311
    mode = _mode_key(mode)
    nmin = max(2, int(nmin))
    nmax = max(nmin, int(nmax))

    if mode == "Drug-like":
        return _one_druglike(nmin, nmax, rng)
    if mode == "Cyclic":
        return _one_cyclic(nmin, nmax, rng)
    return _one_chain(nmin, nmax, rng)


def generate_batch(mode="Drug-like", nmin=6, nmax=16, count=3, rng=None):
    """Generate `count` distinct valid molecules (list of SMILES)."""
    mode = _mode_key(mode)
    # Intentionally non-cryptographic: sampling molecules, not secrets.
    rng = rng or random.Random()  # nosec B311
    out, seen = [], set()
    for _ in range(max(1, int(count))):
        for _ in range(20):
            s = generate_random_molecule(mode, nmin, nmax, rng)
            if s is not None and s not in seen:
                seen.add(s)
                out.append(s)
                break
    return out


# ---------------------------------------------------------------------------
# Qt UI (singleton options dialog) -- imported lazily enough to test headless
# ---------------------------------------------------------------------------

if _HAS_QT:

    class GeneratorDialog(QDialog):
        """Options + batch generation dialog (kept alive via register_window)."""

        def __init__(self, context, parent=None):
            super().__init__(parent)
            self.context = context
            self.setWindowTitle("Random Molecule Generator")
            self.resize(480, 380)

            cfg = {}
            try:
                cfg = (context.get_setting("last_config") or {})
            except Exception:
                cfg = {}

            form = QFormLayout()
            self.mode_combo = QComboBox()
            self.mode_combo.addItems(MODES)
            self.mode_combo.setCurrentText(str(cfg.get("mode", MODES[0])))

            self.min_spin = QSpinBox()
            self.min_spin.setRange(2, 30)
            self.min_spin.setValue(int(cfg.get("nmin", 8)))

            self.max_spin = QSpinBox()
            self.max_spin.setRange(3, 80)
            self.max_spin.setValue(int(cfg.get("nmax", 18)))

            self.count_spin = QSpinBox()
            self.count_spin.setRange(1, 25)
            self.count_spin.setValue(int(cfg.get("count", 3)))

            form.addRow("Mode", self.mode_combo)
            form.addRow("Min heavy atoms", self.min_spin)
            form.addRow("Max heavy atoms", self.max_spin)
            form.addRow("Number of molecules", self.count_spin)
            form.addRow(QLabel("<span style='color:#b4b4b4'>"
                               "Size is approximate: whole fragments may "
                               "overshoot the target by a few atoms.</span>"))

            btns = QHBoxLayout()
            self.generate_btn = QPushButton("Generate & Add to Canvas")
            self.generate_btn.setDefault(True)
            self.generate_btn.clicked.connect(self._run)
            self.copy_btn = QPushButton("Copy SMILES")
            self.copy_btn.clicked.connect(self._copy)
            btns.addWidget(self.generate_btn)
            btns.addWidget(self.copy_btn)

            self.result_box = QPlainTextEdit()
            self.result_box.setReadOnly(True)
            self.result_box.setPlaceholderText(
                "Generated SMILES will appear here...")

            layout = QVBoxLayout()
            layout.addLayout(form)
            layout.addLayout(btns)
            layout.addWidget(self.result_box)
            self.setLayout(layout)

        def _current_settings(self):
            return {
                "mode": self.mode_combo.currentText(),
                "nmin": self.min_spin.value(),
                "nmax": max(self.min_spin.value(), self.max_spin.value()),
                "count": self.count_spin.value(),
            }

        def _run(self):
            st = self._current_settings()
            try:
                self.context.set_setting("last_config", st)
            except Exception:
                logging.debug(
                    "RandomMolGen: could not persist dialog settings",
                    exc_info=True)
            smis = _add_molecules(
                self.context, st["mode"], st["nmin"], st["nmax"], st["count"])
            if smis:
                self.result_box.setPlainText("\n".join(smis))

        def _copy(self):
            text = self.result_box.toPlainText()
            if not text:
                return
            app = QApplication.instance()
            if app is not None:
                app.clipboard().setText(text)
                self.context.show_status_message(
                    f"Copied {len(text.splitlines())} SMILES string(s).", 2500)

else:  # pragma: no cover - only when Qt bindings are not importable

    class GeneratorDialog:
        def __init__(self, *args, **kwargs):
            raise RuntimeError(
                "PyQt6 is not available in this host; the options "
                "dialog cannot be opened. The quick-menu actions still "
                "work.")


# ---------------------------------------------------------------------------
# Plugin entry points
# ---------------------------------------------------------------------------

def _add_molecules(context, mode, nmin, nmax, count):
    """Shared action body: generate, add to canvas, checkpoint, refresh."""
    try:
        smis = generate_batch(mode, nmin, nmax, count)
    except Exception:
        logging.exception("RandomMolGen: generation failed")
        context.show_status_message(
            "Generation failed - see the log for details.", 5000)
        return []

    if not smis:
        context.show_status_message(
            "No valid molecule generated - try a larger size range.", 5000)
        return []

    ok = 0
    for s in smis:
        try:
            context.load_from_smiles(s)
            ok += 1
        except Exception:
            logging.exception("RandomMolGen: failed to add %s to canvas", s)

    if ok:
        context.push_undo_checkpoint()
        try:
            context.fit_2d_view()
        except Exception:
            logging.debug("RandomMolGen: fit_2d_view unavailable",
                          exc_info=True)
        context.show_status_message(f"Added {ok} random molecule(s).", 3000)
    logger.info("RandomMolGen: added %d/%d molecules (mode=%s)",
                ok, len(smis), mode)
    return smis


def _load_saved_config(context, defaults):
    cfg = dict(defaults)
    try:
        saved = context.get_setting("last_config") or {}
        for key in ("mode", "nmin", "nmax", "count"):
            if key in saved:
                cfg[key] = saved[key]
    except Exception:
        logging.debug("RandomMolGen: could not read saved config",
                      exc_info=True)
    return cfg


def _quick(context, mode):
    cfg = _load_saved_config(context,
                             {"nmin": 8, "nmax": 18, "count": 1})
    _add_molecules(context, mode,
                   int(cfg["nmin"]), int(cfg["nmax"]), int(cfg["count"]))


def _open_dialog(context):
    # Singleton first: a previously created window is shown even if the
    # Qt import below would fail (e.g. host without PyQt6).
    win = context.get_window("generator_dialog")
    if win is not None:
        win.show()
        win.raise_()
        return

    if not _HAS_QT:
        context.show_status_message(
            "Qt bindings not found - dialog unavailable "
            "(quick menu actions still work).", 5000)
        return

    dlg = GeneratorDialog(context, context.get_main_window())
    context.register_window("generator_dialog", dlg)
    dlg.show()


def initialize(context):
    """MoleditPy entry point (Plugin API v4.0)."""
    logger.info("Random Molecule Generator v%s initialized (%s)",
                PLUGIN_VERSION, "rdkit ok" if Chem else "RDKIT MISSING")

    context.add_menu_action(
        "Tools/Random Molecule/Generator...",
        lambda: _open_dialog(context),
        shortcut="Ctrl+Shift+R")
    context.add_menu_action(
        "Tools/Random Molecule/Quick: Drug-like",
        lambda: _quick(context, "Drug-like"))
    context.add_menu_action(
        "Tools/Random Molecule/Quick: Cyclic",
        lambda: _quick(context, "Cyclic"))
    context.add_menu_action(
        "Tools/Random Molecule/Quick: Chain",
        lambda: _quick(context, "Chain"))
    context.add_toolbar_action(
        lambda: _quick(context, "Drug-like"),
        "Random Mol",
        tooltip="Add one random drug-like molecule "
                "(settings saved in the dialog)")


# ---------------------------------------------------------------------------
# Standalone demo: python random_mol_generator.py
# ---------------------------------------------------------------------------
if __name__ == "__main__":  # pragma: no cover
    rng = random.Random(42)  # nosec B311 (demo reproducibility, not security)
    for mode in MODES:
        print(f"== {mode} ==")
        for s in generate_batch(mode, 8, 16, 3, rng=rng):
            print("  ", s)
