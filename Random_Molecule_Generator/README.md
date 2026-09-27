# Random Molecule Generator — MoleditPy Plugin
by Qwen3.8:27B

Generates random organic molecules and adds them to the MoleditPy 2D
canvas. Written for the **MoleditPy Plugin API v4.0** (single-file plugin,
no external assets).

| | |
|---|---|
| **File** | `random_mol_generator.py` (the plugin) |
| **Version** | 1.0.0 |
| **Author** | Author |
| **Requires** | MoleditPy ≥ 4.0, RDKit (host-provided), PyQt6 (host-provided) |
| **Tests** | `test_random_mol_generator.py` (11 tests, all passing) |

---

## What it does

Three generation strategies, selectable from the dialog:

| Mode | Strategy |
|---|---|
| **Drug-like** | Common drug-core rings (benzene, piperidine, piperazine, indole, naphthalene, quinoline, …) grown with 0–6 curated, weighted functional groups (methyl, amide, ester, amine, phenyl, heterocycles, …) attached at random positions. |
| **Cyclic** | A freshly generated ring — 3–8 membered, 0–3 heteroatoms (N/O/S), optionally aromatic — grown with random fragments. |
| **Chain** | Monte-Carlo random walk over the SMILES grammar: branched chains (≤3 deep), accidental rings (size ≥3), terminal carbonyls and halogens. |

Every candidate is **validated by RDKit**; invalid attempts are retried
transparently (a 20-molecule batch typically finishes in well under a
second). Output sizes respect the requested min/max within a small
overshoot margin (whole fragments are attached, not clipped).

## Installation

Copy or drop `random_mol_generator.py` into your plugin directory:

- **Windows**: `C:\Users\<You>\.moleditpy\plugins\`
- **Linux / macOS**: `~/.moleditpy/plugins/`

Then restart MoleditPy or choose **Plugins ▸ Reload All Plugins**.

## Usage

### Menu entries

- `Tools ▸ Random Molecule ▸ Generator…` — options dialog (shortcut **Ctrl+Shift+R**)
- `Tools ▸ Random Molecule ▸ Quick: Drug-like / Cyclic / Chain` — one-click generation using the size range saved from the dialog
- Plugin toolbar button **“Random Mol”** — quick drug-like generation

### Options dialog

- **Mode** — one of the three strategies above
- **Min / Max heavy atoms** — target size range (size is approximate)
- **Number of molecules** — 1–25 distinct structures per batch
- **Generate & Add to Canvas** — creates the molecules (undoable)
- **Copy SMILES** — copies the last batch to the clipboard

Molecules are added with `context.load_from_smiles()`, so you can select,
edit, or delete them in the canvas as usual. Each batch is one undo step.
Your last dialog settings are persisted in the MoleditPy settings
(namespaced under the plugin name).

### Example output

```
== Drug-like ==
   O=C(O)C1CCCCC1
   O=C(O)C1NCCCN1C(=O)O
   c1ccc2[nH]ccc2c1
== Cyclic ==
   CCC1CCCCO1
   CCSCN(C)C1CCCCC1CC
   O=CCC1CCCCCC1
== Chain ==
   COC1CSCN[SH](N)N(N)C1F
   CNN1C[SH](Cl)N(C)[SH](SNOC)CO1
   CCN1C(CCOCCCSC)SC[SH]1N
```

## How it works

- **Fragment attachment** is done at the **molecular-graph level**
  (`RWMol`), not by splicing SMILES text — this is robust across RDKit
  versions (some builds, e.g. 2026.03.x, no longer accept atom-map
  tokens in `MolFromSmiles`).
- **Random-walk mode** keeps the SMILES grammar valid by construction
  (balanced branches, single-use ring digits, ring size ≥3); RDKit
  sanitization filters the rest.
- The generator core is **pure Python + RDKit with no Qt**, so it runs
  and is testable headless; Qt is used only for the optional dialog.
- All callbacks run on the GUI thread and complete in milliseconds.

## Development & testing

The test suite mocks `PluginContext` (per the v4.0 manual §12) and covers
generation validity, size bounds, deduplication, determinism, the
`initialize()` wiring, and the undo/checkpoint contract:

```bash
pip install rdkit pytest
pytest test_random_mol_generator.py -q        # 11 passed

python random_mol_generator.py                # standalone demo (no Qt needed)
```

## Security notes

- **No file, network, or subprocess I/O** — the plugin never opens files
  (settings go through the host's `get_setting`/`set_setting`), makes no
  network calls, and never executes external code.
- **No `eval`/`exec`/`pickle`/`ctypes`**, no dynamic imports.
- **No embedded secrets, tokens, or credentials**; no hardcoded private
  paths (the docstring mentions the standard `~/.moleditpy/plugins/`
  install location, as documented in the official manual).
- Static scan results: `bandit 1.9.4` → **0 issues** (the two
  `try/except` import fallbacks log at DEBUG level; the 3 `random.Random`
  usages are intentionally non-cryptographic and marked `# nosec B311`
  with an explanatory comment).
- The only imports are the standard library (`logging`, `random`),
  `rdkit` and `PyQt6` — both provided by the host application.

## Plugin API v4.0 surface used

`initialize(context)` · `add_menu_action` · `add_toolbar_action` ·
`show_status_message` · `load_from_smiles` · `push_undo_checkpoint` ·
`fit_2d_view` · `get_setting` / `set_setting` · `register_window` /
`get_window` (singleton dialog) · `get_main_window`

No direct `MainWindow` access beyond the stable `PluginContext` proxy.
