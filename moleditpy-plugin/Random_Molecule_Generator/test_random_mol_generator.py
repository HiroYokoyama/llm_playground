"""Tests for the Random Molecule Generator plugin (MoleditPy v4.0 API).

Run with: pytest test_random_mol_generator.py  (or: python -m unittest -v)
"""

import random
import unittest

from unittest.mock import MagicMock

from rdkit import Chem

import random_mol_generator as g


class TestGenerationCore(unittest.TestCase):
    def test_all_modes_produce_valid_molecules(self):
        for mode in g.MODES:
            with self.subTest(mode=mode):
                smis = g.generate_batch(mode, 4, 14, 10)
                self.assertGreater(len(smis), 0, mode)
                for s in smis:
                    m = Chem.MolFromSmiles(s)
                    self.assertIsNotNone(m, f"invalid SMILES: {s}")
                    self.assertGreaterEqual(m.GetNumHeavyAtoms(), 4)

    def test_no_duplicates_in_batch(self):
        smis = g.generate_batch("Drug-like", 6, 14, 20)
        self.assertEqual(len(smis), len(set(smis)))

    def test_size_bounds_respected(self):
        for s in g.generate_batch("Cyclic", 3, 6, 10):
            m = Chem.MolFromSmiles(s)
            self.assertGreaterEqual(m.GetNumHeavyAtoms(), 3)

    def test_deterministic_with_seed(self):
        a = g.generate_batch("Chain", 6, 10, 5, rng=random.Random(7))  # nosec B311 (determinism test, not security)
        b = g.generate_batch("Chain", 6, 10, 5, rng=random.Random(7))  # nosec B311 (determinism test, not security)
        self.assertEqual(a, b)

    def test_mode_aliases(self):
        self.assertEqual(g._mode_key("drug"), "Drug-like")
        self.assertEqual(g._mode_key("ring"), "Cyclic")
        self.assertEqual(g._mode_key("random walk"), "Chain")

    def test_attach_fragment_roundtrip(self):
        rng = random.Random(1)  # nosec B311 (test determinism, not security)
        s = _attach_once("c1ccccc1", rng)
        self.assertIsNotNone(s)
        m = Chem.MolFromSmiles(s)
        self.assertGreater(m.GetNumHeavyAtoms(), 6)

    def test_attach_fragment_ring_digit_remap(self):
        # core already uses ring digit 1 -> fragment must not clash
        ok = 0
        for i in range(20):
            rng = random.Random(i)  # nosec B311 (test determinism, not security)
            out = g._attach_fragment("C1CCCCC1", "C2=CC=CC=C2", rng)
            if out is not None:
                ok += 1
                self.assertIsNotNone(Chem.MolFromSmiles(out), out)
        self.assertGreater(ok, 0)


def _attach_once(seed, rng):
    return g._attach_fragment(seed, rng.choice([x[0] for x in g._FRAGMENTS]), rng)


class TestPluginApi(unittest.TestCase):
    """Verify actions behave per the v4.0 manual (mocked PluginContext)."""

    def test_quick_action_calls_plugin_api(self):
        context = MagicMock()
        g._quick(context, "Drug-like")
        # Molecules must be added to the canvas...
        self.assertTrue(context.load_from_smiles.call_count >= 1)
        # ...and the edit must be checkpointed onto the undo stack.
        context.push_undo_checkpoint.assert_called_once()
        # And the view should refit.
        context.fit_2d_view.assert_called_once()

    def test_status_message_on_failure_paths(self):
        context = MagicMock()
        context.get_setting.side_effect = RuntimeError("boom")
        # Saved config read failed -> defaults used; should not raise.
        g._quick(context, "Chain")
        context.show_status_message.assert_called()

    def test_initialize_registers_actions(self):
        context = MagicMock()
        g.initialize(context)
        paths = [c.args[0] for c in context.add_menu_action.call_args_list]
        self.assertTrue(any("Generator..." in p for p in paths))
        self.assertTrue(any("Quick: Drug-like" in p for p in paths))
        self.assertTrue(any("Quick: Cyclic" in p for p in paths))
        self.assertTrue(any("Quick: Chain" in p for p in paths))
        context.add_toolbar_action.assert_called_once()

    def test_open_dialog_singleton(self):
        context = MagicMock()
        fake_win = MagicMock()
        context.get_window.return_value = fake_win
        g._open_dialog(context)
        fake_win.show.assert_called_once()
        context.register_window.assert_not_called()


if __name__ == "__main__":
    unittest.main(verbosity=2)
