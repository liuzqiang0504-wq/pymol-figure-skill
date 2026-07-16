"""Regression test for tetrahedral stereochemistry in the 2D renderer."""

import unittest

from rdkit import Chem

import generate_2d_interactions as renderer


class StereochemistryRenderingTest(unittest.TestCase):
    def test_prepared_ligand_has_solid_and_hashed_wedges(self):
        ligand = Chem.MolFromSmiles(
            "N[C@@H](C)C(=O)N[C@@H](Cc1ccccc1)C(=O)O")
        self.assertIsNotNone(ligand)

        drawing_mol, _ = renderer._prepare_ligand(ligand)
        centers = Chem.FindMolChiralCenters(
            drawing_mol, includeUnassigned=True, includeCIP=True)
        directions = {
            bond.GetBondDir()
            for bond in drawing_mol.GetBonds()
            if bond.GetBondDir() != Chem.BondDir.NONE
        }

        self.assertEqual([(1, "S"), (6, "S")], centers)
        self.assertIn(Chem.BondDir.BEGINWEDGE, directions)
        self.assertIn(Chem.BondDir.BEGINDASH, directions)

        report = renderer.verify_stereochemistry(ligand, drawing_mol)
        self.assertTrue(report["passed"], report["errors"])
        self.assertEqual({1: "S", 6: "S"}, report["assigned_centers"])


class LayoutQualityTest(unittest.TestCase):
    def test_quality_gate_rejects_overlapping_annotations(self):
        report = renderer._verify_annotation_layout(
            residue_boxes=[(10, 10, 90, 40), (70, 20, 150, 50)],
            distance_boxes=[(48, 48, 92, 68)],
            ligand_box=(40, 40, 120, 120),
            fan_boxes=[],
            width=200,
            height=160,
        )
        self.assertFalse(report["passed"])
        self.assertTrue(any("overlap" in error for error in report["errors"]))


if __name__ == "__main__":
    unittest.main()
