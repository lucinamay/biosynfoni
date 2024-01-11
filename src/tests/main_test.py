import unittest

from rdkit import Chem
from biosynfoni.concerto_fp import MolsCollection, Biosynfoni
from biosynfoni.concerto_fp import defaultVersion as DefaultVersion

general_test_mols = [
    Chem.MolFromSmiles("CC(=O)CC=O"),
    Chem.MolFromSmiles("C1CCCCC1"),
]
empty_fp = [0 for _ in DefaultVersion]


class testMolsCollection(unittest.TestCase):
    pass


class testRepresentationConversion(unittest.TestCase):
    def test_smiles(self):
        smiles_list = [""]
        self.assertEqual(
            MolsCollection.from_smiles(["C", "CC"]).to_smiles(), ["C", "CC"]
        )


class testBiosynfoni(unittest.TestCase):
    def test_substructure_detection(self):
        for mol in general_test_mols:
            all_detected = Biosynfoni(mol, DefaultVersion).overlap_matches
            # hand_added_biosynfoni = []
            assert all_detected != empty_fp
            # assert all_detected == hand_added_biosynfoni

    def test_coverage_calculation(self):
        for mol in general_test_mols:
            fullblock = Biosynfoni(mol, DefaultVersion)
            coverage_fullblock = fullblock.coverage
            if fullblock.fingerprint == empty_fp:
                assert coverage_fullblock == 0

            assert coverage_fullblock <= 1
            assert coverage_fullblock >= 0

            coverage_noblock = Biosynfoni(
                mol, DefaultVersion, intrasub_overlap=True, intersub_overlap=True
            ).coverage
            assert coverage_noblock >= 0
            assert coverage_noblock >= coverage_fullblock

        mol = Chem.MolFromSmiles("PC")
        assert Biosynfoni(mol, DefaultVersion).coverage == 0.5
        assert Biosynfoni(mol, DefaultVersion, intersub_overlap=True).coverage == 0.5
        assert (
            Biosynfoni(
                mol, DefaultVersion, intrasub_overlap=True, intersub_overlap=True
            ).coverage
            == 0.5
        )

        mol = Chem.MolFromSmiles("PCC(=O)CP")
        self.assertEqual(Biosynfoni(mol, DefaultVersion).coverage, 0.5)
        self.assertEqual(
            Biosynfoni(mol, DefaultVersion, intrasub_overlap=True).coverage, 0.5
        )

        mol = Chem.MolFromSmiles("C1CCCCC1")
        assert Biosynfoni(mol, DefaultVersion).coverage == 1.0


class testOverlapFiltering(unittest.TestCase):
    def test_intrasub_overlap_filtering(self):
        ethyl_acetyl_overlap = Chem.MolFromSmiles("CC(=O)CC=O")
        mol = ethyl_acetyl_overlap
        fullblock = Biosynfoni(
            mol,
            fp_version=DefaultVersion,
            intersub_overlap=False,
            intrasub_overlap=False,
        )
        interoverlap = Biosynfoni(
            mol,
            fp_version=DefaultVersion,
            intersub_overlap=True,
            intrasub_overlap=False,
        )
        alloverlap = Biosynfoni(
            mol,
            fp_version=DefaultVersion,
            intersub_overlap=True,
            intrasub_overlap=True,
        )

        assert fullblock.fingerprint != interoverlap.fingerprint
        assert alloverlap.fingerprint != interoverlap.fingerprint


if __name__ == "__main__":
    unittest.main()
