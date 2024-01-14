import unittest, logging
from enum import Enum
import subprocess

# set logging level to debug for this test
# logging.basicConfig(level=logging.DEBUG)

from rdkit import Chem
import numpy as np

from biosynfoni.concerto_fp import MolsCollection, Biosynfoni
from biosynfoni.subkeys import defaultVersion, get_smarts
from biosynfoni.moldrawing import draw
from biosynfoni import get_highlight_mapping, draw_with_highlights

general_test_smiles = [
    "CC(=O)CC=O",
    "C1CCCCC1",
    "COc1cc(O)c2c(c1)oc(cc2=O)-c1ccc(OC)c(c1)-c1c(O)cc(O)c2c1oc(cc2=O)-c1ccc(O)cc1",
    "[H][C@]1(CC[C@@H](O)[C@@H](C1)OC)C[C@@H](C)[C@]1([H])CC(=O)[C@H](C)\C=C(C)\[C@@H](O)[C@@H](OC)C(=O)[C@H](C)C[C@H](C)\C=C\C=C\C=C(C)\[C@H](C[C@]2([H])CC[C@@H](C)[C@@](O)(O2)C(=O)C(=O)N2CCCC[C@@]2([H])C(=O)O1)OC",
    "[C@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)(NC(=O)*)CO[C@@H]1O[C@H](CO)[C@H]([C@@H]([C@H]1O)O)O[C@@H]2O[C@H](CO[C@]3(O[C@]([C@@H]([C@H](C3)O)NC(C)=O)([C@@H]([C@@H](CO)O)O)[H])C(=O)O)[C@@H]([C@H](O)[C@H]2O)O",
]
general_test_mols = [Chem.MolFromSmiles(smiles) for smiles in general_test_smiles]
empty_fp = [0 for _ in get_smarts(version=defaultVersion)]


class testCommandLine(unittest.TestCase):
    def test_single_smiles(self):
        for smile in general_test_smiles:
            cmd = ["biosynfoni", smile]
            result = subprocess.run(cmd, capture_output=True)
            # make sure it has no errors
            self.assertEqual(result.returncode, 0)

    def test_multiple_smiles(self):
        cmd = ["biosynfoni"] + general_test_smiles
        result = subprocess.run(cmd, capture_output=True)
        # make sure it has no errors
        self.assertEqual(result.returncode, 0)


# easy functions for testing
def dimensions(lst: list):
    return len(np.array(lst).shape)


def core_type(lst: list):
    return np.dtype(np.array(lst))


class testMolsCollection(unittest.TestCase):
    pass


# class testRepresentationConversion(unittest.TestCase):
#     def test_smiles(self):
#         smiles_list = [""]
#         self.assertEqual(
#             MolsCollection.from_smiles(["C", "CC"]).to_smiles(), ["C", "CC"]
#         )


class testHelperFunctions(unittest.TestCase):
    def test_count_listitems(self):
        a = [0, 1, 2, 3, 4, 5]
        b = [[1, 2, 43, 0], [1], [24], [], [25]]
        c = [[1, 2, 43, 0], [1], [24], [], [25], ""]
        d = [[1, 2, 43, 0], [1], [24], [[], [], [0]], [25], ""]

        from biosynfoni.concerto_fp import _count_listitems as _count

        self.assertEqual(_count(a), 6)
        self.assertEqual(_count(b), 7)
        self.assertEqual(_count(c), 7)
        self.assertEqual(_count(d), 8)


class testBiosynfoni(unittest.TestCase):
    def test_subs_set_retrieval(self):
        for i, sub in enumerate(Biosynfoni(general_test_mols[0]).substructure_set):
            assert isinstance(sub, Chem.Mol), f"sub {i} is not a Chem.Mol"

    def test_substructure_detection(self):
        # check right functioning for easy molecule
        easy_mol = Chem.MolFromSmiles("CC")
        results = Biosynfoni(easy_mol, version_name=defaultVersion).overlap_matches
        # should be list[list[list[int]]]
        assert isinstance(results, list)
        # assert dimensions(results) == 3 #functions not yet ready
        # assert core_type(results) == np.dtype(int)

        for mol in general_test_mols:
            all_detected = Biosynfoni(mol, version_name=defaultVersion).overlap_matches
            # hand_added_biosynfoni = []
            assert all_detected != empty_fp
            # assert all_detected == hand_added_biosynfoni

    def test_substructure_detection_with_overlap(self):
        pass

    def test_substructure_detection_visually(self):
        for i, mol in enumerate(general_test_mols):
            info = get_highlight_mapping(mol)
            svg_str = draw(mol, highlight_atoms_bonds_mappings=info)
            svg_str2 = draw_with_highlights(mol)

            self.assertEqual(svg_str, svg_str2)
            # add mol's smiles to the svg as text
            svg_str = svg_str.replace(
                "</svg>",
                f'<text x="30" y="30" font-size="20" font-family="montserrat">{Chem.MolToSmiles(mol)}</text></svg>',
            )
            with open(f"substructure_detection_test{i}.svg", "w") as f:
                f.write(svg_str)

    def test_coverage_calculation(self):
        for mol in general_test_mols:
            fullblock = Biosynfoni(mol, version_name=defaultVersion)
            coverage_fullblock = fullblock.get_coverage()
            coverage_fullblock = fullblock.coverage
            if fullblock.fingerprint == empty_fp:
                assert coverage_fullblock == 0

            assert coverage_fullblock <= 1
            assert coverage_fullblock >= 0

            coverage_noblock = Biosynfoni(
                mol,
                version_name=defaultVersion,
                intrasub_overlap=True,
                intersub_overlap=True,
            ).get_coverage()
            assert coverage_noblock >= 0
            assert coverage_noblock >= coverage_fullblock

    def test_coverage_calculation_handcontrol_easy(self):
        smiles = "PC"
        mol = Chem.MolFromSmiles(smiles)
        no_overlap = Biosynfoni(mol, version_name=defaultVersion)
        inter_overlap = Biosynfoni(
            mol, version_name=defaultVersion, intersub_overlap=True
        )
        overlap = Biosynfoni(
            mol,
            version_name=defaultVersion,
            intrasub_overlap=True,
            intersub_overlap=True,
        )
        self.assertEqual(no_overlap.get_coverage(), 0.5, f"no_over{smiles}")
        self.assertEqual(inter_overlap.get_coverage(), 0.5, f"inter_over{smiles}")
        self.assertEqual(overlap.get_coverage(), 0.5, f"overlap{smiles}")

    def test_coverage_calculation_handcontrol_medium(self):
        smiles = "PCC(=O)CP"  # P not detected, CC(=O) and CC
        mol = Chem.MolFromSmiles(smiles)

        # # set logging level for concertofp to debug
        # logging.getLogger("biosynfoni.concerto_fp").setLevel(logging.DEBUG)

        no_overlap = Biosynfoni(mol)

        inter_overlap = Biosynfoni(mol, intersub_overlap=True)
        # logging.debug(f"inter_overlap: {inter_overlap.get_coverage()}")
        # logging.debug(f"inter_overlap fp: {inter_overlap.fingerprint}")

        overlap = Biosynfoni(
            mol,
            intrasub_overlap=True,
            intersub_overlap=True,
        )
        self.assertEqual(no_overlap.get_coverage(), (3 / 6), f"no_over{smiles}")
        self.assertEqual(inter_overlap.get_coverage(), (5 / 6), f"inter_over{smiles}")
        self.assertEqual(overlap.get_coverage(), (10 / 6), msg=f"overlap{smiles}")
        # set logging level for concertofp back to normal
        logging.getLogger("biosynfoni.concerto_fp").setLevel(logging.WARNING)

    def test_no_chirality_difference(self):
        for mol in general_test_mols:
            no_chiral = mol.__copy__()
            Chem.RemoveStereochemistry(no_chiral)

            if no_chiral:
                nc_fp = Biosynfoni(no_chiral).fingerprint
                c_fp = Biosynfoni(mol).fingerprint
                self.assertEqual(len(nc_fp), len(c_fp))
                self.assertEqual(
                    Biosynfoni(no_chiral).fingerprint, Biosynfoni(mol).fingerprint
                )
            else:
                logging.warning(f"no_chiral is None for {Chem.MolToSmiles(mol)}")


class testOverlapFiltering(unittest.TestCase):
    def test_intrasub_overlap_filtering(self):
        ethyl_acetyl_overlap = Chem.MolFromSmiles("CC(=O)CC=O")
        mol = ethyl_acetyl_overlap
        fullblock = Biosynfoni(
            mol,
            version_name=defaultVersion,
            intersub_overlap=False,
            intrasub_overlap=False,
        )
        interoverlap = Biosynfoni(
            mol,
            version_name=defaultVersion,
            intersub_overlap=True,
            intrasub_overlap=False,
        )
        alloverlap = Biosynfoni(
            mol,
            version_name=defaultVersion,
            intersub_overlap=True,
            intrasub_overlap=True,
        )

        assert fullblock.fingerprint != interoverlap.fingerprint
        assert alloverlap.fingerprint != interoverlap.fingerprint

    def test_smart_overlap_detection(self):
        mol = Chem.MolFromSmiles("PC1CCCCC1P")  # change to difficult molecule
        # set logging level for concertofp to debug
        logging.getLogger("biosynfoni.concerto_fp").setLevel(logging.DEBUG)
        fp = Biosynfoni(mol).fingerprint
        self.assertEqual(Biosynfoni(mol).get_coverage(), 6 / 8, msg=f"{fp}")
        self.assertEqual(Biosynfoni(mol, intrasub_overlap=True).get_coverage(), 12 / 8)
        # self.assertEqual(Biosynfoni(mol, intersub_overlap=True).get_coverage(), 1.0)


if __name__ == "__main__":
    unittest.main()
