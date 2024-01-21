import unittest

import pandas as pd
import numpy as np

from experiments.biosynthetic_distance import pairs_in_chain
from experiments.metacyc_extract import (
    _row_precursors,
    _row_products,
    _is_next_rxn,
    get_overlap_mol,
    get_beginnings,
    get_child_rxns,
    make_tree,
    tree_to_chains,
    __choose_first_mols,
    _pw_to_chains_mols,
    chains_per_pathway,
)

one_pathway = [
    {
        "reaction_id": "GLUC1PURIDYLTRANS-RXN",
        "pathway_id": "PWY-6527",
        "left": "GLC-1-P",
        "direction": "L2R",
        "right": "CPD-12575",
    },
    {
        "reaction_id": "GALACTOKIN-RXN",
        "pathway_id": "PWY-6527",
        "left": "ALPHA-D-GALACTOSE",
        "direction": "L2R",
        "right": "GALACTOSE-1P",
    },
    {
        "reaction_id": "GALACTURIDYLYLTRANS-RXN",
        "pathway_id": "PWY-6527",
        "left": "GALACTOSE-1P CPD-12575",
        "direction": "L2R",
        "right": "GLC-1-P CPD-14553",
    },
    {
        "reaction_id": "RXN-11502",
        "pathway_id": "PWY-6527",
        "left": "CPD-1099",
        "direction": "L2R",
        "right": "ALPHA-D-GALACTOSE",
    },
    {
        "reaction_id": "RXN-11501",
        "pathway_id": "PWY-6527",
        "left": "CPD-170",
        "direction": "L2R",
        "right": "ALPHA-D-GALACTOSE CPD-1099",
    },
    {
        "reaction_id": "UDPGLUCEPIM-RXN",
        "pathway_id": "PWY-6527",
        "left": "CPD-14553",
        "direction": "L2R",
        "right": "CPD-12575",
    },
    {
        "reaction_id": "UTPHEXPURIDYLYLTRANS-RXN",
        "pathway_id": "PWY-6527",
        "left": "GALACTOSE-1P",
        "direction": "L2R",
        "right": "CPD-14553",
    },
]
df = pd.DataFrame.from_records(one_pathway)


class TestMetacycExtract(unittest.TestCase):
    class TestPathwayToChain:
        def test_row_precursors(self):
            self.assertEqual(
                _row_precursors(df[df["reaction_id"] == "GALACTURIDYLYLTRANS-RXN"]),
                ["GALACTOSE-1P", "CPD-12575"],
            )

        def test_row_products(self):
            self.assertEqual(
                _row_products(df[df["reaction_id"] == "GALACTURIDYLYLTRANS-RXN"]),
                ["GLC-1-P", "CPD-14553"],
            )

        def test_is_next_rxn(self):
            # assert true
            self.assertTrue(
                _is_next_rxn(
                    df[df["reaction_id"] == "GALACTOKIN-RXN"],
                    df[df["reaction_id"] == "GALACTURIDYLYLTRANS-RXN"],
                )
            )

            self.assertTrue(
                _is_next_rxn(
                    df[df["reaction_id"] == "GALACTURIDYLYLTRANS-RXN"],
                    df[df["reaction_id"] == "UDPGLUCEPIM-RXN"],
                )
            )

            # assert false
            self.assertFalse(
                _is_next_rxn(
                    df[df["reaction_id"] == "GALACTURIDYLYLTRANS-RXN"],
                    df[df["reaction_id"] == "RXN-11502"],
                )
            )

        def test_get_overlap_mol(self):
            self.assertEqual(
                get_overlap_mol(
                    df[df["reaction_id"] == "GALACTURIDYLYLTRANS-RXN"],
                    df[df["reaction_id"] == "UDPGLUCEPIM-RXN"],
                ),
                "CPD-14553",
            )
            self.assertEqual(
                get_overlap_mol(
                    df[df["reaction_id"] == "UDPGLUCEPIM-RXN"],
                    df[df["reaction_id"] == "GALACTURIDYLYLTRANS-RXN"],
                ),
                "CPD-12575",
            )

        def test_get_beginnings(self):
            beginnings = get_beginnings(df)
            self.assertEqual(beginnings, ["RXN-11501"])

        def test_get_child_rxns(self):
            # case with child rxns
            parent_rxn = "RXN-11501"
            unused_rxns = df["reaction_id"].tolist()
            unused_rxns.remove(parent_rxn)
            child_rxns = get_child_rxns("RXN-11501", df, unused_rxns=unused_rxns)
            self.assertEqual(child_rxns, ["GALACTOKIN-RXN", "RXN-11502"])
            # case with no child rxns left
            parent_rxn = "RXN-11502"
            unused_rxns = df["reaction_id"].tolist()
            unused_rxns.remove(parent_rxn)
            unused_rxns.remove("GALACTOKIN-RXN")
            child_rxns = get_child_rxns("RXN-11502", df, unused_rxns=unused_rxns)
            self.assertEqual(child_rxns, [])

        def test_make_tree(self):
            # start lower down tree
            tree = make_tree("RXN-11501", df)
            self.assertEqual(
                tree["RXN-11501"]["children"], ["GALACTOKIN-RXN", "RXN-11502"]
            )
            self.assertEqual(tree["RXN-11502"]["children"], [])


class TestBiosyntheticDistance:
    def test_pairs_in_chain(self):
        # test that pairs_in_chain returns the correct pairs
        # for a given chain
        chain = [
            "CHEBI:1234",
            "CHEBI:5678",
            "CHEBI:9012",
            "CHEBI:3456",
            "CHEBI:7890",
        ]
        expected_pairs = [
            ("CHEBI:1234", "CHEBI:5678"),
            ("CHEBI:5678", "CHEBI:9012"),
            ("CHEBI:9012", "CHEBI:3456"),
            ("CHEBI:3456", "CHEBI:7890"),
        ]
        pairs = pairs_in_chain(chain, separation=1)
        self.assertEqual(pairs, expected_pairs)

        chain = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        expected_pairs_2 = [
            [1, 3],
            [2, 4],
            [3, 5],
            [4, 6],
            [5, 7],
            [6, 8],
            [7, 9],
            [8, 10],
        ]
        pairs_2 = pairs_in_chain(chain, separation=2)
        self.assertEqual(pairs_2, expected_pairs_2)
