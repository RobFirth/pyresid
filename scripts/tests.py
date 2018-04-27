# !/usr/bin/env python3
"""

"""

import os
import unittest
import json
import io
import pyresid as pyre

from collections import OrderedDict
import spacy
nlp = spacy.load('en')

if __name__ == "__main__":
    this_dir = os.path.abspath(os.path.join(__file__, os.pardir))
else:
    this_dir = os.curdir

def deep_sort(obj):
    """
    Recursively sort list or dict nested lists - https://stackoverflow.com/a/27949519
    To enable the comparison of OrderedDicts
    """

    if isinstance(obj, dict):
        _sorted = {}
        for key in sorted(obj):
            _sorted[key] = deep_sort(obj[key])

    elif isinstance(obj, list):
        new_list = []
        for val in obj:
            new_list.append(deep_sort(val))
        _sorted = sorted(new_list)

    else:
        _sorted = obj

    return _sorted

class TestClass(unittest.TestCase):
    """
    Class for testing pyresid
    """

    ## Residue identification pyre._identify_residues

    slashed_res_sent = str(
        "Indeed, the guanidyl side chain of arginine is wedged between two conserved aromatic residues (Tyr33/39 and " +
        "Trp71/77 in OccJ/NocT) and points toward the opening of the cleft by making six hydrogen bonds with the" +
        "conserved side chains of residues Glu30/36 and Gln159/165 and the carbonyl of Ala88/94.")

    mixed_res_sent = str(
        "Its carboxyl moiety makes a salt-bridge with the conserved Arg96/102 and three hydrogen " +
        "bonds with the Ser91 side chain (the corresponding residue in NocT is Gly97) and the amide " +
        "NH protons of Ser91/Gly97 and Thr163/Ser169")

    mixed_res_parenthesis_sent = str(
        "Its carboxyl moiety makes a salt-bridge with the conserved Arg(96)/(102) and three hydrogen " +
        "bonds with the Ser(91) side chain (the corresponding residue in NocT is Gly(97)) and the amide " +
        "NH protons of Ser(91)/Gly(97) and Thr(163)/Ser(169)")

    mixed_res_hyphen_sent = str(
        "Its carboxyl moiety makes a salt-bridge with the conserved Arg-96 and three hydrogen " +
        "bonds with the Ser-91 side chain (the corresponding residue in NocT is Gly-97) and the amide " +
        "NH protons of Ser-91/Gly-97 and Thr-163/Ser-169")

    grammatical_sent = str(
        "Although the mutation of a single amino acid in NocT and LAO can modify their ligand " +
        "selectivity to transform them into octopine selective and non-selective PBPs, respectively, " +
        "our work shows that OccJ is evolved for binding octopine and octopinic acid with high " +
        "affinity in nanomolar range mainly due to the presence of a serine at position 91.")

    grammatical_sent_end = str(
        "Although the mutation of a single amino acid in NocT and LAO can modify their ligand " +
        "selectivity to transform them into octopine selective and non-selective PBPs, respectively, " +
        "our work shows that OccJ is evolved for binding octopine and octopinic acid with high " +
        "affinity in nanomolar range mainly due to the presence of a serine at position 91")

    grammatical_sent_start = str(
        "aspartate at position 202 in the OccJ-octopine complex was modelled and its carboxylate group is " +
        "free and positioned far from the pyruvate to be compatible with a bound octopine.")

    grammatical_sent_multiple = str(
        "There are 27 strictly conserved residues, including a cysteine residue at position 275, earlier identified " +
        "by chemical modification as the expected catalytic residue of the second half-reaction (conversion of " +
        "UDP-aldehydoglucose to UDP-glucuronic acid), and 2 lysine residues, at positions 219 and 338, one of which " +
        "may be the expected catalytic residue for the first half-reaction (conversion of UDP-glucose to " +
        "UDP-aldehydoglucose) . ")

    long_dashed_sent = str(
        "This is a really long residue Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202 " + \
        "it is really long")

    long_dashed_sent_start = str(
        "Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202 is a really long residue, " + \
        "it is really long and it is at the start of the string.")

    long_dashed_pos_dashed_sent = str(
        "a really long residue mention is Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202," + \
        "it is really long and it is in the middle of the string.")

    long_dashed_pos_dashed_sent_start = str(
        "Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202 is a really long residue, " + \
        "it is really long and it is at the start of the string.")

    def test_identify_residues_on_slashed_testsent(self):
        matches = pyre.identify_residues_ravikumar(self.slashed_res_sent)

        test_set = [match.__dict__ for match in matches]
        truth = [{'aminoacid': 'Tyrosine',
                  'end': 103,
                  'position': ['33', '39'],
                  'start': 95,
                  'string': 'Tyr33/39',
                  'threeletter': 'Tyr'},
                 {'aminoacid': 'Tryptophan',
                  'end': 116,
                  'position': ['71', '77'],
                  'start': 108,
                  'string': 'Trp71/77',
                  'threeletter': 'Trp'},
                 {'aminoacid': 'Glutamic acid (Glutamate)',
                  'end': 253,
                  'position': ['30', '36'],
                  'start': 245,
                  'string': 'Glu30/36',
                  'threeletter': 'Glu'},
                 {'aminoacid': 'Glutamine',
                  'end': 268,
                  'position': ['159', '165'],
                  'start': 258,
                  'string': 'Gln159/165',
                  'threeletter': 'Gln'},
                 {'aminoacid': 'Alanine',
                  'end': 297,
                  'position': ['88', '94'],
                  'start': 289,
                  'string': 'Ala88/94',
                  'threeletter': 'Ala'}]

        # self.assertEqual(deep_sort(test_set), deep_sort(truth))
        self.assertEqual(truth, test_set)

    def test_identify_residues_on_slashed_testsent_decompose(self):
        matches = pyre.identify_residues_ravikumar(self.slashed_res_sent)

        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]
        truth = [{'aminoacid': 'Tyrosine',
                  'end': 103,
                  'position': '33',
                  'residue': 'Tyr33',
                  'start': 95,
                  'string': 'Tyr33/39',
                  'threeletter': 'Tyr'},
                 {'aminoacid': 'Tyrosine',
                  'end': 103,
                  'position': '39',
                  'residue': 'Tyr39',
                  'start': 95,
                  'string': 'Tyr33/39',
                  'threeletter': 'Tyr'},
                 {'aminoacid': 'Tryptophan',
                  'end': 116,
                  'position': '71',
                  'residue': 'Trp71',
                  'start': 108,
                  'string': 'Trp71/77',
                  'threeletter': 'Trp'},
                 {'aminoacid': 'Tryptophan',
                  'end': 116,
                  'position': '77',
                  'residue': 'Trp77',
                  'start': 108,
                  'string': 'Trp71/77',
                  'threeletter': 'Trp'},
                 {'aminoacid': 'Glutamic acid (Glutamate)',
                  'end': 253,
                  'position': '30',
                  'residue': 'Glu30',
                  'start': 245,
                  'string': 'Glu30/36',
                  'threeletter': 'Glu'},
                 {'aminoacid': 'Glutamic acid (Glutamate)',
                  'end': 253,
                  'position': '36',
                  'residue': 'Glu36',
                  'start': 245,
                  'string': 'Glu30/36',
                  'threeletter': 'Glu'},
                 {'aminoacid': 'Glutamine',
                  'end': 268,
                  'position': '159',
                  'residue': 'Gln159',
                  'start': 258,
                  'string': 'Gln159/165',
                  'threeletter': 'Gln'},
                 {'aminoacid': 'Glutamine',
                  'end': 268,
                  'position': '165',
                  'residue': 'Gln165',
                  'start': 258,
                  'string': 'Gln159/165',
                  'threeletter': 'Gln'},
                 {'aminoacid': 'Alanine',
                  'end': 297,
                  'position': '88',
                  'residue': 'Ala88',
                  'start': 289,
                  'string': 'Ala88/94',
                  'threeletter': 'Ala'},
                 {'aminoacid': 'Alanine',
                  'end': 297,
                  'position': '94',
                  'residue': 'Ala94',
                  'start': 289,
                  'string': 'Ala88/94',
                  'threeletter': 'Ala'}]
        self.assertEqual(truth, test_set)

    def test_identify_residues_on_mixed_testsent(self):
        test_set = [match.__dict__ for match in pyre.identify_residues_ravikumar(self.mixed_res_sent)]
        truth = [{'aminoacid': 'Arginine',
                  'end': 68,
                  'position': ['96', '102'],
                  'start': 59,
                  'string': 'Arg96/102',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Serine',
                  'end': 108,
                  'position': ['91'],
                  'start': 103,
                  'string': 'Ser91',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Glycine',
                  'end': 163,
                  'position': ['97'],
                  'start': 158,
                  'string': 'Gly97',
                  'threeletter': 'Gly'},
                 {'aminoacid': ['Ser', 'Gly'],
                  'end': 204,
                  'position': ['91', '97'],
                  'start': 193,
                  'string': 'Ser91/Gly97',
                  'threeletter': ['Ser', 'Gly']},
                 {'aminoacid': ['Thr', 'Ser'],
                  'end': 222,
                  'position': ['163', '169'],
                  'start': 209,
                  'string': 'Thr163/Ser169',
                  'threeletter': ['Thr', 'Ser']}]
        self.assertEqual(truth, test_set)

    def test_identify_residues_on_mixed_testsent_decompose(self):
        matches = pyre.identify_residues_ravikumar(self.mixed_res_sent)
        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]
        truth = [{'aminoacid': 'Arginine',
                  'end': 68,
                  'position': '96',
                  'residue': 'Arg96',
                  'start': 59,
                  'string': 'Arg96/102',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Arginine',
                  'end': 68,
                  'position': '102',
                  'residue': 'Arg102',
                  'start': 59,
                  'string': 'Arg96/102',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Serine',
                  'end': 108,
                  'position': 91,
                  'residue': 'Ser91',
                  'start': 103,
                  'string': 'Ser91',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Glycine',
                  'end': 163,
                  'position': 97,
                  'residue': 'Gly97',
                  'start': 158,
                  'string': 'Gly97',
                  'threeletter': 'Gly'},
                 {'aminoacid': 'Serine',
                  'end': 204,
                  'position': '91',
                  'residue': 'Ser91',
                  'start': 193,
                  'string': 'Ser91/Gly97',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Glycine',
                  'end': 204,
                  'position': '97',
                  'residue': 'Gly97',
                  'start': 193,
                  'string': 'Ser91/Gly97',
                  'threeletter': 'Gly'},
                 {'aminoacid': 'Threonine',
                  'end': 222,
                  'position': '163',
                  'residue': 'Thr163',
                  'start': 209,
                  'string': 'Thr163/Ser169',
                  'threeletter': 'Thr'},
                 {'aminoacid': 'Serine',
                  'end': 222,
                  'position': '169',
                  'residue': 'Ser169',
                  'start': 209,
                  'string': 'Thr163/Ser169',
                  'threeletter': 'Ser'}]

        self.assertEqual(truth, test_set)

    def test_identify_residues_on_mixed_parenthesis_testsent_decompose(self):
        matches = pyre.identify_residues_ravikumar(self.mixed_res_parenthesis_sent)

        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]

        truth = [{'aminoacid': 'Arginine',
                  'end': 72,
                  'position': '96',
                  'residue': 'Arg96',
                  'start': 59,
                  'string': 'Arg(96)/(102)',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Arginine',
                  'end': 72,
                  'position': '102',
                  'residue': 'Arg102',
                  'start': 59,
                  'string': 'Arg(96)/(102)',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Serine',
                  'end': 114,
                  'position': 91,
                  'residue': 'Ser91',
                  'start': 107,
                  'string': 'Ser(91)',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Glycine',
                  'end': 171,
                  'position': 97,
                  'residue': 'Gly97',
                  'start': 164,
                  'string': 'Gly(97)',
                  'threeletter': 'Gly'},
                 {'aminoacid': 'Serine',
                  'end': 216,
                  'position': '91',
                  'residue': 'Ser91',
                  'start': 201,
                  'string': 'Ser(91)/Gly(97)',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Glycine',
                  'end': 216,
                  'position': '97',
                  'residue': 'Gly97',
                  'start': 201,
                  'string': 'Ser(91)/Gly(97)',
                  'threeletter': 'Gly'},
                 {'aminoacid': 'Threonine',
                  'end': 238,
                  'position': '163',
                  'residue': 'Thr163',
                  'start': 221,
                  'string': 'Thr(163)/Ser(169)',
                  'threeletter': 'Thr'},
                 {'aminoacid': 'Serine',
                  'end': 238,
                  'position': '169',
                  'residue': 'Ser169',
                  'start': 221,
                  'string': 'Thr(163)/Ser(169)',
                  'threeletter': 'Ser'}]

        self.assertEqual(truth, test_set)

    def test_identify_residues_on_mixed_parenthesis_testsent(self):
        matches = pyre.identify_residues_ravikumar(self.mixed_res_parenthesis_sent)

        test_set = [match.__dict__ for match in matches]
        truth = [{'aminoacid': 'Arginine',
                  'end': 72,
                  'position': ['96', '102'],
                  'start': 59,
                  'string': 'Arg(96)/(102)',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Serine',
                  'end': 114,
                  'position': ['91'],
                  'start': 107,
                  'string': 'Ser(91)',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Glycine',
                  'end': 171,
                  'position': ['97'],
                  'start': 164,
                  'string': 'Gly(97)',
                  'threeletter': 'Gly'},
                 {'aminoacid': ['Ser', 'Gly'],
                  'end': 216,
                  'position': ['91', '97'],
                  'start': 201,
                  'string': 'Ser(91)/Gly(97)',
                  'threeletter': ['Ser', 'Gly']},
                 {'aminoacid': ['Thr', 'Ser'],
                  'end': 238,
                  'position': ['163', '169'],
                  'start': 221,
                  'string': 'Thr(163)/Ser(169)',
                  'threeletter': ['Thr', 'Ser']}]
        self.assertEqual(truth, test_set)

    def test_identify_residues_on_mixed_hyphen_testsent(self):
        matches = pyre.identify_residues_ravikumar(self.mixed_res_hyphen_sent)

        test_set = [match.__dict__ for match in matches]

        truth = [{'aminoacid': 'Arginine',
                  'end': 65,
                  'position': ['96'],
                  'start': 59,
                  'string': 'Arg-96',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Serine',
                  'end': 106,
                  'position': ['91'],
                  'start': 100,
                  'string': 'Ser-91',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Glycine',
                  'end': 162,
                  'position': ['97'],
                  'start': 156,
                  'string': 'Gly-97',
                  'threeletter': 'Gly'},
                 {'aminoacid': ['Ser', 'Gly'],
                  'end': 205,
                  'position': ['91', '97'],
                  'start': 192,
                  'string': 'Ser-91/Gly-97',
                  'threeletter': ['Ser', 'Gly']},
                 {'aminoacid': ['Thr', 'Ser'],
                  'end': 225,
                  'position': ['163', '169'],
                  'start': 210,
                  'string': 'Thr-163/Ser-169',
                  'threeletter': ['Thr', 'Ser']}]
        self.assertEqual(truth, test_set)

    def test_identify_residues_on_mixed_hyphen_testsent_decompose(self):
        matches = pyre.identify_residues_ravikumar(self.mixed_res_hyphen_sent)

        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]

        truth = [{'aminoacid': 'Arginine',
                  'end': 65,
                  'position': 96,
                  'residue': 'Arg96',
                  'start': 59,
                  'string': 'Arg-96',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Serine',
                  'end': 106,
                  'position': 91,
                  'residue': 'Ser91',
                  'start': 100,
                  'string': 'Ser-91',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Glycine',
                  'end': 162,
                  'position': 97,
                  'residue': 'Gly97',
                  'start': 156,
                  'string': 'Gly-97',
                  'threeletter': 'Gly'},
                 {'aminoacid': 'Serine',
                  'end': 205,
                  'position': '91',
                  'residue': 'Ser91',
                  'start': 192,
                  'string': 'Ser-91/Gly-97',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Glycine',
                  'end': 205,
                  'position': '97',
                  'residue': 'Gly97',
                  'start': 192,
                  'string': 'Ser-91/Gly-97',
                  'threeletter': 'Gly'},
                 {'aminoacid': 'Threonine',
                  'end': 225,
                  'position': '163',
                  'residue': 'Thr163',
                  'start': 210,
                  'string': 'Thr-163/Ser-169',
                  'threeletter': 'Thr'},
                 {'aminoacid': 'Serine',
                  'end': 225,
                  'position': '169',
                  'residue': 'Ser169',
                  'start': 210,
                  'string': 'Thr-163/Ser-169',
                  'threeletter': 'Ser'}]
        self.assertEqual(truth, test_set)

    def test_identify_residues_on_grammatical_testsent(self):
        matches = pyre.identify_residues_ravikumar(self.grammatical_sent)
        test_set = [match.__dict__ for match in matches]
        truth = [{'aminoacid': 'Serine',
                  'end': 344,
                  'position': ['91'],
                  'start': 323,
                  'string': 'serine at position 91',
                  'threeletter': 'Ser'}]

        self.assertEqual(truth, test_set)

    def test_identify_residues_on_grammatical_testsent_decompose(self):
        matches = pyre.identify_residues_ravikumar(self.grammatical_sent)
        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]
        truth = [{'aminoacid': 'Serine',
                  'end': 344,
                  'position': 91,
                  'residue': 'Ser91',
                  'start': 323,
                  'string': 'serine at position 91',
                  'threeletter': 'Ser'}]
        self.assertEqual(truth, test_set)

    def test_identify_residues_on_grammatical_testsent_end(self):
        matches = pyre.identify_residues_ravikumar(self.grammatical_sent_end)
        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]
        truth = [{'aminoacid': 'Serine',
                  'end': 344,
                  'position': 91,
                  'residue': 'Ser91',
                  'start': 323,
                  'string': 'serine at position 91',
                  'threeletter': 'Ser'}]
        self.assertEqual(truth, test_set)

    def test_identify_residues_on_grammatical_testsent_start(self):
        matches = pyre.identify_residues_ravikumar(self.grammatical_sent_start)
        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]
        truth = [{'aminoacid': 'Aspartate',
                  'end': 25,
                  'position': 202,
                  'residue': 'Asp202',
                  'start': 0,
                  'string': 'aspartate at position 202',
                  'threeletter': 'Asp'}]
        self.assertEqual(truth, test_set)

    def test_identify_residues_on_grammatical_testsent_multiple(self):
        matches = pyre.identify_residues_ravikumar(self.grammatical_sent_multiple)

        test_set = [match.__dict__ for match in matches]
        truth = [{'aminoacid': 'Cysteine',
                  'end': 86,
                  'position': ['275'],
                  'start': 54,
                  'string': 'cysteine residue at position 275',
                  'threeletter': 'Cys'},
                 {'aminoacid': 'Lysine',
                  'end': 301,
                  'position': ['219', '338'],
                  'start': 260,
                  'string': 'lysine residues, at positions 219 and 338',
                  'threeletter': 'Lys'}]

        self.assertEqual(truth, test_set)

    def test_identify_residues_on_grammatical_testsent_multiple_decompose(self):
        matches = pyre.identify_residues_ravikumar(self.grammatical_sent_multiple)

        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]

        truth = [{'aminoacid': 'Cysteine',
                  'end': 86,
                  'position': 275,
                  'residue': 'Cys275',
                  'start': 54,
                  'string': 'cysteine residue at position 275',
                  'threeletter': 'Cys'},
                 {'aminoacid': 'Lysine',
                  'end': 301,
                  'position': '219',
                  'residue': 'Lys219',
                  'start': 260,
                  'string': 'lysine residues, at positions 219 and 338',
                  'threeletter': 'Lys'},
                 {'aminoacid': 'Lysine',
                  'end': 301,
                  'position': '338',
                  'residue': 'Lys338',
                  'start': 260,
                  'string': 'lysine residues, at positions 219 and 338',
                  'threeletter': 'Lys'}]

        self.assertEqual(truth, test_set)

    def test_identify_residues_on_long_dashed_pos_dashed_testsent_start(self):
        """
        Need to be careful with this one - doesn't actually make 100% sense until you
        decompose it!
        """
        matches = pyre.identify_residues_ravikumar(self.long_dashed_pos_dashed_sent_start)

        test_set = [match.__dict__ for match in matches]

        truth = [{'aminoacid': ['Glu',
                                'Tyr',
                                'Trp',
                                'Ser',
                                'Arg',
                                'Gln',
                                'Asn',
                                'Thr',
                                'Ala',
                                'Asn'],
                  'end': 74,
                  'position': ['30',
                               '33',
                               '71',
                               '91',
                               '96',
                               '159',
                               '111',
                               '163',
                               '164',
                               '202'],
                  'start': 0,
                  'string': 'Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': ['Glu',
                                  'Tyr',
                                  'Trp',
                                  'Ser',
                                  'Arg',
                                  'Gln',
                                  'Asn',
                                  'Thr',
                                  'Ala',
                                  'Asn']}]
        self.assertEqual(truth, test_set)

    def test_identify_residues_on_long_dashed_pos_dashed_testsent_start_decompose(self):
        """
        Need to be careful with this one - doesn't actually make 100% sense until you
        decompose it!
        """
        matches = pyre.identify_residues_ravikumar(self.long_dashed_pos_dashed_sent_start)

        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]

        truth = [{'aminoacid': 'Glutamic acid (Glutamate)',
                  'end': 74,
                  'position': '30',
                  'residue': 'Glu30',
                  'start': 0,
                  'string': 'Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Glu'},
                 {'aminoacid': 'Tyrosine',
                  'end': 74,
                  'position': '33',
                  'residue': 'Tyr33',
                  'start': 0,
                  'string': 'Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Tyr'},
                 {'aminoacid': 'Tryptophan',
                  'end': 74,
                  'position': '71',
                  'residue': 'Trp71',
                  'start': 0,
                  'string': 'Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Trp'},
                 {'aminoacid': 'Serine',
                  'end': 74,
                  'position': '91',
                  'residue': 'Ser91',
                  'start': 0,
                  'string': 'Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Arginine',
                  'end': 74,
                  'position': '96',
                  'residue': 'Arg96',
                  'start': 0,
                  'string': 'Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Glutamine',
                  'end': 74,
                  'position': '159',
                  'residue': 'Gln159',
                  'start': 0,
                  'string': 'Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Gln'},
                 {'aminoacid': 'Asparagine',
                  'end': 74,
                  'position': '111',
                  'residue': 'Asn111',
                  'start': 0,
                  'string': 'Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Asn'},
                 {'aminoacid': 'Threonine',
                  'end': 74,
                  'position': '163',
                  'residue': 'Thr163',
                  'start': 0,
                  'string': 'Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Thr'},
                 {'aminoacid': 'Alanine',
                  'end': 74,
                  'position': '164',
                  'residue': 'Ala164',
                  'start': 0,
                  'string': 'Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Ala'},
                 {'aminoacid': 'Asparagine',
                  'end': 74,
                  'position': '202',
                  'residue': 'Asn202',
                  'start': 0,
                  'string': 'Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Asn'}]

    def test_identify_residues_on_long_dashed_pos_dashed_testsent(self):
        """
        Need to be careful with this one - doesn't actually make sense until you
        decompose it!
        """
        matches = pyre.identify_residues_ravikumar(self.long_dashed_pos_dashed_sent)

        test_set = [match.__dict__ for match in matches]
        truth = [{'aminoacid': ['Glu',
                                'Tyr',
                                'Trp',
                                'Ser',
                                'Arg',
                                'Gln',
                                'Asn',
                                'Thr',
                                'Ala',
                                'Asn'],
                  'end': 107,
                  'position': ['30',
                               '33',
                               '71',
                               '91',
                               '96',
                               '159',
                               '111',
                               '163',
                               '164',
                               '202'],
                  'start': 32,
                  'string': ' Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': ['Glu',
                                  'Tyr',
                                  'Trp',
                                  'Ser',
                                  'Arg',
                                  'Gln',
                                  'Asn',
                                  'Thr',
                                  'Ala',
                                  'Asn']}]

        self.assertEqual(truth, test_set)

    def test_identify_residues_on_long_dashed_pos_dashed_testsent_decompose(self):
        """
        Need to be careful with this one - doesn't actually make sense until you
        decompose it!
        """
        matches = pyre.identify_residues_ravikumar(self.long_dashed_pos_dashed_sent)

        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]
        truth = [{'aminoacid': 'Glutamic acid (Glutamate)',
                  'end': 107,
                  'position': '30',
                  'residue': 'Glu30',
                  'start': 32,
                  'string': ' Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Glu'},
                 {'aminoacid': 'Tyrosine',
                  'end': 107,
                  'position': '33',
                  'residue': 'Tyr33',
                  'start': 32,
                  'string': ' Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Tyr'},
                 {'aminoacid': 'Tryptophan',
                  'end': 107,
                  'position': '71',
                  'residue': 'Trp71',
                  'start': 32,
                  'string': ' Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Trp'},
                 {'aminoacid': 'Serine',
                  'end': 107,
                  'position': '91',
                  'residue': 'Ser91',
                  'start': 32,
                  'string': ' Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Arginine',
                  'end': 107,
                  'position': '96',
                  'residue': 'Arg96',
                  'start': 32,
                  'string': ' Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Glutamine',
                  'end': 107,
                  'position': '159',
                  'residue': 'Gln159',
                  'start': 32,
                  'string': ' Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Gln'},
                 {'aminoacid': 'Asparagine',
                  'end': 107,
                  'position': '111',
                  'residue': 'Asn111',
                  'start': 32,
                  'string': ' Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Asn'},
                 {'aminoacid': 'Threonine',
                  'end': 107,
                  'position': '163',
                  'residue': 'Thr163',
                  'start': 32,
                  'string': ' Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Thr'},
                 {'aminoacid': 'Alanine',
                  'end': 107,
                  'position': '164',
                  'residue': 'Ala164',
                  'start': 32,
                  'string': ' Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Ala'},
                 {'aminoacid': 'Asparagine',
                  'end': 107,
                  'position': '202',
                  'residue': 'Asn202',
                  'start': 32,
                  'string': ' Glu-30-Tyr-33-Trp-71-Ser-91-Arg-96-Gln-159-Asn-111-Thr-163-Ala-164-Asn-202',
                  'threeletter': 'Asn'}]

        self.assertEqual(truth, test_set)

    def test_identify_residues_on_long_dashed_testsent(self):
        matches = pyre.identify_residues_ravikumar(self.long_dashed_sent)

        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]
        truth = [{'aminoacid': ['Glu',
                                'Tyr',
                                'Trp',
                                'Ser',
                                'Arg',
                                'Gln',
                                'Asn',
                                'Thr',
                                'Ala',
                                'Asn'],
                  'end': 94,
                  'position': ['30',
                               '33',
                               '71',
                               '91',
                               '96',
                               '159',
                               '111',
                               '163',
                               '164',
                               '202'],
                  'start': 29,
                  'string': ' Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202',
                  'threeletter': ['Glu',
                                  'Tyr',
                                  'Trp',
                                  'Ser',
                                  'Arg',
                                  'Gln',
                                  'Asn',
                                  'Thr',
                                  'Ala',
                                  'Asn']}]

    def test_identify_residues_on_long_dashed_testsent_decompose(self):
        matches = pyre.identify_residues_ravikumar(self.long_dashed_sent)

        test_set = [match.__dict__ for match in pyre.decompose_matches(matches)]
        truth = [{'aminoacid': 'Glutamic acid (Glutamate)',
                  'end': 94,
                  'position': '30',
                  'residue': 'Glu30',
                  'start': 29,
                  'string': ' Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202',
                  'threeletter': 'Glu'},
                 {'aminoacid': 'Tyrosine',
                  'end': 94,
                  'position': '33',
                  'residue': 'Tyr33',
                  'start': 29,
                  'string': ' Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202',
                  'threeletter': 'Tyr'},
                 {'aminoacid': 'Tryptophan',
                  'end': 94,
                  'position': '71',
                  'residue': 'Trp71',
                  'start': 29,
                  'string': ' Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202',
                  'threeletter': 'Trp'},
                 {'aminoacid': 'Serine',
                  'end': 94,
                  'position': '91',
                  'residue': 'Ser91',
                  'start': 29,
                  'string': ' Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202',
                  'threeletter': 'Ser'},
                 {'aminoacid': 'Arginine',
                  'end': 94,
                  'position': '96',
                  'residue': 'Arg96',
                  'start': 29,
                  'string': ' Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Glutamine',
                  'end': 94,
                  'position': '159',
                  'residue': 'Gln159',
                  'start': 29,
                  'string': ' Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202',
                  'threeletter': 'Gln'},
                 {'aminoacid': 'Asparagine',
                  'end': 94,
                  'position': '111',
                  'residue': 'Asn111',
                  'start': 29,
                  'string': ' Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202',
                  'threeletter': 'Asn'},
                 {'aminoacid': 'Threonine',
                  'end': 94,
                  'position': '163',
                  'residue': 'Thr163',
                  'start': 29,
                  'string': ' Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202',
                  'threeletter': 'Thr'},
                 {'aminoacid': 'Alanine',
                  'end': 94,
                  'position': '164',
                  'residue': 'Ala164',
                  'start': 29,
                  'string': ' Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202',
                  'threeletter': 'Ala'},
                 {'aminoacid': 'Asparagine',
                  'end': 94,
                  'position': '202',
                  'residue': 'Asn202',
                  'start': 29,
                  'string': ' Glu30-Tyr33-Trp71-Ser91-Arg96-Gln159-Asn111-Thr163-Ala164-Asn202',
                  'threeletter': 'Asn'}]

    ## Old Tests
    def _test_decompose_slashed_residue(self):
        test_string =  "Arg99/101"
        expected = ['Arg99', 'Arg101']
        decomposed = pyre.decompose_slashed_residue(test_string)
        self.assertEqual(sorted(decomposed), sorted(expected))

    def _test_subdivide_compound_residues_good(self):
        good = "Arg101/103"
        expected = ["Arg101", "Arg103"]
        subdivision = pyre.subdivide_compound_residues(good)
        self.assertEqual(sorted(subdivision), sorted(expected))

    def _test_subdivide_compound_residues_bad(self):
        bad = "Thr163/Ser164"
        expected = ["Thr163", "Ser164"]
        subdivision = pyre.subdivide_compound_residues(bad)
        self.assertEqual(sorted(subdivision), sorted(expected))

    def _test_subdivide_compound_residues_ugly(self):
        ugly = "Thr163/Ser164/166"
        expected = ["Thr163", "Ser164", "Ser166"]
        subdivision = pyre.subdivide_compound_residues(ugly)
        self.assertEqual(sorted(subdivision), sorted(expected))

    def _test_subdivide_compound_residues_a(self):
        A = "Tyr33"
        expected = ["Tyr33"]
        subdivision = pyre.subdivide_compound_residues(A)
        self.assertEqual(sorted(subdivision), sorted(expected))

    def _test_subdivide_compound_residues_few_dollars(self):
        few_dollars = "Thr163/164/Ser166"
        expected = ["Thr163", "Thr164", "Ser166"]
        subdivision = pyre.subdivide_compound_residues(few_dollars)
        self.assertEqual(sorted(subdivision), sorted(expected))

    def _test_subdivide_compound_residues_more(self):
        more="Thr163/164/Ser166/167"
        expected = ["Thr163", "Thr164", "Ser166", "Ser167"]
        subdivision = pyre.subdivide_compound_residues(more)
        self.assertEqual(sorted(subdivision), sorted(expected))

    ## New Tests
    def test_decompose_slashed_residue(self):
        test_string = "Arg99/101"

        matches = pyre.identify_residues(test_string)
        decomposed = pyre.decompose_matches(matches)
        match_list = [match.__dict__ for match in decomposed]

        truth = [{'aminoacid': 'Arginine',
                  'end': 9,
                  'position': '99',
                  'residue': 'Arg99',
                  'start': 0,
                  'string': 'Arg99/101',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Arginine',
                  'end': 9,
                  'position': '101',
                  'residue': 'Arg101',
                  'start': 0,
                  'string': 'Arg99/101',
                  'threeletter': 'Arg'}]

        self.assertEqual(match_list, truth)

    def test_decompose_compound_residues_good(self):
        good = "Arg101/103"
        matches = pyre.identify_residues(good)
        decomposed = pyre.decompose_matches(matches)
        match_list = [match.__dict__ for match in decomposed]

        truth = [{'aminoacid': 'Arginine',
                  'end': 10,
                  'position': '101',
                  'residue': 'Arg101',
                  'start': 0,
                  'string': 'Arg101/103',
                  'threeletter': 'Arg'},
                 {'aminoacid': 'Arginine',
                  'end': 10,
                  'position': '103',
                  'residue': 'Arg103',
                  'start': 0,
                  'string': 'Arg101/103',
                  'threeletter': 'Arg'}]

        self.assertEqual(match_list, truth)


    def test_decompose_compound_residues_bad(self):
        bad = "Thr163/Ser164"
        matches = pyre.identify_residues(bad)
        decomposed = pyre.decompose_matches(matches)
        match_list = [match.__dict__ for match in decomposed]

        truth = [{'aminoacid': 'Threonine',
                  'end': 13,
                  'position': '163',
                  'residue': 'Thr163',
                  'start': 0,
                  'string': 'Thr163/Ser164',
                  'threeletter': 'Thr'},
                 {'aminoacid': 'Serine',
                  'end': 13,
                  'position': '164',
                  'residue': 'Ser164',
                  'start': 0,
                  'string': 'Thr163/Ser164',
                  'threeletter': 'Ser'}]

        self.assertEqual(match_list, truth)


    ## TODO - possibly use more inclusive regex followed by "subdivide_compound_residues" within the matchclass
    def _test_subdivide_compound_residues_ugly(self):
        ugly = "Thr163/Ser164/166"
        expected = ["Thr163", "Ser164", "Ser166"]
        subdivision = pyre.subdivide_compound_residues(ugly)
        self.assertEqual(sorted(subdivision), sorted(expected))

    def _test_subdivide_compound_residues_a(self):
        A = "Tyr33"
        expected = ["Tyr33"]
        subdivision = pyre.subdivide_compound_residues(A)
        self.assertEqual(sorted(subdivision), sorted(expected))

    def _test_subdivide_compound_residues_few_dollars(self):
        few_dollars = "Thr163/164/Ser166"
        expected = ["Thr163", "Thr164", "Ser166"]
        subdivision = pyre.subdivide_compound_residues(few_dollars)
        self.assertEqual(sorted(subdivision), sorted(expected))

    def _test_subdivide_compound_residues_more(self):
        more = "Thr163/164/Ser166/167"
        expected = ["Thr163", "Thr164", "Ser166", "Ser167"]
        subdivision = pyre.subdivide_compound_residues(more)
        self.assertEqual(sorted(subdivision), sorted(expected))


    ## Residue locator - these are now redundant
    res_locations_dict = OrderedDict([('Tyr33',
                                      OrderedDict([('locations', [16]),
                                                   ('string', ['Tyr33/39']),
                                                   ('freq', 1)])),
                                     ('Tyr39',
                                      OrderedDict([('locations', [16]),
                                                   ('string', ['Tyr33/39']),
                                                   ('freq', 1)])),
                                     ('Trp71',
                                      OrderedDict([('locations', [17]),
                                                   ('string', ['andTrp71/77']),
                                                   ('freq', 1)])),
                                     ('Trp77',
                                      OrderedDict([('locations', [17]),
                                                   ('string', ['andTrp71/77']),
                                                   ('freq', 1)])),
                                     ('Ala88',
                                      OrderedDict([('locations', [49]),
                                                   ('string', ['Ala88/94']),
                                                   ('freq', 1)])),
                                     ('Ala94',
                                      OrderedDict([('locations', [49]),
                                                   ('string', ['Ala88/94']),
                                                   ('freq', 1)])),
                                     ('Glu30',
                                      OrderedDict([('locations', [42]),
                                                   ('string', ['Glu30/36']),
                                                   ('freq', 1)])),
                                     ('Glu36',
                                      OrderedDict([('locations', [42]),
                                                   ('string', ['Glu30/36']),
                                                   ('freq', 1)])),
                                     ('Gln159',
                                      OrderedDict([('locations', [44]),
                                                   ('string', ['Gln159/165']),
                                                   ('freq', 1)])),
                                     ('Gln165',
                                      OrderedDict([('locations', [44]),
                                                   ('string', ['Gln159/165']),
                                                   ('freq', 1)]))])

    def _test_locate_residues_sample_sent(self):
        unique_residues = pyre._identify_residues(self.slashed_res_sent, return_mentions=True)
        sentdoc = nlp(self.slashed_res_sent)
        sent_tkns = [t.text for t in sentdoc]

        location_dict = pyre._locate_residues(sent_tkns, unique_residues, verbose=False)

        sorted(location_dict.keys())

        # self.assertEqual(location_dict, self.res_locations_dict)
        self.assertEqual( sorted(self.res_locations_dict), sorted(location_dict.keys()))
        self.assertDictEqual(deep_sort(location_dict), deep_sort(self.res_locations_dict))

    def test_locate_residues(self):
        pmc_id = "PMC5740067"

        pyre.process([pmc_id, ], this_dir, filename="test_validate_locations_dict_PMC5740067.json", provider="West-Life",
                                 cifscantype='flex', save=True, overwrite=True,
                                 return_dict=False, decompose=True, verbose=True)

        infile = os.path.join(this_dir, "test_validate_locations_dict_PMC5740067.json")
        with open(infile, "r") as ifile:
            # expected_dict = json.load(ifile)
            test_dict_as_string = ifile.readlines()[0]

        infile = os.path.join(this_dir, "test_locations_dict_PMC5740067.json")
        with open(infile, "r") as ifile:
            # expected_dict = json.load(ifile)
            expected_dict_as_string = ifile.readlines()[0]

        self.maxDiff = None
        self.assertEqual(test_dict_as_string, expected_dict_as_string)

    ## Protein identification pyre.identify_protein_ID

    prot_sent = "Coordinates and structure factors have been deposited at the Protein Data Bank (PDB) under accession " + \
                "codes 5ORE for OccJ structure, 5ORG for OccJ-octopine structure, 5OT8 for NocT-G97S-octopine " + \
                "structure, 5OTA for NocT-octopinic acid structure, 5OTC for NocT-noroctopinic acid, 5OT9 for " + \
                "`NocT-histopine structure and 5OVZ for the high resolution structure of NocT-nopaline complex."

    prot_sent_invalid = "0000 " + prot_sent

    def test_identify_protein_ID_simple(self):
        expected = sorted(['5OTC', '5OT9', '5ORE', '5OTA', '5OVZ', '5OT8', '5ORG', '0000'])

        unique_proteins = sorted(pyre.identify_protein_ID(self.prot_sent_invalid, simple=True))
        self.assertEqual(expected, unique_proteins)

    def test_identify_protein_ID(self):
        expected = sorted(['5OTC', '5OT9', '5ORE', '5OTA', '5OVZ', '5OT8', '5ORG'])

        unique_proteins = sorted(pyre.identify_protein_ID(self.prot_sent_invalid, simple=False))
        self.assertEqual(expected, unique_proteins)

    def test_locate_proteins_PMC5740067(self):
        ext_id = "PMC5740067"

        source = pyre.SourceClass
        source.text_dict = pyre.get_text(ext_id)
        source.fulltext = pyre.reconstruct_fulltext(source.text_dict, tokenise=False)
        protein_locations = pyre.locate_proteins(source.fulltext)

        protein_locations_dictlist = [loc.__dict__ for loc in protein_locations]

        expected = [{'start': 12827, 'end': 12831, 'string': '5ITP'},
                    {'start': 15400, 'end': 15404, 'string': '4POX'},
                    {'start': 17161, 'end': 17165, 'string': '5ITP'},
                    {'start': 17514, 'end': 17518, 'string': '5ITO'},
                    {'start': 18814, 'end': 18818, 'string': '1LAF'},
                    {'start': 20200, 'end': 20204, 'string': '1LAF'},
                    {'start': 39503, 'end': 39507, 'string': '4POW'},
                    {'start': 39630, 'end': 39634, 'string': '5ITP'},
                    {'start': 42443, 'end': 42447, 'string': '5ORE'},
                    {'start': 42468, 'end': 42472, 'string': '5ORG'},
                    {'start': 42502, 'end': 42506, 'string': '5OT8'},
                    {'start': 42541, 'end': 42545, 'string': '5OTA'},
                    {'start': 42581, 'end': 42585, 'string': '5OTC'},
                    {'start': 42614, 'end': 42618, 'string': '5OT9'},
                    {'start': 42652, 'end': 42656, 'string': '5OVZ'}]

        self.assertEqual(expected, protein_locations_dictlist)

    def test_locate_proteins_sent(self):

        protein_locations = pyre.locate_proteins(self.prot_sent)
        protein_locations_dictlist = [loc.__dict__ for loc in protein_locations]

        expected = [{'end': 111, 'start': 107, 'string': '5ORE'},
                    {'end': 136, 'start': 132, 'string': '5ORG'},
                    {'end': 170, 'start': 166, 'string': '5OT8'},
                    {'end': 209, 'start': 205, 'string': '5OTA'},
                    {'end': 249, 'start': 245, 'string': '5OTC'},
                    {'end': 282, 'start': 278, 'string': '5OT9'},
                    {'end': 321, 'start': 317, 'string': '5OVZ'}]

        self.assertEqual(expected, protein_locations_dictlist)

    def test_locate_proteins_sent_invalid(self):

        protein_locations = pyre.locate_proteins(self.prot_sent_invalid)
        protein_locations_dictlist = [loc.__dict__ for loc in protein_locations]

        expected = [{'end': 116, 'start': 112, 'string': '5ORE'},
                    {'end': 141, 'start': 137, 'string': '5ORG'},
                    {'end': 175, 'start': 171, 'string': '5OT8'},
                    {'end': 214, 'start': 210, 'string': '5OTA'},
                    {'end': 254, 'start': 250, 'string': '5OTC'},
                    {'end': 287, 'start': 283, 'string': '5OT9'},
                    {'end': 326, 'start': 322, 'string': '5OVZ'}]

        self.assertEqual(expected, protein_locations_dictlist)


    ## API queries

    def test_ePMC_id_query(self):
        query = "((IN_EPMC:y) AND (HAS_PDB:y)) AND OPEN_ACCESS:Y AND FIRST_PDATE:[1990 TO 2000]"
        pmc_id_list = pyre.extract_ids_from_query(query=query)

        truth = ['PMC2175196',
                 'PMC2792768',
                 'PMC2193205',
                 'PMC2174853',
                 'PMC2174850',
                 'PMC2193177',
                 'PMC2151181',
                 'PMC2195553',
                 'PMC2746121',
                 'PMC2192988',
                 'PMC2132952',
                 'PMC2213406',
                 'PMC2212146',
                 'PMC2196387',
                 'PMC2120550',
                 'PMC2119855',
                 'PMC2289279']

        self.assertEqual(sorted(pmc_id_list), sorted(truth))

    def test_PMCID_to_PMID_conversion(self):
        pmc_id = "PMC5740067"
        expected = "29269740"
        pm_id = pyre.convert_PMCID_to_PMID(ext_id=pmc_id)
        self.assertTrue(pm_id, expected)

    def test_request_fulltextXML_request_success(self):
        pmc_id = "PMC5740067"

        r = pyre.request_fulltextXML(pmc_id)
        self.assertEqual(r.status_code, 200)

    def test_invalid_request_fulltextXML_returns_404(self):
        pmc_id = "PMC0000"

        r = pyre.request_fulltextXML(pmc_id)
        # self.assertIsNone(r)
        self.assertEqual(404, r.status_code)

    def test_request_fulltextXML_returns_XML(self):
        pmc_id = "PMC5740067"

        r = pyre.request_fulltextXML(pmc_id)
        self.assertTrue(isinstance(r.text, str))

    ## mmCIF

    def test_mmCIF_dir_exists(self):
        self.assertTrue(os.path.isdir(pyre.mmCIF_dir))

    ## PDB

    def test_PDB_dir_exists(self):
        self.assertTrue(os.path.isdir(pyre.PDB_dir))


if __name__ == "__main__":
    # test = False
    test = True

    if test:

        print("Running test suite:\n")

        suite = unittest.TestLoader().loadTestsFromTestCase(TestClass)
        unittest.TextTestRunner(verbosity=2).run(suite)
else:
    pass
