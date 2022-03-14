import unittest, importlib
from tempfile import TemporaryDirectory
import io
from contextlib import redirect_stdout
from intersect_callsets import Variant, record_from_cluster, reciprocal_overlap, merge_variants, find_clusters, parse_vcf, coordinates_match

class Test_Variant(unittest.TestCase):

    def test_variant(self):
        variant1 = Variant('chr1', 0, 3, 2, 'DEL', 'test1', 10)
        variant2 = Variant('chr1', 0, 3, 2, 'DEL', 'test1', 10)
        variant3 = Variant('chr1', 0, 3, 2, 'DEL', 'test2', 30)
        self.assertEqual(variant1, variant2)
        self.assertFalse(variant1 == variant3)
        self.assertEqual(len(variant1), 3)

    def test_variant2(self):
        variant1 = Variant('chr1', 1, 3, 2, 'DEL', 'test1', 10)
        variant2 = Variant('chr1', 2, 3, 1, 'DEL', 'test1', 10)
        variant3 = Variant('chr2', 1, 3, 2, 'DEL', 'test1', 10)
        self.assertFalse(variant1 == variant3)
        self.assertTrue(variant1 < variant3)
        self.assertTrue(variant1 < variant2)
        self.assertFalse(variant2 < variant2)


class Test_record_from_cluster(unittest.TestCase):

    def test_record_from_cluster(self):
        cluster = [
         Variant('chr1', 0, 3, 3, 'DEL', 'test1', 10), Variant('chr1', 0, 2, 2, 'DEL', 'test1', 20), Variant('chr1', 1, 4, 3, 'DEL', 'test2', 30)]
        callsets = ['test2', 'test1']
        record = record_from_cluster(cluster, callsets)[0]
        expected = ['chr1-1-DEL-3-test2', 'chr1', '1', '4', '3', 'DEL', '30', 'nan', 'nan', 'True', 'True', 'chr1-1-DEL-3-test2', 'chr1-0-DEL-3-test1;chr1-0-DEL-2-test1']
        self.assertEqual(record, expected)


class Test_reciprocal_overlap(unittest.TestCase):

    def test_reciprocal_overlap(self):
        variant1 = Variant('chr1', 0, 3, 2, 'DEL', 'test1', 10)
        variant2 = Variant('chr1', 0, 5, 2, 'DEL', 'test1', 20)
        variant3 = Variant('chr1', 0, 10, 10, 'DEL', 'test1', 10)
        self.assertTrue(reciprocal_overlap(variant1, variant2, 0.5))
        self.assertTrue(reciprocal_overlap(variant2, variant1, 0.5))
        self.assertFalse(reciprocal_overlap(variant1, variant3, 0.5))
        self.assertFalse(reciprocal_overlap(variant3, variant1, 0.5))


class Test_coordinates_match(unittest.TestCase):

    def test_coordinates_match(self):
        variant1 = Variant('chr1', 100, 200, 100, 'DEL', 'test1', 10)
        variant2 = Variant('chr1', 150, 250, 100, 'DEL', 'test1', 10)
        variant3 = Variant('chr1', 210, 250, 40, 'DEL', 'test1', 10)
        variant4 = Variant('chr1', 210, 310, 100, 'DEL', 'test1', 10)
        self.assertTrue(coordinates_match(variant1, variant2))
        self.assertFalse(coordinates_match(variant1, variant3))
        self.assertFalse(coordinates_match(variant2, variant3))
        self.assertTrue(coordinates_match(variant2, variant4))


class Test_merge_variants(unittest.TestCase):

    def test_merge_variants1(self):
        variant1 = Variant('chr1', 0, 3, 3, 'DEL', 'test1', 20)
        variant2 = Variant('chr1', 0, 5, 5, 'DEL', 'test2', 10)
        variant3 = Variant('chr1', 0, 3, 3, 'DEL', 'test3', 20)
        expected = [
         ['chr1-0-DEL-3-test1', 'chr1', '0', '3', '3', 'DEL', '20', 'nan', 'nan','True', 'True', 'True', 'chr1-0-DEL-3-test1', 'chr1-0-DEL-5-test2', 'chr1-0-DEL-3-test3']]
        self.assertEqual(merge_variants([variant1, variant2, variant3], ['test1', 'test2', 'test3'])[0], expected)
        expected = [
         ['chr1-0-DEL-5-test2', 'chr1', '0', '5', '5', 'DEL', '10', 'nan', 'nan','False', 'True', 'True', 'nan', 'chr1-0-DEL-5-test2', 'chr1-0-DEL-3-test3']]
        self.assertEqual(merge_variants([variant2, variant3], ['test1', 'test2', 'test3'])[0], expected)
        expected = [
         ['chr1-0-DEL-5-test2', 'chr1', '0', '5', '5', 'DEL', '10', 'nan', 'nan', 'False', 'True', 'False', 'nan', 'chr1-0-DEL-5-test2', 'nan']]
        self.assertEqual(merge_variants([variant2], ['test1', 'test2', 'test3'])[0], expected)

    def test_merge_variants2(self):
        variant1 = Variant('chr1', 100, 300, 200, 'DEL', 'test1', 20)
        variant2 = Variant('chr1', 400, 500, 100, 'DEL', 'test2', 10)
        variant3 = Variant('chr1', 100, 130, 30, 'DEL', 'test3', 30)
        expected = [
         ['chr1-100-DEL-200-test1', 'chr1', '100', '300', '200', 'DEL', '20', 'nan', 'nan', 'True', 'False', 'False', 'chr1-100-DEL-200-test1', 'nan', 'nan'],
         ['chr1-400-DEL-100-test2', 'chr1', '400', '500', '100', 'DEL', '10', 'nan', 'nan', 'False', 'True', 'False', 'nan', 'chr1-400-DEL-100-test2', 'nan'],
         ['chr1-100-DEL-30-test3', 'chr1', '100', '130', '30', 'DEL', '30', 'nan', 'nan', 'False', 'False', 'True', 'nan', 'nan', 'chr1-100-DEL-30-test3']
         ]
        self.assertEqual(merge_variants([variant1, variant2, variant3], ['test1', 'test2', 'test3'])[0], expected)

    def test_merge_variants3(self):
        variant1 = Variant('chr1', 100, 300, 200, 'DEL', 'test1', 10, 90, 100)
        variant2 = Variant('chr1', 400, 500, 100, 'DEL', 'test1', 5, 80, 50)
        variant3 = Variant('chr1', 100, 130, 30, 'DEL', 'test1', 20, 40, 50)
        expected = [
         ['chr1-100-DEL-200-test1', 'chr1', '100', '300', '200', 'DEL', '10', '90', '100', 'True', 'False', 'False', 'chr1-100-DEL-200-test1', 'nan', 'nan'],
         ['chr1-400-DEL-100-test1', 'chr1', '400', '500', '100', 'DEL', '5', '80', '50','True', 'False', 'False', 'chr1-400-DEL-100-test1', 'nan', 'nan'],
         ['chr1-100-DEL-30-test1', 'chr1', '100', '130', '30', 'DEL', '20', '40', '50', 'True', 'False', 'False', 'chr1-100-DEL-30-test1', 'nan', 'nan']
         ]
        self.assertEqual(merge_variants([variant1, variant2, variant3], ['test1', 'test2', 'test3'])[0], expected)

    def test_merge_variants4(self):
        variant1 = Variant('chr1', 0, 5, 5, 'DEL', 'test1', 30, 5, 10)
        variant2 = Variant('chr1', 2, 5, 3, 'DEL', 'test1', 20, 4, 30)
        variant3 = Variant('chr1', 1, 6, 5, 'DEL', 'test1', 10, 3, 10)
        expected = [
         ['chr1-0-DEL-5-test1', 'chr1', '0', '5', '5', 'DEL', '30', '5', '10', 'True', 'False', 'False', 'chr1-0-DEL-5-test1;chr1-2-DEL-3-test1;chr1-1-DEL-5-test1', 'nan', 'nan']
         ]
        self.assertEqual(merge_variants([variant1, variant2, variant3], ['test1', 'test2', 'test3'])[0], expected)

    def test_merge_variants5(self):
        variant1 = Variant('chr1', 100, 300, 200, 'DEL', 'test1', 10)
        variant2 = Variant('chr1', 100, 300, 200, 'INV', 'test1', 5)
        variant3 = Variant('chr1', 100, 300, 200, 'DUP', 'test1', 20)
        # variants of different types should not be merged
        expected = [
         ['chr1-100-DEL-200-test1', 'chr1', '100', '300', '200', 'DEL', '10', 'nan', 'nan', 'True', 'False', 'False', 'chr1-100-DEL-200-test1', 'nan', 'nan'],
         ['chr1-100-INV-200-test1', 'chr1', '100', '300', '200', 'INV', '5', 'nan', 'nan', 'True', 'False', 'False', 'chr1-100-INV-200-test1', 'nan', 'nan'],
         ['chr1-100-DUP-200-test1', 'chr1', '100', '300', '200', 'DUP', '20', 'nan', 'nan', 'True', 'False', 'False', 'chr1-100-DUP-200-test1', 'nan', 'nan']
         ]
        self.assertEqual(merge_variants([variant1, variant2, variant3], ['test1', 'test2', 'test3'])[0], expected)

    def test_merge_variants6(self):
        variant1 = Variant('chr1', 0, 5, 5, 'DEL', 'test1', 30)
        variant2 = Variant('chr1', 2, 5, 3, 'DEL', 'test1', 20)
        variant3 = Variant('chr1', 1, 6, 5, 'DEL', 'test1', 10)
        expected = [
         ['chr1-0-DEL-5-test1', 'chr1', '0', '5', '5', 'DEL', '30' , 'nan', 'nan','True', 'False', 'False', 'chr1-0-DEL-5-test1', 'nan', 'nan'],
         ['chr1-2-DEL-3-test1', 'chr1', '2', '5', '3', 'DEL', '20', 'nan', 'nan', 'True', 'False', 'False', 'chr1-2-DEL-3-test1', 'nan', 'nan'],
         ['chr1-1-DEL-5-test1', 'chr1', '1', '6', '5', 'DEL', '10', 'nan', 'nan', 'True', 'False', 'False', 'chr1-1-DEL-5-test1', 'nan', 'nan']
         ]
        self.assertEqual(merge_variants([variant1, variant2, variant3], ['test1', 'test2', 'test3'], True)[0], expected)


class Test_find_clusters(unittest.TestCase):

    def test_find_clusters(self):
        variant1 = Variant('chr1', 0, 3, 3, 'DEL', 'test1', 30)
        variant2 = Variant('chr1', 0, 5, 5, 'DEL', 'test2', 20)
        variant3 = Variant('chr1', 4, 7, 3, 'DEL', 'test3', 10)
        variants = [variant1, variant2, variant3]
        clusters = [c for c in find_clusters(variants)]
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0], variants)
        variants = [variant1, variant3]
        clusters = [c for c in find_clusters(variants)]
        self.assertEqual(len(clusters), 2)
        self.assertEqual(clusters[0], [variant1])
        self.assertEqual(clusters[1], [variant3])

    def test_find_clusters2(self):
        variant1 = Variant('chr1', 0, 3, 3, 'DEL', 'test1', 30)
        variant2 = Variant('chr1', 0, 5, 5, 'DEL', 'test2', 20)
        variant3 = Variant('chr2', 4, 7, 3, 'DEL', 'test3', 10)
        variants = [variant1, variant2, variant3]
        clusters = [c for c in find_clusters(variants)]
        self.assertEqual(len(clusters), 2)
        self.assertEqual(clusters[0], [variant1, variant2])
        self.assertEqual(clusters[1], [variant3])

    def test_find_clusters3(self):
        variant1 = Variant('chr1', 0, 3, 3, 'DEL', 'test1', 30)
        variant2 = Variant('chr1', 5, 7, 2, 'DEL', 'test2', 20)
        variants = [variant1, variant2]
        clusters = [c for c in find_clusters(variants)]
        self.assertEqual(len(clusters), 2)
        self.assertEqual(clusters[0], [variant1])
        self.assertEqual(clusters[1], [variant2])

        clusters = [c for c in find_clusters(variants, 3)]
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0], [variant1, variant2])


class Test_parse_vcf(unittest.TestCase):

    def test_parse_vcf(self):
        variants = []
        vartypes = set([])
        nr_variants = parse_vcf('test.vcf', 'test2', variants, vartypes)
        self.assertEqual(nr_variants, 2)
        expected = [Variant('21', 21492142, 21492648, 506, 'DEL', 'test2', '.', '90', '40'), Variant('21', 21492145, 21492648, 506, 'DEL', 'test2', '.', '40', '30')]
        self.assertEqual(expected, variants)
