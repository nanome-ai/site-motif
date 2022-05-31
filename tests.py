from unittest import TestCase
import os
import tempfile


class MultipleSiteAlignmentTestCase(TestCase):

    def test_atp_subset_alignment(self):
        dirname = os.path.dirname(__file__)
        with tempfile.TemporaryDirectory() as output_dir:
            pairs_list_file = f'{output_dir}/PairList.txt'
            pdb_size_file = f'{output_dir}/PDBSize.txt'
            align_output_file = f'{output_dir}/align_output.txt'
            sites_folder = f'{dirname}/test_data/ATP_SUBSET'

            cmd = f'{dirname}/bin/site-motif {sites_folder} {output_dir}'
            output = os.system(cmd)
            self.assertEqual(output, 0)
            self.assertTrue(os.path.isfile(pairs_list_file))
            self.assertTrue(os.path.isfile(pdb_size_file))       
            self.assertTrue(os.path.isfile(align_output_file))
            with open(pairs_list_file, 'r') as f:
                self.assertNotEqual(f.read(), '')
            with open(pdb_size_file, 'r') as f:
                self.assertNotEqual(f.read(), '')
            with open(align_output_file, 'r') as f:
                self.assertNotEqual(f.read(), '')
