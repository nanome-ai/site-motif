from unittest import TestCase
import os
import subprocess
import tempfile


class MultipleSiteAlignmentTestCase(TestCase):

    # def setUp(self):
    #     # self.temp_dir = tempfile.TemporaryDirectory().name
    #     # Delete intermediate files if they exist
    #     if os.path.isfile(self.pairs_list_file):
    #         os.remove(self.pairs_list_file)
    #     if os.path.isfile(self.pdb_size_file):
    #         os.remove(self.pdb_size_file)
    #     if os.path.isfile(self.align_output_file):
    #         os.remove(self.align_output_file)

    # def tearDown(self):
    #     self.temp_dir.cleanup()

    def test_atp_subset_alignment(self):
        with tempfile.TemporaryDirectory() as output_dir:
            pairs_list_file = f'{output_dir}/PairList.txt'
            pdb_size_file = f'{output_dir}/PDBSize.txt'
            align_output_file = f'{output_dir}/align_output.txt'
            sites_folder = 'test_data/ATP_SUBSET'

            cmd = f'./bin/site-motif {sites_folder} {output_dir}'
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
