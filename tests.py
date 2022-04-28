from unittest import TestCase
import os


class MultipleSiteAlignmentTestCase(TestCase):

    def setUp(self):
        self.pairs_list_file = 'PairList.txt'
        self.pdb_size_file = 'PDBSize.txt'
        self.align_output_file = 'align_output.txt'
        
        # Delete intermediate files if they exist
        if os.path.isfile(self.pairs_list_file):
            os.remove(self.pairs_list_file)
        if os.path.isfile(self.pdb_size_file):
            os.remove(self.pdb_size_file)
        if os.path.isfile(self.align_output_file):
            os.remove(self.align_output_file)
            

    def test_alignment(self):
        self.assertFalse(os.path.isfile(self.pairs_list_file))
        self.assertFalse(os.path.isfile(self.pdb_size_file))       
        self.assertFalse(os.path.isfile(self.align_output_file))

        test_folder = 'test_data/ATP2'
        cmd = f'./run.sh {test_folder}'
        os.system(cmd)
        
        self.assertTrue(os.path.isfile(self.pairs_list_file))
        self.assertTrue(os.path.isfile(self.pdb_size_file))       
        self.assertTrue(os.path.isfile(self.align_output_file))
        with open(self.pairs_list_file, 'r') as f:
            self.assertNotEqual(f.read(), '')
        with open(self.pdb_size_file, 'r') as f:
            self.assertNotEqual(f.read(), '')
        with open(self.align_output_file, 'r') as f:
            self.assertNotEqual(f.read(), '')
