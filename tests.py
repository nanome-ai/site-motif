from unittest import TestCase
import os
import subprocess
import tempfile
from site_motif import write_pairs, write_pdb_size


class PairListTestCase(TestCase):

    def test_write_pairs(self):
        dirname = os.path.dirname(__file__)
        sites_dir = f'{dirname}/test_data/ATP'
        with tempfile.TemporaryDirectory() as output_dir:
            output_file = write_pairs(sites_dir, output_dir)
            self.assertTrue(os.path.exists(output_file))
            with open(output_file, 'r') as f:
                data = f.read()
            site_count = len(os.listdir(sites_dir))
            # Get pair counts and remove blank lines
            pair_count = len([pair for pair in data.split('\n') if pair])
            self.assertEqual(pair_count, site_count ** 2)

    def test_write_pairs_command_no_reference(self):
        """Test that the command line version of write_pairs works."""
        dirname = os.path.dirname(__file__)
        sites_dir = f'{dirname}/test_data/ATP'
        with tempfile.TemporaryDirectory() as output_dir:
            subprocess.run(['python', 'site_motif/Pairs.py', sites_dir, output_dir])
            output_file = os.path.join(output_dir, 'PairList.txt')
            self.assertTrue(os.path.exists(output_file))
            with open(output_file, 'r') as f:
                data = f.read()
            site_count = len(os.listdir(sites_dir))
            # Get pair counts and remove blank lines
            pair_count = len([pair for pair in data.split('\n') if pair])
            self.assertEqual(pair_count, site_count ** 2)

    def test_write_pairs_command_with_reference_pdb(self):
        """Test that the command line version of write_pairs works when providing reference pdb."""
        dirname = os.path.dirname(__file__)
        sites_dir = f'{dirname}/test_data/ATP'
        reference_pdb = os.listdir(sites_dir)[0]
        with tempfile.TemporaryDirectory() as output_dir:
            subprocess.run(['python', 'site_motif/Pairs.py', sites_dir, output_dir, reference_pdb])
            output_file = os.path.join(output_dir, 'PairList.txt')
            self.assertTrue(os.path.exists(output_file))
            with open(output_file, 'r') as f:
                data = f.read()
            site_count = len(os.listdir(sites_dir))
            # Get pair counts and remove blank lines
            pair_count = len([pair for pair in data.split('\n') if pair])
            self.assertEqual(pair_count, (site_count * 2) - 1)


class PDBSizeTestCase(TestCase):

    def test_write_pdb_size(self):
        """Test write_pdb_size function."""
        dirname = os.path.dirname(__file__)
        sites_dir = f'{dirname}/test_data/ATP'
        with tempfile.TemporaryDirectory() as output_dir:
            output_file = write_pdb_size(sites_dir, output_dir)
            self.assertTrue(os.path.exists(output_file))
            with open(output_file, 'r') as f:
                output_lines = f.readlines()
                output_data = [line.strip() for line in output_lines if line]
            self.assertEqual(len(output_data), len(os.listdir(sites_dir)))
    
    def test_write_pdb_size_command(self):
        """Test that the command-line version of write_pdb_size works."""
        dirname = os.path.dirname(__file__)
        sites_dir = f'{dirname}/test_data/ATP'
        with tempfile.TemporaryDirectory() as output_dir:
            subprocess.run(['python', 'site_motif/PDBSize.py', sites_dir, output_dir])
            output_file = os.path.join(output_dir, 'PDBSize.txt')
            self.assertTrue(os.path.exists(output_file))
            with open(output_file, 'r') as f:
                output_lines = f.readlines()
                output_data = [line.strip() for line in output_lines if line]
            self.assertEqual(len(output_data), len(os.listdir(sites_dir)))


class MultipleSiteAlignmentTestCase(TestCase):

    def test_atp_subset_alignment(self):
        dirname = os.path.dirname(__file__)
        with tempfile.TemporaryDirectory() as output_dir:
            pairs_list_file = f'{output_dir}/PairList.txt'
            pdb_size_file = f'{output_dir}/PDBSize.txt'
            align_output_file = f'{output_dir}/align_output.txt'
            sites_dir = f'{dirname}/test_data/ATP_SUBSET'
            
            reference_pdb = os.listdir(sites_dir)[0]
            cmd = f'{dirname}/bin/site-motif {reference_pdb} {sites_dir} {output_dir}'
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
