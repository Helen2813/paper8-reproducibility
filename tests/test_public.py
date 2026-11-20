import sys, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'portable'))
sys.path.insert(0,str(ROOT/'scripts'))
import portable_cli as cli
import public_code_preflight as gate

class PublicCodeTests(unittest.TestCase):
    def test_schema(self):
        conf=cli.ensure_bundle()
        self.assertEqual(len(conf['resources']),9)
        self.assertEqual(set(conf['external_references']),{'s2','s3'})
    def test_private_inputs_absent(self):
        conf=cli.load_conf()
        for k,resource in conf['resources'].items():
            self.assertFalse((ROOT/resource['portable_filename']).exists(),k)
        for spec in conf['external_references'].values():
            self.assertFalse((ROOT/spec['portable_filename']).exists())
    def test_corrupt_input_hold(self):
        conf=cli.load_conf()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'combat_pre.tsv';p.write_bytes(b'bad data')
            with self.assertRaisesRegex(RuntimeError,'HOLD_PORTABLE_INPUT_HASH'):
                cli.check_inputs(conf,Path(d),['combat_pre'])
    def test_missing_reference_hold(self):
        conf=cli.load_conf()
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(RuntimeError,'HOLD_PRIVATE_REFERENCE_MISSING'):
                cli.check_external_references(conf,Path(d))
    def test_corrupt_reference_hold(self):
        conf=cli.load_conf()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'reference_S2.csv';p.write_bytes(b'bad data')
            with self.assertRaisesRegex(RuntimeError,'HOLD_PRIVATE_REFERENCE_HASH'):
                cli.check_external_references(conf,Path(d))
    def test_refuses_git_as_output(self):
        with self.assertRaisesRegex(RuntimeError,'HOLD_UNSAFE'):
            cli.safe_out(ROOT/'results',[ROOT])
    def test_static_gate(self):
        self.assertEqual(gate.audit(),[])
    def test_worker_source_parses(self):
        import ast
        for name in ('portable/portable_cli.py','portable/classifier_replay.py'):
            ast.parse((ROOT/name).read_text())
    def test_gitignore_blocks_private_tables(self):
        x=(ROOT/'.gitignore').read_text()
        for line in ('*.csv','*.tsv','.envs/','.mamba/','/data/'):
            self.assertIn(line,x)

if __name__=='__main__':unittest.main()
