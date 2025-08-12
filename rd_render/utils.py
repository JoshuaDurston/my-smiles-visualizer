from rdkit import Chem
import sys
import os

def safe_mol_from_smiles(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        return None

    # Backup original stdout and stderr
    stdout_fileno = sys.stdout.fileno()
    stderr_fileno = sys.stderr.fileno()

    # Redirect both to null temporarily
    with open(os.devnull, 'w') as devnull:
        old_stdout = os.dup(stdout_fileno)
        old_stderr = os.dup(stderr_fileno)
        os.dup2(devnull.fileno(), stdout_fileno)
        os.dup2(devnull.fileno(), stderr_fileno)

        try:
            mol = Chem.MolFromSmiles(smiles)
        finally:
            # Restore original stdout and stderr
            os.dup2(old_stdout, stdout_fileno)
            os.dup2(old_stderr, stderr_fileno)
            os.close(old_stdout)
            os.close(old_stderr)

    return mol