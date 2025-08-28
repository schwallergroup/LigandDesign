import numpy as np
from oracles.oracle_component import OracleComponent
from oracles.dataclass import OracleComponentParameters
from rdkit.Chem import Mol
from rdkit import Chem



class Neutral(OracleComponent):
    def __init__(self, parameters: OracleComponentParameters):
        super().__init__(parameters)

    def __call__(self, mols: np.ndarray[Mol]) -> np.ndarray[float]:
        return np.vectorize(self._compute_property)(mols)

    def _compute_property(self, mol: Mol) -> float:
        """
        Wrapper function in case of exceptions.
        """
        try:
            smiles = Chem.MolToSmiles(mol)

            positive_charges = smiles.count('+')
            negative_charges = smiles.count('-')

            if positive_charges == negative_charges:
                return 1.0
            return 0.0
        except Exception:
            return 0.0