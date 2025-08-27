import numpy as np
from oracles.oracle_component import OracleComponent
from oracles.dataclass import OracleComponentParameters
from rdkit.Chem import Mol
from rdkit import Chem



class Hydroxypyridine(OracleComponent):
    def __init__(self, parameters: OracleComponentParameters):
        super().__init__(parameters)

    def __call__(self, mols: np.ndarray[Mol]) -> np.ndarray[float]:
        return np.vectorize(self._compute_property)(mols)

    def _compute_property(self, mol: Mol) -> float:
        """
        Wrapper function in case of exceptions.
        """
        try:
            hydroxypyridine_pattern = Chem.MolFromSmarts("[OH]c1ccccn1")
            hydroxy_pattern = Chem.MolFromSmarts('[OH]')
            amine_pattern = Chem.MolFromSmarts('[NH2]')
            if (len(mol.GetSubstructMatches(hydroxypyridine_pattern)) == 1 and len(mol.GetSubstructMatches(hydroxy_pattern)) == 1 and (not mol.HasSubstructMatch(amine_pattern))):
                return 1.0
            return 0.0
        except Exception:
            return 0.0