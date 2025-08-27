import numpy as np
from oracles.oracle_component import OracleComponent
from oracles.dataclass import OracleComponentParameters
from rdkit.Chem import Mol
from rdkit import Chem
from oracles.fujiwara.pipeline_activation_energy import calculate_energy


class ActivationEnergy(OracleComponent):
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
            return (calculate_energy(smiles))
        except Exception:
            return 1000