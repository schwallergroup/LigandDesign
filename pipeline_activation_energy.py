from AaronTools.geometry import Geometry
from AaronTools.atoms import Atom
from AaronTools.theory import*
from AaronTools.fileIO import FileReader
from openbabel import openbabel
from openbabel import pybel
import subprocess
from pathlib import Path
import re
import click
import shutil

# Constants to be used
benzene_electronic_energy = -15.90075765 # Calculated at XTB level of theory with DDCOSMO(benzene)
benzene_free_energy = -15.83967938 # Calculated at XTB level of theory with DDCOSMO (benzene)
T = 383 # In K
R = 1.98720425864083 / 1000 # In kcal/mol
conversion_eh_kcal = 627.5094740631 # From hartree to kcal/mol
k_b = 1.380649e-23 # In SI units
h = 6.62607015e-34
frequency_factor = k_b * T / h


def tag_catalyst_atoms (catalyst : Geometry, catalyst_type : str ) -> None :
    '''
    Adds tags to the geometry object corresponding to a catalyst.

    This function receives a geometry object representing a catalyst
    and the type of catalyst according to ACS Catal. 2025, 15, 5503−5510.
    It adds tags but does not return anything.

    Parameters:
    ----------
        catalyst : AaronTools.geometry.Geometry
            Geometry object representing the catalyst
        catalyst_type : str
            Type of catalyst according to ACS Catal. 2025, 15, 5503−5510

    Returns:
    ----------
        None
            This function modifies the geometry object corresponding to the 
            catalyst and does not return any value

    '''

    catalyst.detect_components() # Adds tags 'center', 'key' for atoms bonded to metal, 'backbone' and 'sub-x'.
    atom_list_cat = catalyst.atoms # Defines the list of atoms

    # Tags the atoms in the type 0 catalyst

    if catalyst_type == '0': # Corresponds to catalyst type 0 in the paper

        # Add tag corresponding to the each one of the two ligands

        for i in range(1,len(atom_list_cat)): # 0 is Pd
            if i <= 11: # Atoms that belong to L1
                atom_list_cat[i].add_tag('L1')
            else: # Atoms that belong to L2
                atom_list_cat[i].add_tag('L2')
    
    # Tags the atoms in the type ts (transition state) catalyst

    else: # Corresponds to catalyst type ts in the paper

        # Add tag corresponding to the each one of the two ligands and to the benzene substrate

        for i in range(1,len(atom_list_cat)): # 0 is Pd
            if i <= 11: # Atoms that belong to L2
                atom_list_cat[i].add_tag('L2')
            elif 12 <= i <= 22: # Atoms that belong to L1
                atom_list_cat[i].add_tag('L1')
            else: # Atoms belong to the substrate
                atom_list_cat[i].add_tag('benzene_substrate')

        atom_list_cat[13].add_tag('non_bonded_key') # Ads a tag on the non-coordinating O atom of L1

def old_key_atoms (catalyst : Geometry, catalyst_type : str) -> tuple[list[Atom], list[Atom]]:
    '''
    This function returns the key atoms needed for the substitution
    of the two original ligands.

    This function receives a geometry object representing a catalyst
    and the type of catalyst according to the ACS Catal. 2025, 15, 5503−5510.
    The catalyst must already have the 'key' and the 'non_bonded_key' tags.
    It returns two lists containing the key N and O atoms corresponding
    to the original catalyst. These atoms are needed for an effective
    ligand substitution

    Parameters:
    ----------
        catalyst : AaronTools.geometry.Geometry
            Geometry object representing the catalyst
        catalyst_type : str
            Type of catalyst according to ACS Catal. 2025, 15, 5503−5510

    Returns:
    ----------
        tuple (atom_old_L1 , atom_old_L2)
            A tuple containing two lists

                atom_old_L1 : list[AaronTools.atoms.Atom]
                    List containing key atoms needed for the subtitution of ligand L1
                atom_old_L2 : list[AaronTools.atoms.Atom]
                    List containing key atoms needed for the subtitution of ligand L2
            
    '''

    atom_list_cat = catalyst.atoms # Gets the list of atoms from the analyzed geometry object

    if catalyst_type == '0': # For the catalyst of type 0 only the 'key' tag matters

        atom_old_L1 = [atom_list_cat[i] for i in range(len(atom_list_cat)) 
                   if 'L1' in atom_list_cat[i].tags and 'key' in atom_list_cat[i].tags]
        atom_old_L2 = [atom_list_cat[i] for i in range(len(atom_list_cat)) 
                    if 'L2' in atom_list_cat[i].tags and 'key' in atom_list_cat[i].tags]

    else: # For the TS type of catalysts, also the 'non_bonded_key' is important
        atom_old_L1 = [atom_list_cat[i] for i in range(len(atom_list_cat)) 
                   if 'L1' in atom_list_cat[i].tags and ('key' in atom_list_cat[i].tags
                                                                  or 'non_bonded_key' in atom_list_cat[i].tags)]
        atom_old_L2 = [atom_list_cat[i] for i in range(len(atom_list_cat)) 
                    if 'L2' in atom_list_cat[i].tags and 'key' in atom_list_cat[i].tags]

    # In order to obtain the correct trans isomer of the catalyst, the lists need to be reverted
    atom_old_L1.reverse()
    atom_old_L2.reverse()

    return atom_old_L1, atom_old_L2


def monodeprotonate_donor_atom (ligand : Geometry, atom : Atom) -> None :
    '''
    This function monodeprotonates a given atom of a ligand

    This function receives a geometry object representing a
    ligand and an atom object representing a certain atom of
    the ligand. It monodeprotonates the atom if this one is bound
    to at least one H atom

    Parameters:
    ----------
        ligand : AaronTools.geometry.Geometry 
            Geometry object representing a ligand
        atom : AaronTools.atoms.Atom
            Atom object representing an atom part of the ligand
    
    Returns:
    ----------
        None
            This function modifies the geometry object corresponding to the 
            ligand and does not return any value

    '''

    for neighbor in atom.connected: # Iterates through the atoms bonded to the coordinating atom
        if neighbor.element == 'H': # Verifies if it is an H atom and removes it if so
            ligand.atoms.remove(neighbor)
            break

def tag_ligand_atoms (ligand : Geometry, atom : Atom) -> None:
    '''
    Adds tags to the geometry object corresponding to a ligand.

    This function receives a geometry object representing a ligand
    and the atom on which it has to add tags. Those tags are
    necessary for a proper ligand substitution. It adds tags and
    does not return anything.

    Parameters:
    ----------
        ligand : AaronTools.geometry.Geometry 
            Geometry object representing a ligand
        atom : AaronTools.atoms.Atom
            Atom object representing an atom part of the ligand
    
    Returns:
    ----------
        None
            This function modifies the geometry object corresponding to the 
            ligand and does not return any value

    '''

    atom.add_tag('key', 'new') # Adds tags on the key atom of the ligand
    monodeprotonate_donor_atom(ligand, atom) # Monodeprotonates the atom in case it is bound to H


def attach_ligand (SMILES : str, N_atom_index : int) -> None:
    '''
    Substitutes the ligand in catalysts from ACS Catal. 2025, 15, 5503−5510 with a new ligand.

    Given the SMILES of a ligand, this function replaces
    the original ligand in the catalyst structure with the 
    new ligand and writes the modified geometry to an .xyz file.

    Parameters:
    ----------
        SMILES : str
            SMILES string of the ligand to be added
        N_atom_index : int
            Integer number representing the index (according to AaronTools)
            of the coordinating N atom in the catalyst

    Returns:
    ----------
        None
            This function attached the SMILES of the 
            proposed ligand to the catalyst, writes the new
            geometry in an .xyz file and doesn't output anything

    '''

    original_geometry_directory = Path(r'/home/octavian/original_geometries') # Path to the direcotry with original geometries

    # Creates a directory for the SMILES analyzed, in which all the calculations will be performed
    new_geometry_directory = Path(r'/data/octavian/', f'{SMILES}') # Path to the directory with new geometries
    if not(new_geometry_directory.exists()):
        new_geometry_directory.mkdir(exist_ok = True)

    
    for geometry_file in original_geometry_directory.iterdir(): # Iterates through the geometry files
        
        original_initial_cat = Geometry(f'{geometry_file}') # Gets catalyst from the paper

        filename = geometry_file.name # Gets the name of the file where the geometry will be saved, '0.xyz' for example
        type_of_compound = geometry_file.stem.strip() # Gets the type of catalyst, '0' or 'ts'
        
        tag_catalyst_atoms(original_initial_cat, type_of_compound) # Adds tags to the atoms
        atom_old_L1, atom_old_L2 = old_key_atoms(original_initial_cat, type_of_compound) # Gets old keys for each type of copy

        # map_ligand function modifies the geometry object corresponding to
        # the new ligand; as two substitutions with the same ligand are 
        # needed, this one is defined twice 
        new_ligand_1 = Geometry.from_string(SMILES)
        new_ligand_2 = Geometry.from_string(SMILES)

        # Tags are added, for the moment knowing exactly on which atoms; also, certain atoms may be deprotonated
        tag_ligand_atoms(new_ligand_1, new_ligand_1.atoms[0])
        tag_ligand_atoms(new_ligand_1, new_ligand_1.atoms[N_atom_index])
        tag_ligand_atoms(new_ligand_2, new_ligand_2.atoms[0])
        tag_ligand_atoms(new_ligand_2, new_ligand_2.atoms[N_atom_index])

        # Creates catalyst with L1 and L2 in original order
        original_initial_cat.map_ligand(new_ligand_1, atom_old_L1, minimize = False)
        original_initial_cat.map_ligand(new_ligand_2, atom_old_L2, minimize = False)
        
        # Saves the geometries in new files
        outfile = new_geometry_directory / filename 
        original_initial_cat.write(outfile = outfile)

def replace_cpcm_with_ddcosmo(input_file : str) -> None:
    '''
    This function replaces the CPCM solvent model
    with DDCOSMO in an ORCA input file.

    This function takes as an argument the path to
    an ORCA input file, that contains as solvation
    model CPCM. It replaces CPCM in the file with 
    DDCOSMO and doesn't return anything.

    Parameters:
    ----------
        input_file : str
            String path to the ORCA input file that is modified

    Returns:
    ----------
        None
            This function replaces CPCM with DDCOSMO in an
            ORCA input file and doesn't return anything.

    '''
    # Open the input file for reading
    with open(input_file, 'r') as file:
        content = file.read()
    
    # Replace 'CPCM' with 'DDCOSMO'
    content = content.replace("CPCM", "DDCOSMO")
    
    # Open the input file for writing
    with open(input_file, 'w') as file:
        file.write(content)

def extract_K_indices(filename : str) -> list[str]:
    '''
    This function returns the list of key atoms
    indicated in an .xyz file 

    This function takes as argument the string path
    to an .xyz file. It returns a list of strings,
    where each string represents the index of a 
    key atom (coordinating, or not coordinating such as
    the O atom in the ts). The key indices are marked 
    in the 2nd line of the .xyz file by K

    Parameters:
    ----------
        filename : str
            String path to the .xyz file from which indices are extracted

    Returns:
    ----------
        indices : list[str]
            List of strings, where each string represents the index of a 
            key atom in the catalyst analyzed

    '''
    with open(filename, 'r') as f:
        for line in f:
            if 'K:' in line:
                # Find the K: portion and everything after it
                # until the next tag or end of line
                match = re.search(r'K:([^\nLC]*)', line) # Finds the K field (until it encounters L or C or endspace and saves in group(1))
                if match:
                    k_string = match.group(1)
                    chunks = k_string.strip().split(';') # Splits string into chunks according to ;
                    indices = []
                    for chunk in chunks:
                        chunk = chunk.strip() # Removes whitespaces in case there are still
                        if not chunk: # Handles case for no chunk
                            continue
                        if ',' in chunk: # Handles cases like 19,20
                            subparts = chunk.split(',')
                            for sp in subparts:
                                indices.append(sp)
                        elif '-' in chunk: # Handles cases like 19-21
                            start, end = map(int, chunk.split('-'))
                            indices.extend([str(i) for i in range(start, end + 1)])
                        else:
                            indices.append(chunk)
                    return indices
    return []

def automated_files_XTB (SMILES : str) -> None:
    '''
    This function creates ORCA input files and 
    launches ORCA calculations for a given ligand.

    This function receives as argument a SMILES string,
    creates ORCA input files for conformer search and
    geometry optimization (uses lowest energy conformer 
    as input for the geometry optimization) for the 
    0 and ts compounds corresponding to the SMILES and 
    launches the ORCA calculations. It doesn't output
    anything.

    Parameters:
    ----------
        SMILES : str
            String of the SMILES for which ORCA calculations are performed

    Returns:
    ----------
        None
            This function creates ORCA input files, launches ORCA calculations
            and doesn't output anything

    '''
    smiles_path = Path(r'/data/octavian/',f'{SMILES}') # Goes to the direcotry containing the initial geometry
    files = [f for f in smiles_path.iterdir() if f.is_file()] # List of geometry files corresponding to compound 0 and to ts

    for file in files: # Iterates through compounds 0 and ts

        #Creates a directory for the compound for which calculations are performed if the file does not already exist
        computation_directory_path = Path(smiles_path, f'{file.stem}_comput')
        if not(computation_directory_path.exists()):
            computation_directory_path.mkdir(exist_ok= True)

        catalyst = Geometry(f'{file}') # Creates a geometry object for the compound being analyzed 
        index = file.stem.strip() # Type of compound being analyzed

        # Define theory level for geometry optimization
        method = 'XTB2'
        solvent = ImplicitSolvent('CPCM', 'benzene')

        if index != 'ts':
                
            # Fixes the coordinating atoms (the key atoms) and the Pd atom for the conformational search
            fixed_atoms = extract_K_indices(f'{file}')
            fixed_atoms.append('1')
            job_list_conformer = [ConformerSearchJob(constraints = {'atoms': fixed_atoms})]

            # Defines the theory level of the conformational search
            catalyst_theory_level_conformer = Theory(
                method = method,
                job_type = job_list_conformer
            )
            ORCA_BLOCKS_normal = {"geom": ["MaxIter 500"]} 

            # Defines the input and output files for the conformational search and writes the inp file
            outfile_input_conformer = Path(computation_directory_path, 'conformer.inp')
            outfile_output_conformer = Path(computation_directory_path, 'conformer.out')
            catalyst.write(
                outfile = f'{outfile_input_conformer}',
                theory = catalyst_theory_level_conformer,
                processors = 8,
                memory = 12,
                ORCA_BLOCKS = ORCA_BLOCKS_normal
            )

            #Launch ORCA
            with open(outfile_output_conformer, 'w') as out_f:
                orca_run = subprocess.run(f"/data/octavian/software/orca_6_1_0_linux_x86-64_shared_openmpi418/orca conformer.inp", cwd = computation_directory_path, 
                                        shell=True, stdout=out_f, stderr=subprocess.STDOUT)
                if orca_run.returncode != 0:
                    raise RuntimeError(f"ORCA failed with return code {orca_run.returncode}")

            # Creates the directory where the geometry optimization takes place
            computation_directory_path_lowest_conf = Path(computation_directory_path, 'energy')
            if not(computation_directory_path_lowest_conf.exists()):
                computation_directory_path_lowest_conf.mkdir(exist_ok=True)

            # Gets the geometry object for the most stable conoformer
            most_stable_conformer_path = Path(computation_directory_path, 'conformer.globalminimum.xyz')
            most_stable_conformer = Geometry(f'{most_stable_conformer_path}')

            # Theory level for the geometry optimization; nothing is constrained
            job_list = [OptimizationJob(), FrequencyJob(temperature = 383)]
            
            catalyst_theory_level = Theory(
                method = method,
                solvent = solvent,
                job_type = job_list
            )
            ORCA_BLOCKS_normal = {"geom": ["MaxIter 500"]}

            #Create ORCA input file
            outfile_input = Path(computation_directory_path_lowest_conf, "catalyst.inp")
            outfile_output = Path(computation_directory_path_lowest_conf, "catalyst.out")
            most_stable_conformer.write(
                outfile = f'{outfile_input}',
                theory = catalyst_theory_level,
                processors = 8,
                memory = 12,
                ORCA_BLOCKS = ORCA_BLOCKS_normal
            )
            replace_cpcm_with_ddcosmo(f'{outfile_input}') # Replaces CPCM with DDCOSMO

            #Launch ORCA
            with open(outfile_output, 'w') as out_f:
                orca_run = subprocess.run(f"/data/octavian/software/orca_6_1_0_linux_x86-64_shared_openmpi418/orca catalyst.inp", cwd = computation_directory_path_lowest_conf, 
                                        shell=True, stdout=out_f, stderr=subprocess.STDOUT)
                if orca_run.returncode != 0:
                    raise RuntimeError(f"ORCA failed with return code {orca_run.returncode}")
                
        else: # Handles the transition state

            # Extracts the tags of the key atoms
            # (coordinating atoms, including coordinating C atom from benzene, + O atom to which H transfer takes place)
            # from the geometry files
            fixed_atoms = extract_K_indices(f'{file}')

            # Apart from the key atoms, the Pd, the benzene C atoms (apart the one coordinating) 
            # and the H atom transferred need fixed 
            for i in range(1,9): 
                if not(str(i) in fixed_atoms): 
                    fixed_atoms.append(str(i))

            # Defines the theory level for the conformational search of the ts
            conformer_search_theory_level = Theory(
                method = method,
                job_type = [ConformerSearchJob(constraints = {'atoms' : fixed_atoms})]
            )
            ORCA_BLOCKS_normal = {"geom": ["MaxIter 500"]} 

            # Create ORCA input and output files for the conformational search and write the inp
            outfile_input_conformer = Path(computation_directory_path, 'conformer.inp')
            outfile_output_conformer = Path(computation_directory_path, 'conformer.out')
            catalyst.write(
                outfile = f'{outfile_input_conformer}',
                theory = conformer_search_theory_level,
                processors = 8,
                memory = 12,
                ORCA_BLOCKS = ORCA_BLOCKS_normal
            )

            #Launch ORCA
            orca_run = subprocess.run(f"/data/octavian/software/orca_6_1_0_linux_x86-64_shared_openmpi418/orca conformer.inp", cwd = computation_directory_path,
                                    shell=True, stdout=open(outfile_output_conformer, 'w'), stderr=subprocess.STDOUT)
            if orca_run.returncode != 0:
                raise RuntimeError(f"ORCA failed with return code {orca_run.returncode}")

            # Creates the directory where the geometry optimization takes place
            computation_directory_path_lowest_conf = Path(computation_directory_path, 'energy')
            if not(computation_directory_path_lowest_conf.exists()):
                computation_directory_path_lowest_conf.mkdir(exist_ok=True)

            # Gets the geometry object for the most stable conoformer
            most_stable_conformer_path = Path(computation_directory_path, 'conformer.globalminimum.xyz')
            most_stable_conformer = Geometry(f'{most_stable_conformer_path}')

            # Defines the theory level for the constrained geometry optimization of the ts
            ts_constrained_opt_theory_level = Theory(
                method = method,
                solvent = solvent,
                job_type = [OptimizationJob(constraints = {"atoms": fixed_atoms}), FrequencyJob(temperature = 383)]
            )
            ORCA_BLOCKS_normal = {"geom": ["MaxIter 500"]}

            # Create ORCA input file for the constrained opt of the ts and run the optimization
            constrained_outfile_input = Path(computation_directory_path_lowest_conf, "ts_constrained.inp")
            constrained_outfile_output = Path(computation_directory_path_lowest_conf, "ts_constrained.out")
            most_stable_conformer.write(
                outfile = f'{constrained_outfile_input}',
                theory = ts_constrained_opt_theory_level,
                processors = 8,
                memory = 12,
                ORCA_BLOCKS = ORCA_BLOCKS_normal
            )
            replace_cpcm_with_ddcosmo(f'{constrained_outfile_input}') # Replaces CPCM with DDCOSMO in the ORCA input file

            #Launch ORCA
            orca_run = subprocess.run(f"/data/octavian/software/orca_6_1_0_linux_x86-64_shared_openmpi418/orca ts_constrained.inp", cwd = computation_directory_path_lowest_conf,
                                    shell=True, stdout=open(constrained_outfile_output, 'w'), stderr=subprocess.STDOUT)
            if orca_run.returncode != 0:
                raise RuntimeError(f"ORCA failed with return code {orca_run.returncode}")

def verify_big_changes(SMILES: str, type_of_compound: str) -> bool:
    '''
    This function verifies if big changes such as
    bond breaks have taken place during the 
    optimization process.

    This function takes as argument a SMILES string
    corresponding to the ligand that has been optimized
    and the type of compound that has been optimized (either
    type 0 or type ts). It returns True if bond breaks have
    taken place and False if no major changes have taken
    place.

    Parameters:
    ----------
        SMILES : str
            String of the SMILES for which the optimzation has been performed
        type_of_compound : str
            String of the type of compound that is being analyzed. Needs to be
            either 0.xyz or ts.xyz

    Returns:
    ----------
        True : bool
            If major changes, such as bond breaking, have taken place
        False : bool
            If no major changes have taken place
    '''

    SMILES_dir = Path(r'/data/octavian/',f'{SMILES}') # Path to the SMILES directory
    original_geometry_file = Path(SMILES_dir, type_of_compound) # Path to the file of the unoptimized geometry
    # Path to the file of the optimized geometry
    if type_of_compound == '0.xyz':
        optimized_geometry_file = Path(SMILES_dir, '0_comput/energy/catalyst.xyz')
    elif type_of_compound == 'ts.xyz':
        optimized_geometry_file = Path(SMILES_dir, 'ts_comput/energy/ts_constrained.xyz')

    # Geometry objects are defined for the unoptimized and optimized compounds
    original_geometry = Geometry(f'{original_geometry_file}')
    optimized_geometry = Geometry(f'{optimized_geometry_file}')

    # The list of atoms for each of the two geometries are extracted
    atoms_original = original.atoms
    atoms_optimized = optimized.atoms

    # The geometry files of the original and optimized geometries
    # should have the atoms indicated in the same order, and therefore
    # the atoms should have the same indices.
    # Thus, atoms with the same indices in both lists should have the 
    # same neighbors. So, by sorting the list of indices of the neighbors
    # of two atoms with the same index, the list should be the same if
    # no major chages have taken place.

    for i in range(len(atoms_original)): # Iterates through the atom lists

        neighbors_original = atoms_original[i].connected # Gets the neighbor list of the atom analyzed in the original geometry
        original_indices = [] # List that will contain the indices of the neighbors
        for neighbor_atom_original in neighbors_original: # Iterates through the neighbor atoms
            original_indices.append(int(neighbor_atom_original.name)) # Appends the index of the neighbor atom analyzed
        original_indices.sort() # Sorts the neighbor list in an increasing order

        neighbors_optimized = atoms_optimized[i].connected # Gets the neighbor list of the atom analyzed in the optimized geometry
        optimized_indices = [] # List that will contain the indices of the neighbors
        for neighbor_atom_optimized in neighbors_optimized: # Iterates through the neighbor atoms
            optimized_indices.append(int(neighbor_atom_optimized.name)) # Appends the index of the neighbor atom analyzed
        optimized_indices.sort() # Sorts the neighbor list in an increasing order
        
        if original_indices != optimized_indices: # Compares the two lists
            return(True) # Big changes have taken place
    return (False)

                   
def extract_energies(SMILES : str) -> float:
    '''
    This function returns the activation energy of a 
    ligand in the reaction described in
    ACS Catal. 2025, 15, 5503−5510.

    This function takes as argument a string of a SMILES
    corresponding to a ligand. It parses the ORCA output
    files corresponding to that ligand, writes the
    coordinates of the precomplex (compound 0) and the ts
    in a .log file and outputs the Gibbs free activation
    energy [kcal/mol] corresponding to that ligand.

    Parameters:
    ----------
        SMILES : str
            String of the SMILES for which the activation energy is returned

    Returns:
    ----------
        activation_energy : float
            Gibbs free activation energy corresponding to the ligand for
            which the SMILES is given as input

    '''

    # Path to the directory that contains the energy calculations for the given SMILES
    # and path to the .log file where the coordinates need to be output
    SMILES_dir = Path(r'/data/octavian/',f'{SMILES}') 
    output_file = '/data/octavian/saturn.log'

    with open (output_file, 'a') as f:

        compound_0_file = Path(SMILES_dir, '0_comput/energy/catalyst.out') # Path to ORCA output file of compound 0
        compound_ts_file = Path(SMILES_dir, 'ts_comput/energy/ts_constrained.out') # Path to ORCA output file of compound ts
        compound_0_geometry = Path(SMILES_dir, '0_comput/energy/catalyst.xyz') # Path to optimized geometry of compound 0
        compound_ts_geometry = Path(SMILES_dir, 'ts_comput/energy/ts_constrained.xyz') # Path to optimized geometry of compound ts

        # Extracts the content from the output files of the compounds
        compound_0_output_content = FileReader(f'{compound_0_file}', just_geom = False)
        compound_ts_output_content = FileReader (f'{compound_ts_file}', just_geom = False)

        if not(compound_0_output_content['finished']): # Verrifies whether compound 0 finished
            f.write(f'{SMILES} Compound 0 did not finish the optimization \n')
            f.write(f'----------------------------- \n')
        elif compound_0_output_content['error']: # Verrifies whether compound 0 encountered any error during the calculations
            f.write(f'{SMILES} Calculation of compound 0 encountered an error \n')
            f.write(f'----------------------------- \n')
        elif not(compound_ts_output_content['finished']): # Verrifies whether compound ts finished
            f.write(f'{SMILES} Compound ts did not finish the optimization \n')
            f.write(f'-----------------------------\n')
        elif compound_ts_output_content['error']: # Verrifies whether compound ts encountered any error during the calculations
            f.write(f'{SMILES} Calculation of compound ts encountered an error \n')
            f.write(f'-----------------------------\n')
        elif verify_big_changes(f'{SMILES}', '0.xyz'): # Verifies whether compound 0 suffered major changes during optimization
            f.write(f'{SMILES} Calculation of compound 0 resulted in major geometry changes \n)')
        elif verify_big_changes(f'{SMILES}', 'ts.xyz'): # Verifies whether compound ts suffered major changes during optimization
            f.write(f'{SMILES} Calculation of compound ts resulted in major geometry changes \n)')    
        else:

            # Extracts the Gibbs free energy of compound 0 and of the ts (both in Eh)
            free_energy_0 = compound_0_output_content['free_energy']
            free_energy_ts = compound_ts_output_content['free_energy']

            # Calculates the Gibbs free activation energy [kcal/mol]
            activation_energy = (free_energy_ts - free_energy_0 - benzene_free_energy) * conversion_eh_kcal

            # Writes the activation energy in the .log file
            f.write(f'Activation energy for the ligand {SMILES} = {activation_energy} kcal/mol \n')
            f.write(f'----------------------------- \n')

        # Writes the coordinates of the precatalyst in the .log file 
        with open (f'{compound_0_geometry}', 'r') as h:
            for line in h:
                f.write(line)
        f.write(f'----------------------------- \n')

        # Writes the coordinates of the ts in the .log file 
        with open (f'{compound_ts_geometry}', 'r') as g:
            for line in g:
                f.write(line)  
        f.write(f'----------------------------- \n')
        f.write(f'----------------------------- \n')

    shutil.rmtree(f'{SMILES_dir}') # Removes the entire computation directory 
    return (activation_energy)   

def rearrange_smiles (SMILES : str) -> str:
    '''
    This function rearranges a 2-hydroxypyridine SMILES
    so that it begins with the O atom.

    This function takes as argument a SMILES string corresponding
    to a molecule containing the 2-hydroxypyridine pattern. It 
    outputs a SMILES string of the same molecule but which
    begins with the O atom corresponding to the 2-hydroxy group.

    Parameters:
    ----------
        SMILES : str
            String of the SMILES which is modified

    Returns:
    ----------
        arranged_smiles : str
            SMILES of the same molecule but which begins
            with the O atom of the 2-hydroxy group

    '''

    mol = pybel.readstring('smi', SMILES) # Reads smiles
    
    try:

        pattern = pybel.Smarts('[OH]c1ncccc1') # Pattern of 2-hydroxypyridine that is recognizable to the pipeline

        # Pattern findall returns a list of tuples
        # In the tuple, the indices of the atoms that form the pattern are found
        # Search for the 2-hydroxypyridine pattern and save the index of the O atom 
        # Of the 2-hydroxy group
        pattern_idx = pattern.findall(mol)[0][0] 

        rearranger = openbabel.OBConversion() # Defines the rearranger
        rearranger.SetInAndOutFormats('smi', 'smi') # Sets input and output format as SMILES

        # 'f' ensures that first atom is the one wanted, in this case O of the 2-hydroxy group
        rearranger.AddOption('f', openbabel.OBConversion.OUTOPTIONS, str(pattern_idx)) 

        outmol = openbabel.OBMol() # Defines an output molecule

        rearranger.ReadString(outmol, SMILES) # Rearranges the SMILES into the output molecule
        arranged_smiles = rearranger.WriteString(outmol).strip() # Extracts the rearranged SMILES

        return(arranged_smiles) # Returns the rearranged SMILES
    except:
        return 0

def get_N_index (SMILES : str) -> int:
    '''
    This function function returns the index of the pyridine N
    atom (according to AaronTools) in a 2-hydroxypyridine

    This function takes as input a string corresponding to
    the ordered SMILES (begins with the O atom of the
    2-hydroxy group) of a 2-hydroxypyridine and returns the
    index of the pyridine N atom (for the AaronTools Geometry 
    objects). It is based on the fact that the indices of the 
    atoms in the AaronTools Geometry objects are the same to 
    the order in which atoms appear in the SMILES string 
    corresponding to the molecule.

    Parameters:
    ----------
        SMILES : str
            String of the SMILES for which the index of the pyridine N atom is found

    Returns:
    ----------
        i : int
            Index of the pyridine N atom according to AaronTools 

    '''

    mol_geometry = Geometry.from_string(SMILES)

    # Since the molecule starts with the O atom (index 0), 
    # the second one atom (index 1) is necessarily the C
    # atom bound to the pyridine N atom
    for neighbour in mol_geometry.atoms[1].connected:
        if neighbour.element == 'N': # Finds the pyridine N atom
            for i, atom in enumerate(mol_geometry.atoms): # Searches for the index of the pyridine N atom
                if atom == neighbour: # Finds the index of the pyridine N atom
                    return(i) # Returns the index of the pyridine N atom 


@click.group()
def cli():
    """A script with multiple functions accessible via command line."""
    pass
@cli.command()
@click.argument('smiles')
def calculate_energy(smiles : str) -> float:
    '''
    This function attaches a ligand to Pd,
    computes its activation energy in the reaction
    described in ACS Catal. 2025, 15, 5503−5510
    and return the Gibbs free activation energy

    This function takes as input a SMILES string 
    corresponding to a 2-hydroxypyridine. The ligand
    is attached (in two steps, from one side and from
    the other) to Pd in order to give the trans isomer
    of the catalyst described in ACS Catal. 2025, 15, 5503−5510.
    ORCA calculations are performed for the ligand, geometry of the
    initial state of the catalyst and the catalyst in the ts corresponding
    to the RDS are written in a .log file and the Gibbs free activation
    energy is returned 
    
    Parameters:
    ----------
        smiles : str
            String of the SMILES for which the Gibbs free activation
            energy is returned

    Returns:
    ----------
        activation_energy : float
            Gibbs free activation energy corresponding to the analyzed ligand 
    '''

    arranged_smiles = rearrange_smiles(smiles) # Rearranges the SMILES such that it begins with the O atom of the 2-hydroxy group
    N_atom_index = get_N_index(arranged_smiles) # Determines the index of the pyridine N atom
    attach_ligand(arranged_smiles, N_atom_index) # Attaches two ligands of the same kind to a Pd atom and saves the geometries
    automated_files_XTB(arranged_smiles) # Performs ORCA calculations on the geometries (saved above) of the precomplex and ts
    activation_energy = extract_energies(arranged_smiles) # Extracts the Gibbs free activation energy and writes the coordinates in a .log file
    return(activation_energy) # Returns the activation energy 

if __name__ == "__main__":
    cli()