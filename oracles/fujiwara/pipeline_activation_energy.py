from AaronTools.geometry import Geometry
from AaronTools.atoms import Atom
from AaronTools.theory import*
from AaronTools.fileIO import FileReader
from openbabel import openbabel
from openbabel import pybel
import numpy as np
import subprocess
from pathlib import Path
import re
import click
import shutil

# Constants to be used
benzene_electronic_energy = -15.88624946567 # Calculated at XTB level of theory without solvent
benzene_free_energy = -15.82458903 # Calculated at XTB level of theory without solvent
T = 383 # In K
R = 1.98720425864083 / 1000 # In kcal/mol
conversion_eh_kcal = 627.5094740631 # From hartree to kcal/mol
k_b = 1.380649e-23 # In SI units
h = 6.62607015e-34
frequency_factor = k_b * T / h


def tag_catalyst_atoms (catalyst, catalyst_type):
    '''
    This function receives a geometry object representing a catalyst
    and the type of catalyst according to the Fujiwara paper. It adds
    tags but does not return anything.

    Parameters:
        catalyst (geom): Geometry object representing the catalyst
        catalyst_type (int or string): Type of catalyst according to the paper

    Outputs:
        This function outputs nothing, it just adds tags to the atoms

    '''
    catalyst.detect_components() # Adds tags 'center', 'key' for atoms bonded to metal, 'backbone' and 'sub-x'.
    atom_list_cat = catalyst.atoms
    if catalyst_type == '0': # Corresponds catalyst type 0 in the paper
        # Add tag corresponding to the belonging ligand
        for i in range(1,len(atom_list_cat)): # 0 is Pd
            if i <= 11: # Atoms that belong to L1
                atom_list_cat[i].add_tag('L1')
            else: # Atoms that belong to L2
                atom_list_cat[i].add_tag('L2')
    elif catalyst_type == '1':
        # Add tag corresponding to the belonging ligand
        for i in range(1,len(atom_list_cat)): # 0 is Pd
            if i <= 11: # Atoms that belong to L1
                atom_list_cat[i].add_tag('L1')
            elif 12 <= i <= 22: # Atoms that belong to L2
                atom_list_cat[i].add_tag('L2')
            else: # Atoms belong to the substrate
                atom_list_cat[i].add_tag('benzene_substrate')
    elif catalyst_type == '3':
        # Add tag corresponding to the belonging ligand
        for i in range(1,len(atom_list_cat)): # 0 is Pd
            if i <= 11: # Atoms that belong to L2
                atom_list_cat[i].add_tag('L2')
            elif 12 <= i <= 22: # Atoms that belong to substrate
                atom_list_cat[i].add_tag('benzene_substrate')
            else: # Atoms belong to L1
                atom_list_cat[i].add_tag('L1')
        atom_list_cat[24].add_tag('non_bonded_key') # Ads a tag on the non-coordinating O atom of L1
    else:
        for i in range(1,len(atom_list_cat)): # 0 is Pd
            if i <= 11: # Atoms that belong to L2
                atom_list_cat[i].add_tag('L2')
            elif 12 <= i <= 22: # Atoms that belong to L1
                atom_list_cat[i].add_tag('L1')
            else: # Atoms belong to the substrate
                atom_list_cat[i].add_tag('benzene_substrate')
        atom_list_cat[13].add_tag('non_bonded_key') # Ads a tag on the non-coordinating O atom of L1

def old_key_atoms (catalyst, catalyst_type):
    '''
    This function receives a geometry object representing a catalyst
    and the type of catalyst according to the Fujiwara paper. The catalyst
    must already have the 'key' and the 'non_bonded_key' tags. It returns
    a list containing the key atoms corresponding to the original catalysts.

    Parameters:
        catalyst (geom): Geometry object representing the catalyst
        catalyst_type (int or string): Type of catalyst according to the paper

    Outputs:
        atom_old_L1 (list): List containing the key atoms for the substitution of L1
        atom_old_L2 (list): List containing the key atoms for the substitution of L2
    '''
    atom_list_cat = catalyst.atoms # Gets the list of atoms
    if catalyst_type == '0' or catalyst_type == '1': # For first two types of catalyst, only the 'key' tag matters
        atom_old_L1 = [atom_list_cat[i] for i in range(len(atom_list_cat)) 
                   if 'L1' in atom_list_cat[i].tags and 'key' in atom_list_cat[i].tags]
        atom_old_L2 = [atom_list_cat[i] for i in range(len(atom_list_cat)) 
                    if 'L2' in atom_list_cat[i].tags and 'key' in atom_list_cat[i].tags]
    else: # For the other types of catalysts, also the 'non_bonded_key' is important
        atom_old_L1 = [atom_list_cat[i] for i in range(len(atom_list_cat)) 
                   if 'L1' in atom_list_cat[i].tags and ('key' in atom_list_cat[i].tags
                                                                  or 'non_bonded_key' in atom_list_cat[i].tags)]
        atom_old_L2 = [atom_list_cat[i] for i in range(len(atom_list_cat)) 
                    if 'L2' in atom_list_cat[i].tags and 'key' in atom_list_cat[i].tags]
    return atom_old_L1, atom_old_L2

def return_old_key_atoms (catalyst, catalyst_type, catalyst_copy):
    '''
    This function receives a geometry object representing a catalyst, the 
    type of catalyst according to the Fujiwara paper and the type of copy. 
    The catalyst must already have the 'key' and the 'non_bonded_key' tags.
    It returns a list containing the key atoms needed for the substitution of L1
    and L2 depending on the type of copy used.

    Parameters:
        catalyst (geom): Geometry object representing the catalyst
        catalyst_type (int or string): Type of catalyst according to the paper
        catalyst_copy (int): Type of copy; for a copy of 0 the lists are 
        not modified, for copy 1 the list for L2 is reversed, for copy 2 L1 is
        reversed and for copy 3 both are reversed.

    Outputs:
        atom_old_L1 (list): List containing the key atoms for the substitution of L1
        atom_old_L2 (list): List containing the key atoms for the substitution of L2
    '''
    atom_old_L1, atom_old_L2 = old_key_atoms(catalyst, catalyst_type) # Extracts the key atoms
    if catalyst_copy == 1: # For copy type 1, the key atoms for L2 are reversed compared to the extracted order
        atom_old_L2.reverse()
        #return atom_old_L1, atom_old_L2
    if catalyst_copy == 2: # For copy type 2, the key atoms for L1 are reversed compared to the extracted order
        atom_old_L1.reverse()
        #return atom_old_L1, atom_old_L2
    if catalyst_copy == 3: # For copy type 3, both the key atoms for L1 and L2 are reversed compared to the extracted order
        atom_old_L1.reverse()
        atom_old_L2.reverse()
        #return atom_old_L1, atom_old_L2
    return atom_old_L1, atom_old_L2

def monodeprotonate_donor_atom (ligand, atom):
    '''
    This function receives a geometry object representing a
    ligand and an atom object representing a certain atom of
    the ligand. It monodeprotonates the atom if this one is bound
    to H.

    Parameters:
        ligand: geometry object representing a ligand
        atom: atom object representing an atom part of the ligand
    '''
    for neighbor in atom.connected: # Iterates through the atoms bonded to the coordinating atom
        if neighbor.element == 'H': # Verifies if it is an H atom and removes it if so
            ligand.atoms.remove(neighbor)
            break

def tag_ligand_atoms (ligand, atom, let_H = False):
    H_neighbors = atom.get_invariant()[-1]
    atom.add_tag('key', 'new', str(H_neighbors))
    '''
    for atom_old in atoms_old_ligand:
        if atom_old.element == atom.element and ('non_bonded_key' in atom_old.tags):
            atom.add_tag('non_bonded_key')
    if not('non_bonded_key' in atom.tags):
    '''
    if not(let_H):
        monodeprotonate_donor_atom(ligand, atom)

def substitute_ligand (catalyst, new_ligand, old_key_atoms, remove_H = False):

    catalyst.map_ligand(new_ligand, old_key_atoms, minimize = False )
    if remove_H:
        for atom in catalyst.atoms:
            if 'new' in atom.tags:
                H_count = 0
                is_Pd = False
                #print(atom.connected)
                for neighbor_atom in atom.connected:
                    if neighbor_atom.element == 'Pd':
                        is_Pd = True
                    if neighbor_atom.element == 'H':
                        H_count += 1
                        H_atom = neighbor_atom
                        
                if is_Pd:
                    if str(H_count) in atom.tags:
                        catalyst.atoms.remove(H_atom)
                        break

def attach_ligand (SMILES, N_atom_index):
    """
    Substitutes the ligand in catalysts from the Fujiwara paper with a new ligand.

    Given the SMILES of a ligand, this function replaces the original ligand 
    in the catalyst structure with the new ligand and writes the modified geometry 
    to an .xyz file.

    Parameters:
        SMILES (str): SMILES string of the ligand to be added.

    Output:
        Writes a new .xyz file with the geometry of the modified catalyst.
    """
    original_geometry_directory = Path(r'/home/octavian/original_geometries') # Path to the direcotry with original geometries
    new_geometry_directory = Path(r'/data/octavian/', f'{SMILES}') # Path to the directory with new geometries
    if not(new_geometry_directory.exists()):
        new_geometry_directory.mkdir(exist_ok = True)
    for j, geometry_file in enumerate(original_geometry_directory.iterdir()): # Iterates through the geometry files
        
        original_initial_cat = Geometry(f'{geometry_file}') # Gets catalyst from the paper

        filename = geometry_file.name
        type_of_compound = geometry_file.stem.strip()
        
        # Creates copies of the original catalyst
        original_cat_copy_1 = original_initial_cat.copy()
        original_cat_copy_2 = original_initial_cat.copy()
        original_cat_copy_3 = original_initial_cat.copy()

        cat_list = [original_initial_cat, original_cat_copy_1, original_cat_copy_2, original_cat_copy_3]

        for i in range(len(cat_list)):
            
            if i == 3:
                tag_catalyst_atoms(cat_list[i], type_of_compound) # Adds tags to the atoms
                atom_old_L1, atom_old_L2 = return_old_key_atoms(cat_list[i], type_of_compound, i) # Gets old keys for each type of copy

                # map_ligand modified the geometry object corresponding to
                # the new ligand; as two substitutions with the same ligand
                # are needed, this one is defined twice 
                new_ligand_1 = Geometry.from_string(SMILES)
                new_ligand_2 = Geometry.from_string(SMILES)

                # Tags are added, for the moment knowing exactly on which atoms; also, certain atoms may be deprotonated
                let_H = False if type_of_compound != 3 else True
                tag_ligand_atoms(new_ligand_1, new_ligand_1.atoms[0], let_H)
                tag_ligand_atoms(new_ligand_1, new_ligand_1.atoms[N_atom_index], let_H)
                tag_ligand_atoms(new_ligand_2, new_ligand_2.atoms[0], False)
                tag_ligand_atoms(new_ligand_2, new_ligand_2.atoms[N_atom_index], False)

                # Creates catalyst with L1 and L2 in original order
                if type_of_compound != 3:
                    substitute_ligand(cat_list[i], new_ligand_1, atom_old_L1)
                    substitute_ligand(cat_list[i], new_ligand_2, atom_old_L2)
                    #cat_list[i].display()
                else:
                    substitute_ligand(cat_list[i], new_ligand_1, atom_old_L1, True)
                    substitute_ligand(cat_list[i], new_ligand_2, atom_old_L2)
                    #cat_list[i].display()

                # Saves the geometries in new files; 0 for first variant of catalysts, 1 for second, 2 for third, 3 for fourth
                new_geometry_directory_type = new_geometry_directory / str(i)
                new_geometry_directory_type.mkdir(exist_ok = True)
                outfile = new_geometry_directory_type / filename 
                cat_list[i].write(outfile = outfile)

def replace_cpcm_with_ddcosmo(input_file):
    # Open the input file for reading
    with open(input_file, 'r') as file:
        content = file.read()
    
    # Replace 'CPCM' with 'ALPB'
    content = content.replace("CPCM", "DDCOSMO")
    
    # Open the output file for writing
    with open(input_file, 'w') as file:
        file.write(content)

def extract_K_indices(filename):
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

def automated_files_XTB (SMILES):
    smiles_path = Path(r'/data/octavian/',f'{SMILES}') # Goes to the direcotry containing the initial geometry
    for variant in smiles_path.iterdir(): # Iterates through the variants 0-3
        if str(variant.stem) == '3':
            files = [f for f in variant.iterdir() if f.is_file()]
            for file in files: # Iterates through each intermediate 0, 1, 2, 3 and ts
                #Creates a directory for the compound for which calculations are performed if the file does not already exist
                computation_directory_path = Path(variant, f'{file.stem}_comput')
                if not(computation_directory_path.exists()):
                    computation_directory_path.mkdir(exist_ok= True)
                catalyst = Geometry(f'{file}')
                index = file.stem.strip() # Type of compound being analyzed
                # Define theory level for relaxed geometry optimization
        
                method = 'XTB2'
                solvent = ImplicitSolvent('CPCM', 'benzene')
                if index != 'ts':

                    # Fixes the coordinating atoms and the Pd atom for the conformational search
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
                    job_list = [OptimizationJob(), FrequencyJob(temperature = 373)]
                    
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
                    replace_cpcm_with_ddcosmo(f'{outfile_input}')

                    #Launch ORCA
                    with open(outfile_output, 'w') as out_f:
                        orca_run = subprocess.run(f"/data/octavian/software/orca_6_1_0_linux_x86-64_shared_openmpi418/orca catalyst.inp", cwd = computation_directory_path_lowest_conf, 
                                                shell=True, stdout=out_f, stderr=subprocess.STDOUT)
                        if orca_run.returncode != 0:
                            raise RuntimeError(f"ORCA failed with return code {orca_run.returncode}")
                        
                else: # Handles the transition state
                    fixed_atoms = extract_K_indices(f'{file}') # Extracts the key tags from the geometry files
                    for i in range(1,9): # Apart from the key atoms, the Pd, the benzene C atoms and the H atom transferred need fixed
                        if not(str(i) in fixed_atoms): 
                            fixed_atoms.append(str(i))
                    ts_constrained_opt_theory_level = Theory(
                        method = method,
                        solvent = solvent,
                        job_type = [OptimizationJob(constraints = {"atoms": fixed_atoms}), FrequencyJob(temperature = 383)]
                    )
                    ORCA_BLOCKS_constrained = {"geom": ["MaxIter 500"]}
                    # Create ORCA input file for the constrained opt of the TS and run the optimization
                    constrained_outfile_input = Path(computation_directory_path, "ts_constrained.inp")
                    constrained_outfile_output = Path(computation_directory_path, "ts_constrained.out")
                    catalyst.write(
                        outfile = f'{constrained_outfile_input}',
                        theory = ts_constrained_opt_theory_level,
                        processors = 8,
                        memory = 12,
                        ORCA_BLOCKS = ORCA_BLOCKS_constrained
                    )
                    replace_cpcm_with_ddcosmo(f'{constrained_outfile_input}')

                    #Launch ORCA
                    orca_run = subprocess.run(f"/data/octavian/software/orca_6_1_0_linux_x86-64_shared_openmpi418/orca ts_constrained.inp", cwd = computation_directory_path,
                                            shell=True, stdout=open(constrained_outfile_output, 'w'), stderr=subprocess.STDOUT)
                    if orca_run.returncode != 0:
                        raise RuntimeError(f"ORCA failed with return code {orca_run.returncode}")
                   
def extract_energies(SMILES):

    SMILES_dir = Path(r'/data/octavian/',f'{SMILES}') # Path to the directory that contains the energy calculations for the given SMILES
    output_file = '/data/octavian/saturn.log'

    with open (output_file, 'a') as f:

        compound_0_file = Path(SMILES_dir, '3/0_comput/energy/catalyst.out')
        compound_ts_file = Path(SMILES_dir, '3/ts_comput/ts_constrained.out')
        compound_0_geometry = Path(SMILES_dir, '3/0_comput/energy/catalyst.xyz')
        compound_ts_geometry = Path(SMILES_dir, '3/ts_comput/ts_constrained.xyz')

        compound_0_output_content = FileReader(f'{compound_0_file}', just_geom = False)
        compound_ts_output_content = FileReader (f'{compound_ts_file}', just_geom = False)

        if not(compound_0_output_content['finished']):
            f.write(f'{SMILES} Compound 0 did not finish the optimization \n')
            f.write(f'----------------------------- \n')
        elif compound_0_output_content['error']:
            f.write(f'{SMILES} Calculation of compound 0 encountered an error \n')
            f.write(f'----------------------------- \n')
        elif not(compound_ts_output_content['finished']):
            f.write(f'{SMILES} Compound ts did not finish the optimization \n')
            f.write(f'-----------------------------\n')
        elif compound_ts_output_content['error']:
            f.write(f'{SMILES} Calculation of compound ts encountered an error \n')
            f.write(f'-----------------------------\n')
        else:

            free_energy_0 = compound_0_output_content['free_energy']
            free_energy_ts = compound_ts_output_content['free_energy']

            activation_energy = (free_energy_ts - free_energy_0 - benzene_free_energy) * conversion_eh_kcal
            f.write(f'Activation energy for the ligand {SMILES} = {activation_energy} kcal/mol \n')
            f.write(f'----------------------------- \n')

        with open (f'{compound_0_geometry}', 'r') as h:
            for line in h:
                f.write(line)
        f.write(f'----------------------------- \n')

        with open (f'{compound_ts_geometry}', 'r') as g:
            for line in g:
                f.write(line)  # prevent extra blank lines
        f.write(f'----------------------------- \n')
        f.write(f'----------------------------- \n')

    shutil.rmtree(f'{SMILES_dir}')
    return (activation_energy)    

def rearrange_smiles (smiles):

    mol = pybel.readstring('smi', smiles) # Reads smiles
    
    try:
        pattern = pybel.Smarts('[OH]c1ncccc1') # Pattern of 2-hydroxypyridine that is recognizable to the pipeline

        # Pattern findall returns a list of tuples
        # In the tuple, the indices of the atoms that form the pattern are found
        # Search for the 2-hydroxypyridine pattern and save the index of the O atom in hydroxy
        pattern_idx = pattern.findall(mol)[0][0] 

        rearranger = openbabel.OBConversion() # Defines the rearranger
        rearranger.SetInAndOutFormats('smi', 'smi') # Input and output in SMILES

        # 'f' ensures that first atom is the one wanted, in this case O
        rearranger.AddOption('f', openbabel.OBConversion.OUTOPTIONS, str(pattern_idx)) 

        outmol = openbabel.OBMol() # Defines an output molecule

        rearranger.ReadString(outmol, smiles) # Rearranges the SMILES into the output molecule
        arranged_smiles = rearranger.WriteString(outmol).strip()

        return(arranged_smiles) # Returns the rearranged SMILES
    except:
        return 0

def get_N_index (smiles):

    # The in AaronTools, the indexing in the molecule is the same as the order in the SMILES

    mol_geometry = Geometry.from_string(smiles)

    # Since the molecule starts with the O atom, the second one (index 1) is necessarily the C bonded to N
    for neighbour in mol_geometry.atoms[1].connected:
        if neighbour.element == 'N':
            for i, atom in enumerate(mol_geometry.atoms):
                if atom == neighbour:
                    return(i) # Returns the index of the N atom

#@click.group()
def cli():
    """A script with multiple functions accessible via command line."""
    pass
#@cli.command()
#@click.argument('smiles')
def calculate_energy(smiles):
    arranged_smiles = rearrange_smiles(smiles)
    N_atom_index = get_N_index(arranged_smiles)
    attach_ligand(arranged_smiles, N_atom_index)
    automated_files_XTB(arranged_smiles)
    activation_energy = extract_energies(arranged_smiles)
    return(activation_energy)

#if __name__ == "__main__":
    #cli()