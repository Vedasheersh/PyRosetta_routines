Inputs:-
  1. PDB file containing an antigen-antibody complex
  2. Epitope residues on the antigen specified along with chain as per PDB numbering;
      For eg:- EPITOPES = [('A',12),('A',14)] # where A is the antigen chain
  3. Antigen and Antibody chain names
      For eg:-  ANTIGEN = 'A' # chain A corresponds to antigen
                ANTIBODY = 'HK' # both chains H and K correspond to antibody
  4. Number of independent simulations to run (~10k is preferred)
  
 Outputs:-
  1. Output pdbs numbered according to job number
  2. Output score file containing all rosetta score-terms (also interaction energy)
  3. Log file containing designed sequences vs (total_energy , interaction_energy) in a csv format
  
 Protocol description:- 
  1. Input PDB will be minimized and docked according to J. Gray et al. 2003 to calculate WT- total_energy and interaction_energy
  2. MonteCarlo annealing based RosettaDesign is carried out for all residues within contact distance from Epitope residues
  3. Designed sequence is re-docked to calculate total_energy and interaction_energy
  
  
