import os

import json

from pyrosetta import *
init("-ignore_zero_occupancy false -use_input_sc -ex1 -ex2 -ex2aro")

from pyrosetta import Pose, Vector1, pose_from_file, create_score_function, PyJobDistributor
from pyrosetta.rosetta import core, protocols
from pyrosetta import MonteCarlo
from pyrosetta.rosetta.protocols.rigid import RigidBodyPerturbMover

from pyrosetta.rosetta.core.simple_metrics.metrics import InteractionEnergyMetric

from pyrosetta.rosetta.protocols import minimization_packing as pack_min
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation

from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.select.residue_selector import NotResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.core.select.residue_selector import InterGroupInterfaceByVectorSelector
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector

from pyrosetta.rosetta.utility import vector1_std_string

class AffMut(object):
    """docstring for ReDocker"""
    def __init__(self, PDB,EPITOPES = [], ANTIBODY = 'HK', ANTIGEN = 'A',pymol=True):
        pose = Pose() 
        self.pose = pose_from_file(PDB)
        self.native_pose = self.pose.clone()
        self.scorefxn = create_score_function('ref2015')
        self.partners = '{0}_{1}'.format(ANTIGEN,ANTIBODY)
        self.antibody = ANTIBODY
        self.antigen = ANTIGEN
        
        self.epitopes_pdb = EPITOPES
        # store corresponding pose residue numbers
        self.epitopes_pose = []
        
        info = self.pose.pdb_info()
        for ep in EPITOPES:
            chain,res = ep
            self.epitopes_pose.append(info.pdb2pose(chain,res))
        
        if pymol:
            pymover = PyMOLMover() 
            self.pymol_mover = pymover
            pymover.update_interval(1)
            pymover.keep_history(True)
            self.send_to_pymol()
              
    def redock(self,hires=True,perturb=False,TRANS_PERT=0.05,ROT_PERT=1):
        # setup the docking FoldTree
        # using this method, the jump number 1 is automatically set to be the
        #    inter-body jump
        dock_jump = 1
        protocols.docking.setup_foldtree(self.pose, self.partners, Vector1([dock_jump]))
        
        # create ScoreFunctions for centroid and fullatom docking
        self.scorefxn = create_score_function('ref2015')
        
        # setup the high resolution (fullatom) docking protocol (DockMCMProtocol)

        hiresdocker = protocols.docking.DockMCMProtocol()
        hiresdocker.set_scorefxn(self.scorefxn)
        hiresdocker.set_first_cycle(4)
        hiresdocker.set_second_cycle(20)
        self.hiresdocker = hiresdocker
        
        # Setup perturber
        dockperturb = RigidBodyPerturbMover(dock_jump,TRANS_PERT,ROT_PERT)
        self.perturber = dockperturb
        
        if hires:
            self.hiresdocker.apply(self.pose)
        elif perturb:
            self.perturber.apply(self.pose)
    
    def relax(self,ncycles=5):
        relaxer = pyrosetta.rosetta.protocols.relax.FastRelax(self.scorefxn,ncycles)
        relaxer.apply(self.pose)
    
    def minimize(self):
        min_mover = pack_min.MinMover()
        mm = MoveMap()
        mm.set_bb(True)
        mm.set_chi(True)
        mm.set_jump(True)
        min_mover.movemap(mm)
        min_mover.score_function(self.scorefxn)
        min_mover.apply(self.pose)
        
    def send_to_pymol(self): self.pymol_mover.apply(self.pose)
        
    def calc_energy(self):
        
        min_mover = pack_min.MinMover()
        mm = MoveMap()
        mm.set_bb(True)
        mm.set_chi(True)
        mm.set_jump(True)
        min_mover.movemap(mm)
        min_mover.score_function(self.scorefxn)
        min_mover.apply(self.pose)
        
        return self.pose.scores['total_score']
    
    def calc_interaction(self):
        ie = InteractionEnergyMetric()
        
        c1 = ChainSelector()
        c2 = ChainSelector()
        vec1 = vector1_std_string()
        vec2 = vector1_std_string()
        for chain in self.antibody:
            vec1.append(chain)
        for chain in self.antigen:
            vec2.append(chain)

        c2.set_chain_strings(vec2)
        c1.set_chain_strings(vec1)
        ie.set_residue_selectors(c1,c2)
        
        ie.apply(self.pose)

        return self.pose.scores['interaction_energy']
        
    def design(self,antibody=True,antigen=False,pack_only=False):
        self.pose = self.native_pose.clone()
        info = self.pose.pdb_info()

        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())
        
        epi_res = ''
        for res in self.epitopes_pose:
            epi_res+='{0},'.format(res)
        epi_selector = ResidueIndexSelector(epi_res)
        
        antibody_selector = ChainSelector()
        vec = vector1_std_string()
        for chain in self.antibody:
            vec.append(chain)
            vec.append(chain)
        antibody_selector.set_chain_strings(vec)
        interface_res_selector = InterGroupInterfaceByVectorSelector(epi_selector,antibody_selector)
        interface_antibody_selector = AndResidueSelector(interface_res_selector,antibody_selector)

        if pack_only: tf.push_back(operation.RestrictToRepacking())
        
        else:
            if antigen: design_selector = epi_selector
        
            elif antibody: design_selector = interface_antibody_selector
        
            # FIRST select designable residues and prevent everything else from design and packing
            no_design_selector = NotResidueSelector(design_selector)
            prevent_repacking_rlt = operation.PreventRepackingRLT()
            #restrict_topack = operation.RestrictToRepackingRLT()
            prevent_design = operation.OperateOnResidueSubset(prevent_repacking_rlt,no_design_selector,False)
            tf.push_back(prevent_design)
            
        print(tf.create_task_and_apply_taskoperations(self.pose))
        
        packer = pack_min.PackRotamersMover('ref2015')
        packer.task_factory(tf)

        packer.apply(self.pose)

def protocol(input_pdb,epitopes,antibody,antigen,njobs=100,outname='output',logfile='log',design_antigen=False):
    
    fscore = open('{0}.sc'.format(outname),'w')
    flog = open(logfile,'w')
    flog.write('SEQUENCE,ENERGY,INTERACTION\n')
    
    affmut = AffMut(PDB=input_pdb,EPITOPES=epitopes,ANTIBODY=antibody,ANTIGEN=antigen)
    
    scores = {}
    energies = {'wt_total':0,'wt_inter':0,'mut_total':[],'mut_inter':[]}

    affmut.minimize()
    energies['wt_total'] = affmut.calc_energy()
    affmut.redock()
    energies['wt_inter'] = affmut.calc_interaction()
    flog.write('{0},{1},{2}\n'.format(affmut.pose.sequence(),energies['wt_total'],energies['wt_inter']))
    scores['WT'] = affmut.pose.scores
    
    for job in range(njobs):  
        if design_antigen: affmut.design(antibody=False,antigen=True)
        else: affmut.design()
        affmut.redock()
        ener = affmut.calc_energy()
        energies['mut_total'].append(ener)
        inter = affmut.calc_interaction()
        energies['mut_inter'].append(inter)
        flog.write('{0},{1},{2}\n'.format(affmut.pose.sequence(),ener,inter))
        
        affmut.pose.pdb_info().name(outname + '_' + str(job))
        affmut.pose.dump_pdb(outname + '_' + str(job)+'.pdb')
        scores[outname + '_' + str(job)] = affmut.pose.scores
        affmut.pose.clear()
    
    flog.close()
    fscore.write(json.dumps(scores))
    fscore.close()
    
    return affmut, energies

if __name__ == '__main__':
    # run 100 times for 10,000 structures
    for j in range(1,100):
    	INPUT_PDB = 'WT1.pdb'

    	EPITOPES = [('A',346),('A',347),('A',348),('A',351),('A',352),('A',354),('A',355)]

    	ANTIGEN = 'A'
    	ANTIBODY = 'HK'

    	affmut,energies = protocol(INPUT_PDB,EPITOPES,ANTIBODY,ANTIGEN,njobs=100,\
    	outname='WT1_Set{}'.format(j),logfile='WT1_Set{}.log'.format(j),design_antigen=False)
