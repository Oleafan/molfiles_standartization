# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 12:51:44 2023

@author: Oleafan
"""
from  indigo import Indigo
indigo = Indigo()
indigo.setOption('molfile-saving-mode', 3000)

def _get_num_arom_bonds(mol):
    num_arom = 0
    for bond in mol.iterateBonds():
        if bond.bondOrder() == 4:
            num_arom += 1
    return num_arom

def _get_relevant_aromatic_molfile(mol, max_tau = 1000):
    #find tautomer with maximum aromatic bonds
    #max_tau - breaks iterator if there are too much tautomers
    mol.foldHydrogens()
    tau_aromatic_dict = {}
    
    for idx, tautomer in enumerate(indigo.iterateTautomers(mol, 'RSMARTS')):
        
        mol_tau = tautomer.clone()
        mol_tau.aromatize()
        tau_aromatic_dict.update({mol_tau: _get_num_arom_bonds(mol_tau)})
        if idx > max_tau:
            break
    tau = max(tau_aromatic_dict, key = tau_aromatic_dict.get)
    tau.dearomatize() 

    return tau.molfile()

MOLFILE_WRONG_AMIDE = '''
  -INDIGO-04232323302D

  0  0  0  0  0  0  0  0  0  0  0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 28.125 -10.925 0.0 0
M  V30 2 N 28.991 -11.425 0.0 0
M  V30 3 O 28.125 -9.92501 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 3
M  V30 2 2 2 1
M  V30 END BOND
M  V30 END CTAB
M  END'''

query_wrong_amide = indigo.loadQueryMolecule(MOLFILE_WRONG_AMIDE) 
query_wrong_amide.aromatize()


MOLFILE_WRONG_SULFAMIDE = """
  -INDIGO-08152311142D

  0  0  0  0  0  0  0  0  0  0  0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 16.8225 -9.08611 0.0 0
M  V30 2 N 17.6886 -9.58612 0.0 0
M  V30 3 S 16.8225 -8.0861 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 3
M  V30 2 2 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
"""

query_wrong_sulfamide = indigo.loadQueryMolecule(MOLFILE_WRONG_SULFAMIDE) 
query_wrong_sulfamide.aromatize()

def standartize_molecule(mol_id, max_tau = 1000):
    #mol_id - any varian of mol representation (inchi, smiles, molfile, etc.)
    #max_tau - breaks iterator if there are too much tautomers
    #crates correct tautomer 
    try:
        mol = indigo.loadMolecule(mol_id)
        if mol is None:
            return None
    except:
        return None
    
    try:
        mol = indigo.loadMolecule(_get_relevant_aromatic_molfile(mol)) 
        #this is strange, but transfer of the mol object between functions normally works only via molfile 
    except:
        pass
    
    try:
        mol.foldHydrogens()
        mol.layout()

        matcher = indigo.substructureMatcher(mol) 
        num_matches_amide = matcher.countMatches(query_wrong_amide)
        num_matches_sulfamide = matcher.countMatches(query_wrong_sulfamide)

        if num_matches_amide > 0 or num_matches_sulfamide > 0:
            tau_dict = {} #{tau_mol: num_matches} 
            
            for idx, tautomer in enumerate(indigo.iterateTautomers(mol, 'RSMARTS')):
                mol_tau = tautomer.clone()
                mol_tau.foldHydrogens()
                mol_tau.unfoldHydrogens()
                
                matcher = indigo.substructureMatcher(mol_tau) 
                num_matches_amide = matcher.countMatches(query_wrong_amide)
                num_matches_sulfamide = matcher.countMatches(query_wrong_sulfamide)    
                tau_dict.update({mol_tau: num_matches_amide + num_matches_sulfamide})
                if idx > max_tau:
                    break
                
            tau = min(tau_dict, key = tau_dict.get) #минимизация количества неправильных амидов и сульфамидов
            tau.foldHydrogens()
            #tau.unfoldHydrogens()
            tau.layout() 
            return tau.molfile()     
                
        else:
            mol.foldHydrogens()
            #mol.unfoldHydrogens()
            mol.layout() 
            return mol.molfile()
    except:
        try:
            mol.foldHydrogens()
            #mol.unfoldHydrogens()
            mol.layout()         
            return mol.molfile()
        except:
            return None
        
        
if __name__ == '__main__':
    inchi_list = [
     'InChI=1S/C8H9NO2/c1-6(10)9-7-2-4-8(11)5-3-7/h2-5,11H,1H3,(H,9,10)',
     'InChI=1S/C7H7NO2/c9-7(8-10)6-4-2-1-3-5-6/h1-5,10H,(H,8,9)',
     'InChI=1S/C8H4ClNO3/c9-6-3-1-2-5(4-6)7-10-13-8(11)12-7/h1-4H', #remains intact
     'InChI=1S/C8H9NO3/c1-12-9-8(11)6-4-2-3-5-7(6)10/h2-5,10H,1H3,(H,9,11)',
     'InChI=1S/C8H7ClN2O4/c1-15-10-8(12)6-3-2-5(11(13)14)4-7(6)9/h2-4H,1H3,(H,10,12)',
     'InChI=1S/C9H8BrNO2/c1-12-11-9-8-3-2-7(10)4-6(8)5-13-9/h2-4H,5H2,1H3/b11-9-', #remains intact
     'InChI=1S/C12H15NO4/c1-12(2,3)16-11(15)17-13-10(14)9-7-5-4-6-8-9/h4-8H,1-3H3,(H,13,14)',
     'InChI=1S/C16H13NO/c1-11-13-9-5-6-10-14(13)16(18)17-15(11)12-7-3-2-4-8-12/h2-10H,1H3,(H,17,18)',#remains intact
     'InChI=1S/C5H5NO/c7-5-3-1-2-4-6-5/h1-4H,(H,6,7)', #remains intact
     'InChI=1S/C16H15NO/c1-11-13-9-5-6-10-14(13)16(18)17-15(11)12-7-3-2-4-8-12/h2-11,15H,1H3,(H,17,18)',
     'InChI=1S/C9H7Cl3N2O2/c10-9(11,12)8(13)16-14-7(15)6-4-2-1-3-5-6/h1-5,13H,(H,14,15)',
    'InChI=1S/C8H13NO2/c1-5-6-9-7(10)11-8(2,3)4/h1H,6H2,2-4H3,(H,9,10)',
    'InChI=1S/C8H9NO2/c1-6(10)9-7-2-4-8(11)5-3-7/h2-5,11H,1H3,(H,9,10)'] 
    for inchi in inchi_list:
        print(standartize_molecule(inchi))