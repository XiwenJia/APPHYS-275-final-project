import os
from labutil.plugins.pwscf import run_qe_pwscf, PWscf_inparam, parse_qe_pwscf_output
from labutil.objects import Struc, Dir, ase2struc, Kpoints, PseudoPotential
from ase.spacegroup import crystal
from ase.io import write
import matplotlib.pyplot as plt
import numpy
import time

from ase.io import write
from ase.build import fcc110
import ase
from ase.build import stack
import math


#generating the bulk (for calculating cohesive energy)
def make_struc(alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    unitcell = crystal('Au', [(0, 0, 0)], spacegroup=225, cellpar=[alat, alat, alat, 90, 90, 90])
    #unitcell.edit()
    natoms_unitcell=(len(unitcell.numbers))
    print('number of atoms in the bulk is', natoms_unitcell)
    structure = Struc(ase2struc(unitcell))
    return structure


#for generating the slab of the perfect FCC Au(110) 1*1 ideal surface
def make_struc():   
    slab_ideal=fcc110('Au',size=(4,2,8), a=4.163737636, vacuum=20)
    #slab_ideal=fcc110('Au',size=(4,8,i), vacuum=20)
    write('slab_ideal.struct',slab_ideal)
    slab_ideal.edit()
    #natoms_100=(len(slab_100.numbers))
    #print ('positions',slab_ideal.positions)
    natoms_ideal=len(slab_ideal.positions)
    #print(natoms_ideal)
    struc_ideal=Struc(ase2struc(slab_ideal))
    return natoms_ideal,struc_ideal



#for generating the slab of the FCC Au(110) 1*2 missing row reconstruction
def make_struc():   
    slab_ideal=fcc110('Au',size=(4,2,8), a=4.163737636, vacuum=20)
    #slab_ideal=fcc110('Au',size=(4,8,i), vacuum=20)
    write('slab_ideal.struct',slab_ideal)
    slab_ideal.edit()
    #natoms_100=(len(slab_100.numbers))
    list_z=[slab_ideal.positions[i][2] for i in range(len(slab_ideal.positions))]
    min_z=min(list_z)
    list_x_top_row=[]
    for i in range(len(slab_ideal.positions)):
        if slab_ideal.positions[i][2]==min_z:
            list_x_top_row.append(slab_ideal.positions[i][0])
    list_x_top_row_concise=list(set(list_x_top_row))
    list_x_top_row_concise.sort()
    select_list_x=list_x_top_row_concise[1::2]    #this slices list_x_top_row_concise using even index
    atoms_to_delete=[]
    atoms_to_delete_index=[]
    for i in range(len(slab_ideal.positions)):
        if slab_ideal.positions[i][0] in select_list_x and slab_ideal.positions[i][2] == min_z:
            atoms_to_delete.append(slab_ideal.positions[i])     
            atoms_to_delete_index.append(i)
    atoms_to_delete_index.sort()
    atoms_to_delete_index.reverse()
    #to delete atoms in slab to create the 1*2 missing row reconstruction
    for i in atoms_to_delete_index:
        slab_ideal.pop(i)
    slab_ideal.edit()
    #print('list of atoms is',slab_ideal.positions)
    #print('list_x_top_row',list_x_top_row)
    #print('list_x_top_row_concise',list_x_top_row_concise)
    #print ('min z is',min_z)
    #print('select_list_x',select_list_x)
    #print('atoms_to_delete',atoms_to_delete)
    print('atoms_to_delete_index',atoms_to_delete_index)
    natoms_ideal=len(slab_ideal.positions)
    print('number of atoms left',natoms_ideal)
    struc_ideal=Struc(ase2struc(slab_ideal))
    return natoms_ideal,struc_ideal



def compute_energy(alat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    start_time=time.time()
    potname = 'Au.pbe-n-kjpaw_psl.1.0.0.UPF'
    pseudopath = os.environ['QE_POTENTIALS']
    potpath = os.path.join(pseudopath, potname)
    pseudopots = {'Au': PseudoPotential(name=potname, path=potpath, ptype='paw', element='Au', functional='PBE')}
    struc = make_struc(alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "final_project/bulk_cohesive_energy_QE/alat_convergence/", str(alat)))
    input_params = PWscf_inparam({
        'CONTROL': {
           # 'calculation': 'scf',
            'calculation':'vc-relax',
            'pseudo_dir': pseudopath,
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none',
        },
        'SYSTEM': {
            'ecutwfc': ecut,
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.5,
            'conv_thr': 1e-7,
        },
        'IONS': {},
        'CELL': {},
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts,ncpu=4)
    output = parse_qe_pwscf_output(outfile=output_file)
    end_time=(time.time()-start_time)
    print('running time is',end_time)
    return output



def compute_energy(ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    start_time=time.time()
    potname = 'Au.pbe-n-kjpaw_psl.1.0.0.UPF'
    pseudopath = os.environ['QE_POTENTIALS']
    potpath = os.path.join(pseudopath, potname)
    pseudopots = {'Au': PseudoPotential(name=potname, path=potpath, ptype='paw', element='Au', functional='PBE')}
    struc = make_struc()[1]
    kpts = Kpoints(gridsize=[5, 7, 1], option='automatic', offset=False)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "final_project/Au110_ideal_QE/scf_v3", str(ecut)))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'scf',
          #  'calculation':'vc-relax',
            'pseudo_dir': pseudopath,
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none',
        },
        'SYSTEM': {
            'ecutwfc': ecut,
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.5,
            'conv_thr': 1e-7,
        },
        'IONS': {},
        'CELL': {},
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts,ncpu=24)
    output = parse_qe_pwscf_output(outfile=output_file)
    end_time=(time.time()-start_time)
    print('running time is',end_time)
    return output



def lattice_scan():
    
    #test convergence of energy with respect to cutoff energies
    nk=4
    alat=4.065
    ecut_list=[41,50,55,65,70,75,80,90]
    print('cutoff energy (Ry) is',ecut_list )
    #calculate the cohesive energy by dividing total energy by the number of atoms (4)
    output=[compute_energy(alat=alat,ecut=i,nk=nk) for i in ecut_list] 
    energy = [output[b]['energy'] for b in range(len(output))]
    print('total energy (eV) is',energy)
    plt.plot(ecut_list,energy)
    plt.xlabel('cutoff energy (Ry)')
    plt.ylabel('total energy (eV)')
    plt.show
   
    
    
    #test convergence of energy with respect to k-points
    ecut=50
    alat=4.065
    nk_list=[4,5,8,9,10,11,12,13,14,15,16]
    print('number of k-points is',nk_list)
    output=[compute_energy(alat=alat,ecut=ecut,nk=int(i)) for i in nk_list] 
    energy = [output[b]['energy'] for b in range(len(output))]
    print('total energy (eV) is',energy)
    plt.plot(nk_list,energy)
    plt.xlabel('number of k-points')
    plt.ylabel('total energy (eV)')
    plt.show
    
    
    #optimize lattice parameter using 'vc-relax'     
    nk = 13
    print(nk)
    ecut = 50
    alat = 4.065
    output = compute_energy(alat=alat, ecut=ecut, nk=nk)
    energy = output['energy']
    print(energy)
   

    
def lattice_scan():  
    ecut=50
    output=compute_energy(ecut=ecut)
    energy=output['energy']
    print('total energy is',energy)
    

    

if __name__ == '__main__':
    # put here the function that you actually want to run
    #make_struc()
    lattice_scan()
    
    
    
    
