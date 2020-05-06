from labutil.plugins.lammps import lammps_run, get_lammps_energy,parse_structure_dump
from labutil.objects import Struc, Dir, ClassicalPotential, ase2struc
from ase.spacegroup import crystal
from ase.build import make_supercell
import numpy, os
import matplotlib.pyplot as plt

from ase.io import write
from ase.build import fcc110
import ase
from ase.build import stack
import math

input_template = """
# ---------- 1. Initialize simulation ---------------------
units metal
atom_style atomic
dimension  3
boundary   p p p
read_data $DATAINPUT

# ---------- 2. Specify interatomic potential ---------------------
pair_style eam
pair_coeff * * $POTENTIAL

#pair_style lj/cut 4.5
#pair_coeff 1 1 0.392 2.620 4.5

# ---------- 3. Run single point calculation  ---------------------
thermo_style custom step pe lx ly lz press pxx pyy pzz
run 0

#--include optimization of the unit cell parameter
fix 1 all box/relax iso 0.0 vmax 0.001

#---enable optimization of atomic positions (and the cell)
min_style cg
minimize 1e-10 1e-10 1000 10000

# ---- 4. Define and print useful variables -------------
variable natoms equal "count(all)"
variable totenergy equal "pe"
variable length equal "lx"

print "Total energy (eV) = ${totenergy}"
print "Number of atoms = ${natoms}"
print "Lattice constant (Angstoms) = ${length}"
        """

'''
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


#for calculating cohesive energy
def compute_energy(alat, template):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potpath = os.path.join(os.environ['LAMMPS_POTENTIALS'],'Au_u3.eam')
    potential = ClassicalPotential(path=potpath, ptype='eam', element=["Au"])
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "final_project/bulk_cohesive_energy/test10_relaxed_cubic_a_and_b/", str(alat)))
    struc = make_struc(alat=alat)
    output_file = lammps_run(struc=struc, runpath=runpath, potential=potential, intemplate=template, inparam={})
    energy, lattice = get_lammps_energy(outfile=output_file)
    return energy, lattice
'''


'''
#for generating the slab of the perfect FCC Au(110) 1*1 ideal surface
def make_struc(i):   
    slab_ideal=fcc110('Au',size=(4,4,i), a=4.0800000183054, vacuum=20)
    #slab_ideal=fcc110('Au',size=(4,8,i), vacuum=20)
    write('slab_ideal.struct',slab_ideal)
    slab_ideal.edit()
    #natoms_100=(len(slab_100.numbers))
    print ('positions',slab_ideal.positions)
    natoms_ideal=len(slab_ideal.positions)
    #print(natoms_ideal)
    struc_ideal=Struc(ase2struc(slab_ideal))
    return natoms_ideal,struc_ideal
'''




#for generating the slab of the FCC Au(110) 1*2 missing row reconstruction
def make_struc(i):   
    slab_ideal=fcc110('Au',size=(16,16,40), a=4.0800000183054, vacuum=i)
    #slab_ideal=fcc110('Au',size=(4,8,i), vacuum=20)
    write('slab_ideal.struct',slab_ideal)
    #slab_ideal.edit()
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
    #slab_ideal.edit()
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



def compute_energy(template, i):
    potpath = os.path.join(os.environ['LAMMPS_POTENTIALS'],'Au_u3.eam')
    potential = ClassicalPotential(path=potpath, ptype='eam', element=["Au"])
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "final_project/Au110_reconstructed/test6_cubic/", str(i)))
    natoms,struc = make_struc(i=i)
    output_file = lammps_run(struc=struc, runpath=runpath, potential=potential, intemplate=template, inparam={})
    energy, lattice = get_lammps_energy(outfile=output_file)
    print('lattice is', lattice)

    surface_area=1/(math.sqrt(2))*lattice*lattice
    
    print(natoms)
    
    energy_per_atom=energy/natoms
    surface_energy=[]
    surface_energy=((1/(2*surface_area))*(energy-(natoms)*(-3.930000000176775)))

    return energy, lattice, energy_per_atom, surface_energy


'''
#for calculating coheisve energy
def lattice_scan():
    alat_list = numpy.linspace(4.06,4.10,21)
    #alat_list=[4.57,4.59]
    energy_list = [compute_energy(alat=a, template=input_template)[0] for a in alat_list]
    print(energy_list)
    plt.plot(alat_list, energy_list)
    plt.xlabel('lattice constant (Angstrom)')
    plt.ylabel('energy (eV)')
    plt.show()
'''


def lattice_scan():
    i_list=numpy.linspace(3,20,18)
    print(i_list)
    surface_energy_list=[compute_energy(template=input_template,i=int(a))[3] for a in i_list]
    print('surface energy is ',surface_energy_list)
    plt.plot(i_list,surface_energy_list)
    plt.xlabel('length (x,y dimension) of the slab')
    plt.ylabel('surface energy (eV/A^2)')
    plt.show
 
'''    
surface_energy_x_and_y=[0.06042452351867576, 0.06010373941484659, 0.060300668033233716, 0.06010553429149964, 
                        0.06026721967220595, 0.060128318222880194, 0.060235136817377105, 0.06014707142368978, 
                        0.06024130703549732, 0.060159461141438296, 0.060256525354933206, 0.060117725747536925,
                        0.060206515150248954, 0.060164755985760966, 0.06019450599291324, 0.0601311175828273, 
                        0.06022481022835219, 0.06018618330366263]


surface_energy_versus_vacuum_size= [0.06016475267199812, 0.06016476928707486, 0.06016476323861836, 
                                      0.06016477638651081, 0.06016477153085311, 0.06016477907147559, 
                                      0.0601647627401785, 0.060164786033552665, 0.060164752409000564, 
                                      0.06016476431691229, 0.06016478361129873, 0.06016477382892909, 
                                      0.06016476090477189, 0.060164783186645567, 0.060164783031863074,
                                      0.06016478034620985, 0.0601647805824661, 0.060164755985760966,
                                      0.060164794514424885, 0.06016476083093836, 0.06016477019223293, 
                                      0.06016475321878943, 0.06016476079161699, 0.06016476229719527,
                                      0.06016477761505024, 0.060164784349495445, 0.06016477466809168, 
                                      0.06016478174059472, 0.06016477663727291, 0.0601647746085878, 
                                      0.06016476792419066, 0.06016478191025574, 0.06016478235315225, 
                                      0.06016475893473109, 0.06016475174489961, 0.06016475024888494,
                                      0.060164793804236644, 0.060164768305967155]
'''




if __name__ == '__main__':
    # put here the function that you actually want to run
    #make_struc(4.123,3)
    #compute_energy(4.082,input_template,1)
    #compute_energy(4.082,input_template,2)
    #compute_energy(4.082,input_template,3)
    #compute_energy(input_template,2)
    #lattice_scan(4.1)
    
   #make_struc(2)
   #compute_energy(2.950,input_template)
    lattice_scan()