"""
Code for the physical model to asses the contribution of longitudinal stiffness
to the furrow formation of Drosophila melanogaster embryos during gastrulation

Part of the paper"
"Embryonic cell populations display highly dynamic mechanical transitions during gastrulation"
Juan Manuel Gomez, Carlo Bevilacqua, Abhisha Thayambath, Maria Leptin, Julio M Belmonte and Robert Prevedel


"Code written by Abhisha Thayambath and Julio Monti Belmonte
Department of Physics, North Carolina State University

[Modified on 06/15/2024]

"""

global cd, bR, L, T, nOrder, Time
global LamV, tV

"""
Parameters varied in the paper:
    lambda_link_basal             (line 56)
    lambda_link_apical_peripheral (line 62)
    lambda_link_apical_central    (line 63)
    n_cells                       (line 85)
    starting_stiffness            (line 70)
"""

# PARAMETERS:
cd = 12                 # typical cell diameter
bR = 6                  # radius (in cd)
L = 2 * (bR + 6) * cd   # size of lattice (LxL)
tot = 47                # total number of cells
#
# POTTS PARAMETERS:
T = 25              # temperature
nOrder = 2          # distance of interaction
#
# VOLUME/SURFACE PARAMETERS:
LamV = 2.0                      #strength of volume constraint for neuroectodermal cells
tV = cd * cd                    #target volume
LamV_bas=15.0                   #strength of volume constraint for basal compartment of mesodermal cells
LamV_core=10.0                   #strength of volume constraint for core compartment of mesodermal cells
LamV_ap=2.0                     #strength of volume constraint for apical compartment of mesodermal cells
lambda_surface_apical = 1       #strength of surface constraint for apical compartment of mesodermal cells
#
# APICAL CONSTRICTION PARAMETERS
lam_d = 1.0                     #minimum value of constriction force
peak_val_multiplier = 30        #multiplier for the myosin profile
target_distance = 1.0           #target distance between neighboring apical compartments' COMs
t_relax = 2000                  #relaxation time after each myosin profile is added
Time = 9*t_relax                #Total time 
#
# CELL STIFFNESS PARAMETERS
# SUB-BASAL STIFFNESS 	
lambda_link_basal = 50          #strength of basal-core links(sub-basal stiffness)

lambda_links_basal_central = [lambda_link_basal]*9      #define sub-basal stiffness for each myosin profile
lambda_links_basal_peripheral = [lambda_link_basal]*9    

# SUB-APICAL STIFFNESS
lambda_link_apical_peripheral = 25                    #strength of peripheral apical-core links(peripheral sub-apical stiffness)
lambda_link_apical_central = 200                      #strength of central apical-core links(central sub-apical stiffness)

#constant(not varying in time) stiffness
#lambda_links_apical_central = [lambda_link_apical_central]*9         
#lambda_links_apical_peripheral = [lambda_link_apical_peripheral]*9

#dynamic stiffness
starting_stiffness = 30         #starting strength of apical-core links(starting sub-apical stiffness)

step_peri = (lambda_link_apical_peripheral - starting_stiffness)/2.0     #define step-size for dynamic stiffness
step_cent = (lambda_link_apical_central - starting_stiffness)/3.0

lambda_links_apical_central = [starting_stiffness,starting_stiffness,starting_stiffness,starting_stiffness,
    starting_stiffness+step_cent,starting_stiffness+step_cent,
    starting_stiffness+(2*step_cent),starting_stiffness+(2*step_cent),lambda_link_apical_central]
    
lambda_links_apical_peripheral = [starting_stiffness,starting_stiffness,starting_stiffness+step_peri,
        lambda_link_apical_peripheral,lambda_link_apical_peripheral,lambda_link_apical_peripheral,
        lambda_link_apical_peripheral,lambda_link_apical_peripheral,lambda_link_apical_peripheral]

#NUMBER OF CENTRAL CELLS
#n_cells=4 implies (2n+1) = 9 central cells
n_cells = 4

def configure_simulation():
    from cc3d.core.XMLUtils import ElementCC3D
    
    CompuCell3DElmnt = ElementCC3D("CompuCell3D", {"Revision": "0", "Version": "4.3.2"})
    
    MetadataElmnt = CompuCell3DElmnt.ElementCC3D("Metadata")
    
    # Basic properties of simulation
    MetadataElmnt.ElementCC3D("NumberOfProcessors", {},1)
    MetadataElmnt.ElementCC3D("DebugOutputFrequency", {},0)
    PottsElmnt = CompuCell3DElmnt.ElementCC3D("Potts")
    
    # Basic properties of CPM (GGH) algorithm
    PottsElmnt.ElementCC3D("Dimensions", {"x": L*1.4, "y": L/2 + 30, "z":1})
    PottsElmnt.ElementCC3D("Steps", {}, Time)
    PottsElmnt.ElementCC3D("Temperature", {}, int(T))
    PottsElmnt.ElementCC3D("NeighborOrder", {}, nOrder)
    
    PluginElmnt = CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "CellType"})
    # Listing all cell types in the simulation
    PluginElmnt.ElementCC3D("CellType", {"TypeId": "0", "TypeName": "Medium"})
    PluginElmnt.ElementCC3D("CellType", {"TypeId": "1", "TypeName": "basal"})
    PluginElmnt.ElementCC3D("CellType", {"TypeId": "2", "TypeName": "core"})
    PluginElmnt.ElementCC3D("CellType", {"TypeId": "3", "TypeName": "apical"})
    PluginElmnt.ElementCC3D("CellType", {"TypeId": "4", "TypeName": "basal_prime"})
    PluginElmnt.ElementCC3D("CellType", {"TypeId": "5", "TypeName": "core_prime"})

    #  List of plugins used
    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "VolumeLocalFlex"})
    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "SurfaceLocalFlex"})
    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "MomentOfInertia"})
    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "ConnectivityGlobal"})    
    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "CenterOfMass"})
    CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "NeighborTracker"})   

    # External Contact Energies
    # medium
    contact = CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "Contact"})    
    contact.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "Medium"}, 0)
    contact.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "apical"}, 8)
    contact.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "core"}, 100)
    contact.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "basal"}, 8)    
    contact.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "core_prime"}, 100)
    contact.ElementCC3D("Energy", {"Type1": "Medium", "Type2": "basal_prime"}, 8)
    # apical
    contact.ElementCC3D("Energy", {"Type1": "apical", "Type2": "apical"}, 10)
    contact.ElementCC3D("Energy", {"Type1": "apical", "Type2": "core_prime"}, 100)
    contact.ElementCC3D("Energy", {"Type1": "apical", "Type2": "core"}, 100)
    contact.ElementCC3D("Energy", {"Type1": "apical", "Type2": "basal_prime"}, 100)
    contact.ElementCC3D("Energy", {"Type1": "apical", "Type2": "basal"}, 100)
    # core 
    contact.ElementCC3D("Energy", {"Type1": "core", "Type2": "core_prime"}, 10)
    contact.ElementCC3D("Energy", {"Type1": "core", "Type2": "core"}, 100) 
    contact.ElementCC3D("Energy", {"Type1": "core_prime", "Type2": "core_prime"}, 100) 
    contact.ElementCC3D("Energy", {"Type1": "core", "Type2": "basal_prime"}, 100)
    contact.ElementCC3D("Energy", {"Type1": "core_prime", "Type2": "basal_prime"}, 100)
    contact.ElementCC3D("Energy", {"Type1": "core", "Type2": "basal"}, 100)
    contact.ElementCC3D("Energy", {"Type1": "core_prime", "Type2": "basal"}, 100)
    # basal
    contact.ElementCC3D("Energy", {"Type1": "basal", "Type2": "basal_prime"}, 5)
    contact.ElementCC3D("Energy", {"Type1": "basal", "Type2": "basal"}, 100)
    contact.ElementCC3D("Energy", {"Type1": "basal_prime", "Type2": "basal_prime"}, 100)
    # -neighbor order
    contact.ElementCC3D("NeighborOrder", {}, 4)
    
    #Internal Contact Energies
    contactIn = CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "ContactInternal"})
    contactIn.ElementCC3D("Energy", {"Type1": "basal", "Type2": "core"}, 2.0)
    contactIn.ElementCC3D("Energy", {"Type1": "apical", "Type2": "basal"}, 100.0)
    contactIn.ElementCC3D("Energy", {"Type1": "apical", "Type2": "core"}, 2.0)
    contactIn.ElementCC3D("Energy", {"Type1": "basal_prime", "Type2": "core_prime"}, 2.0)
    contactIn.ElementCC3D("Energy", {"Type1": "apical", "Type2": "basal_prime"}, 100.0)
    contactIn.ElementCC3D("Energy", {"Type1": "apical", "Type2": "core_prime"}, 2.0)
    # -neighbor order
    contactIn.ElementCC3D("NeighborOrder", {}, 4)
    
    PluginElmnt_4 = CompuCell3DElmnt.ElementCC3D("Plugin", {"Name": "FocalPointPlasticity"})
    PluginElmnt_4.ElementCC3D("Local")
    
    CompuCellSetup.setSimulationXMLDescription(CompuCell3DElmnt)

from cc3d import CompuCellSetup

configure_simulation()

from Gastrulation_simplifiedSteppables import *

#INITIAL CONDITIONS
CompuCellSetup.register_steppable(steppable=Draw(10, _tV=tV, _LamV=LamV, _LamV_bas=LamV_bas, _LamV_core=LamV_core, _LamV_ap=LamV_ap, _Radius=bR * cd))

#SURFACE CONSTRAINT
CompuCellSetup.register_steppable(steppable=surface_constraint(10, _lambda_surface_apical=lambda_surface_apical))

#SUB-APICAL AND SUB-BASAL LONGITUDINAL STIFFNESS CONSTRAINT
CompuCellSetup.register_steppable(steppable=cell_internal_links(10, tot, t_relax, _lambda_links_basal_central=lambda_links_basal_central, _lambda_links_basal_peripheral=lambda_links_basal_peripheral,
_lambda_links_apical_central=lambda_links_apical_central, _lambda_links_apical_peripheral=lambda_links_apical_peripheral, _n_cells=n_cells))

#SPRING-LIKE INTERCELLULAR FORCE CONSTRAINT FOR APICAL CONSTRICTION
CompuCellSetup.register_steppable(steppable=apical_constriction(10, lam_d, peak_val_multiplier, target_distance, t_relax))

#METRIC-FURROW DEPTH
CompuCellSetup.register_steppable(steppable=metric(10))

#PIFF FILE GENERATOR
CompuCellSetup.register_steppable(steppable=piff_generator(2000, Time))

CompuCellSetup.run()
