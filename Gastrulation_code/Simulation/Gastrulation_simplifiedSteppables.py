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

from cc3d.core.PySteppables import *
from math import *
import numpy as np

#*******************************************************************************************************#
#INITIAL CONDITIONS
class Draw(SteppableBasePy):
    """
    Draws the initial condition of the 2D model with 47 cells(19 mesodermal cells and 
    28 neuroectodermal cells) of the ventral half of the Drosophila embryo. Each cell is composed of
    apical, core and basal compartments. Volume constraint for each of these compartments are also defined here.

    Arguments:
      _frequency -- determines how often this Steppable is executed
      _tV        -- target volume of each cell compartment
      _LamV      -- strength of volume constraint for neurectodermal cells
      _LamV_bas  -- strength of volume constraint for basal compartment of mesodermal cells
      _LamV_core -- strength of volume constraint for core compartment of mesodermal cells
      _LamV_ap   -- strength of volume constraint for apical compartment of mesodermal cells
      _Radius    -- radius of the ventral half of the Drosophila embryo in terms of cell diameter(cd)

    """
    def __init__(self,_frequency,_tV,_LamV,_LamV_bas, _LamV_core, _LamV_ap, _Radius):
        SteppableBasePy.__init__(self,_frequency)
        self.tV=_tV; self.LamV=_LamV;
        self.Radius=_Radius; self.cd=sqrt(_tV);
        self.LamV_bas=_LamV_bas; self.LamV_core=_LamV_core; self.LamV_ap=_LamV_ap;
        
    def start(self):
        tV = self.tV                            #target volume
        cd = sqrt(tV)                           #cell diameter
        Lx = self.dim.x/2                       #half-width of coretice
        Ly = self.dim.y                         #height of the coretice
        
        nLayers = int(self.Radius/cd)           #6 layers
        n = 47                                  #number of cells
        dAng = pi/n                             #angular separation between cells
        rMax = int(cd/2)+0.5                    #maximum radius for the first layer
        
        #Apical, Basal and Core cell types are created in sequential layers 
        for Layer in range(nLayers-3):
            rMin = rMax
            rMax = int(sqrt(rMin*rMin + 2*tV/dAng))
        for Layer in range(nLayers-3,nLayers):
            cellList=[]
            if (Layer==nLayers-3):  
                type=self.BASAL
                rMin=rMax
                rMax=int(sqrt(rMin*rMin + tV*0.5/dAng))
            elif (Layer==nLayers-2): 
                type=self.CORE
                rMin=rMax
                rMax=int(sqrt(rMin*rMin+tV*5.5/dAng))
            if (Layer==nLayers-1):  
                type=self.APICAL
                rMin=rMax
                rMax=int(sqrt(rMin*rMin+tV*0.5/dAng)) 
            
            for i in range(n):
                if(Layer==nLayers-3):
                    if(i%2 == 1):
                        cell = self.new_cell(self.BASAL)          #Basal, Basal Prime and Core, Core Prime
                    elif(i%2 == 0):                               #have same physical properties; two cell-types are made to
                        cell = self.new_cell(self.BASAL_PRIME)    #prevent cells from fragmenting
                elif(Layer==nLayers-2):
                    if(i%2 == 1):
                        cell = self.new_cell(self.CORE)
                    elif(i%2 == 0):
                        cell = self.new_cell(self.CORE_PRIME)
                else:
                    cell = self.new_cell(type)
                cellList.append(cell)
                
            #assigns cells to lattice sites based on angular position and calculated distance from the center for each layer
            for x in range(int(Lx-rMax),int(Lx+rMax)):
                for y in range(int(Ly-rMax),int(Ly)):
                    d=sqrt((x-Lx)**2+(y-Ly)**2)
                    if (d<rMax and d>=rMin):
                        Ang=acos((x-Lx)/d)
                        if ((y-Ly)<0): Ang=pi-Ang
                        self.cell_field[x, y, 0] = cellList[int(Ang/dAng)]     
            
            self.mid = 24               #cell-id at the ventral midline
            self.tot = 47  
            
            #Assign volume properties for each cell type
            bas_list = []
            bas_prime_list = []
            for cell in self.cellListByType(self.BASAL): #BASAL
                cell.targetVolume=cell.volume   
                if (self.mid-9<= cell.id <=self.mid+9):
                    cell.lambdaVolume=self.LamV_bas
                    bas_list.append(cell.id)
                elif (cell.id == 1) or (cell.id == self.tot):
                    cell.lambdaVolume=self.LamV*1000
                    bas_list.append(cell.id)
                elif (cell.id == 2) or (cell.id == self.tot-1):
                    cell.lambdaVolume=self.LamV
                    bas_list.append(cell.id)
                elif (2 < cell.id < self.mid-9) or (self.mid+9 < cell.id < self.tot-1):
                    cell.lambdaVolume=self.LamV*50 
                    bas_list.append(cell.id)
                    
            for cell in self.cellListByType(self.BASAL_PRIME): #BASAL_PRIME
                cell.targetVolume=cell.volume   
                if (self.mid-9<= cell.id <=self.mid+9):
                    cell.lambdaVolume=self.LamV_bas
                    bas_prime_list.append(cell.id)
                elif (cell.id == 1) or (cell.id == self.tot):
                    cell.lambdaVolume=self.LamV*1000
                    bas_prime_list.append(cell.id)
                elif (cell.id == 2) or (cell.id == self.tot-1):
                    cell.lambdaVolume=self.LamV
                    bas_prime_list.append(cell.id)
                elif (2 < cell.id < self.mid-9) or (self.mid+9 < cell.id < self.tot-1):
                    cell.lambdaVolume=self.LamV*50 
                    bas_prime_list.append(cell.id)
            
            core_list = []
            core_prime_list = []
            for cell in self.cellListByType(self.CORE): #CORE
                cell.targetVolume=cell.volume
                if ((1*self.tot)+self.mid-9<= cell.id <=(1*self.tot)+self.mid+9):
                    cell.lambdaVolume=self.LamV_core
                    core_list.append(cell.id)
                elif (cell.id == (1*self.tot)+1) or (cell.id == 2*self.tot):
                    cell.lambdaVolume=self.LamV*1000
                    core_list.append(cell.id)
                elif (cell.id == (1*self.tot)+2) or (cell.id == (2*self.tot)-1):
                    cell.lambdaVolume=self.LamV
                    core_list.append(cell.id)
                elif ((1*self.tot)+2 < cell.id < (1*self.tot)+self.mid-9) or ((1*self.tot)+self.mid+9 < cell.id < (2*self.tot)-1):
                    cell.lambdaVolume=self.LamV*50
                    core_list.append(cell.id)
                    
            for cell in self.cellListByType(self.CORE_PRIME): #CORE_PRIME
                cell.targetVolume=cell.volume
                if ((1*self.tot)+self.mid-9<= cell.id <=(1*self.tot)+self.mid+9):
                    cell.lambdaVolume=self.LamV_core
                    core_prime_list.append(cell.id)
                elif (cell.id == (1*self.tot)+1) or (cell.id == 2*self.tot):
                    cell.lambdaVolume=self.LamV*1000
                    core_prime_list.append(cell.id)
                elif (cell.id == (1*self.tot)+2) or (cell.id == (2*self.tot)-1):
                    cell.lambdaVolume=self.LamV
                    core_prime_list.append(cell.id)
                elif ((1*self.tot)+2 < cell.id < (1*self.tot)+self.mid-9) or ((1*self.tot)+self.mid+9 < cell.id < (2*self.tot)-1):
                    cell.lambdaVolume=self.LamV*50
                    core_prime_list.append(cell.id)
               
            ap_list = [] 
            for cell in self.cellListByType(self.APICAL): #APICAL
                cell.targetVolume=cell.volume + 0.5
                if ((2*self.tot)+self.mid-9<= cell.id <=(2*self.tot)+self.mid+9):
                    cell.lambdaVolume=self.LamV_ap
                    self.connectivityGlobalPlugin.setConnectivityStrength(cell,20000000)        #used to prevent breakage/fragmentation of apical compartment during constriction
                elif (cell.id == (2*self.tot)+1) or (cell.id == 3*self.tot):
                    cell.lambdaVolume=self.LamV*1000
                elif (cell.id == (2*self.tot)+2) or (cell.id == (3*self.tot)-1):
                    cell.lambdaVolume=self.LamV
                else:
                    cell.lambdaVolume=self.LamV*50
                ap_list.append(cell.id)
                
        #Group each apical, core and basal domains into corresponding cell/cluster
        for cell in self.cellListByType(self.BASAL): #BASAL
            core_cell = self.fetch_cell_by_id(cell.id+n)
            ap_cell = self.fetch_cell_by_id(cell.id+(2*n))
            self.reassign_cluster_id(core_cell, cell.clusterId)
            self.reassign_cluster_id(ap_cell, core_cell.clusterId)
        for cell in self.cellListByType(self.BASAL_PRIME): #BASAL_PRIME
            core_cell = self.fetch_cell_by_id(cell.id+n)
            ap_cell = self.fetch_cell_by_id(cell.id+(2*n))
            self.reassign_cluster_id(core_cell, cell.clusterId)
            self.reassign_cluster_id(ap_cell, core_cell.clusterId)
            

#*******************************************************************************************************#
#SURFACE CONSTRAINT
class surface_constraint(SteppableBasePy):   
    """
    Defines a surface constraint on apical compartments. Penalizes the deviations of apical cell surface
    from the initial surface area.

    Arguments:
      _frequency             -- determines how often this Steppable is executed
      _lambda_surface_apical -- strength of surface constraint for apical compartment of mesodermal cells
    
    """
    def __init__(self,_frequency, _lambda_surface_apical):
        SteppableBasePy.__init__(self,_frequency)
        self.lambda_surface_apical = _lambda_surface_apical
        
    def start(self):
        self.apical_surface_dict = {}
        for cell in self.cell_list_by_type(self.APICAL):
            self.apical_surface_dict[cell.id] = cell.surface    #stores the initial surface area of each apical cell into a dictionary
        
    def step(self, mcs):
        self.surface_constraint_apical()
     
    #surface constraint for apical compartment- defines the target surface to be the initial surface area 
    def surface_constraint_apical(self):  
        for cell in self.cell_list_by_type(self.APICAL):
            cell.targetSurface = self.apical_surface_dict[cell.id]
            cell.lambdaSurface = self.lambda_surface_apical
                    
            
#********************************************************************************************#       
#SUB-APICAL AND SUB-BASAL LONGITUDINAL STIFFNESS CONSTRAINT
class cell_internal_links(SteppableBasePy):        
    """
    Generates spring-like intercellular force linking apical-core and basal-core compartments' centre-of-mass
    representing sub-apical and sub-basal longitudinal stiffness within each cell.

    Arguments:
      _frequency                      -- determines how often this Steppable is executed
      _tot                            -- total number of cells
      t_relax                         -- relaxation time between each myosin profile
      _lambda_links_basal_central     -- strength of central basal-core links(sub-basal stiffness)
      _lambda_links_basal_peripheral  -- strength of peripheral basal-core links(sub-basal stiffness)
      _lambda_links_apical_central    -- strength of peripheral apical-core links(peripheral sub-apical stiffness)
      _lambda_links_apical_peripheral -- strength of central apical-core links(central sub-apical stiffness)
      _n_cells                        -- number of central cells

    """    
    def __init__(self,_frequency, _tot, t_relax, _lambda_links_basal_central, _lambda_links_basal_peripheral,
    _lambda_links_apical_central, _lambda_links_apical_peripheral,_n_cells):
        SteppableBasePy.__init__(self,_frequency)
        self.t_relax = t_relax
        self.lambda_links_basal_central = _lambda_links_basal_central
        self.lambda_links_basal_peripheral = _lambda_links_basal_peripheral
        self.lambda_links_apical_central = _lambda_links_apical_central
        self.lambda_links_apical_peripheral = _lambda_links_apical_peripheral
        self.tot = _tot
        self.n_cells = _n_cells
        
    def start(self):        
        self.profile_counter = 0                 #initialize the counter for myosin profile
        self.mid = ((self.tot-1)/2.0)+1          #basal cell at the ventral midline
        CELL1 =  self.fetch_cell_by_id(1)
        CELL2 =  self.fetch_cell_by_id(48)
        CELL3 =  self.fetch_cell_by_id(95)
        self.bas_core_distance = self.distance_between_cells(CELL1, CELL2)      #initial distance between basal-core
        self.ap_core_distance = self.distance_between_cells(CELL2, CELL3)       #and apical-core domains
        
        self.create_internal_links(self.CORE, self.BASAL)                       #creates internal links 
        self.create_internal_links(self.CORE_PRIME, self.BASAL_PRIME)
                                                
        
    def step(self,mcs):
        self.profile_counter = (mcs // self.t_relax)            #counts the myosin profile and updates cell stiffness
        self.delete_internal_links(self.CORE)                   #by deleting existing internal links and creating
        self.create_internal_links(self.CORE, self.BASAL)       #new links between respective apical-core and basal-core
        self.delete_internal_links(self.CORE_PRIME)             #compartments. The stiffness parameters are upadted at each myosin profile
        self.create_internal_links(self.CORE_PRIME, self.BASAL_PRIME)
   
                    
    def delete_internal_links(self, core_type):                 #function to delete existing internal links 
        for cell in self.cell_list_by_type(core_type):
            cluster_cell_list = self.get_cluster_cells(cell.clusterId)
            for cell_cmpt in cluster_cell_list:
                if (cell_cmpt.type == self.APICAL):
                    self.focalPointPlasticityPlugin.deleteInternalFocalPointPlasticityLink(cell, cell_cmpt)
                    
   
    def create_internal_links(self, core_type, basal_type):     #function to create new internal links between respective cell compartments
        for cell in self.cell_list_by_type(core_type):
        
            #create links for central cells
            if ((1*self.tot)+self.mid-self.n_cells <= cell.id <= (1*self.tot)+self.mid+self.n_cells):
                cluster_cell_list = self.get_cluster_cells(cell.clusterId)
                for cell_cmpt in cluster_cell_list:
                    if (cell_cmpt.type == basal_type):
                        link = self.focalPointPlasticityPlugin.createInternalFocalPointPlasticityLink(cell, cell_cmpt, self.lambda_links_basal_central[int(self.profile_counter)], self.bas_core_distance, 200)
                        
                    if (cell_cmpt.type == self.APICAL):
                        link = self.focalPointPlasticityPlugin.createInternalFocalPointPlasticityLink(cell, cell_cmpt, self.lambda_links_apical_central[int(self.profile_counter)], self.ap_core_distance, 200)
        
            #create links for peripheral cells
            elif ((1*self.tot)+self.mid-12 <= cell.id < (1*self.tot)+self.mid-self.n_cells):
                cluster_cell_list = self.get_cluster_cells(cell.clusterId)
                for cell_cmpt in cluster_cell_list:
                    if (cell_cmpt.type == basal_type):
                        link = self.focalPointPlasticityPlugin.createInternalFocalPointPlasticityLink(cell, cell_cmpt, self.lambda_links_basal_peripheral[int(self.profile_counter)], self.bas_core_distance, 200)
                        
                    if (cell_cmpt.type == self.APICAL):
                        link = self.focalPointPlasticityPlugin.createInternalFocalPointPlasticityLink(cell, cell_cmpt, self.lambda_links_apical_peripheral[int(self.profile_counter)], self.ap_core_distance, 200)

            #create links for peripehral cells
            elif ((1*self.tot)+self.mid+self.n_cells < cell.id <= (1*self.tot)+self.mid+12):
                cluster_cell_list = self.get_cluster_cells(cell.clusterId)
                for cell_cmpt in cluster_cell_list:
                    if (cell_cmpt.type == basal_type):
                        link = self.focalPointPlasticityPlugin.createInternalFocalPointPlasticityLink(cell, cell_cmpt, self.lambda_links_basal_peripheral[int(self.profile_counter)], self.bas_core_distance, 200)
                        
                    if (cell_cmpt.type == self.APICAL):
                        link = self.focalPointPlasticityPlugin.createInternalFocalPointPlasticityLink(cell, cell_cmpt, self.lambda_links_apical_peripheral[int(self.profile_counter)], self.ap_core_distance, 200)


#********************************************************************************************#  
#SPRING-LIKE INTERCELLULAR FORCE CONSTRAINT FOR APICAL CONSTRICTION
class apical_constriction(SteppableBasePy):  
    """
    Generates spring-like intercellular force linking centre-of-mass of apical-apical compartments' of neighboring
    cells(for 19 central mesodermal cells). 

    Arguments:
      _frequency          -- determines how often this Steppable is executed
      lam_d               -- minimum value of constriction force
      peak_val_multiplier -- multiplier for the myosin profile
      target_distance     -- target distance between neighboring apical compartments' COMs
      t_relax             -- relaxation time after each myosin profile is added

    """    
    def __init__(self,_frequency, lam_d, peak_val_multiplier, target_distance, t_relax):
        SteppableBasePy.__init__(self,_frequency)
        self.lam_d = lam_d; self.peak_val_multiplier = peak_val_multiplier; self.target_distance = target_distance;
        self.t_relax = t_relax
        
    def start(self):
        self.mid = 24                   #basal cell at the ventral midline
        self.tot = 47                   #total number of cells
        
        #Parameters to generate time-varying myosin profiles(9 profiles)
        self.width = [8.392684812435327, 6.8575263251319045, 5.875012663000353, 5.344166722429383, 5.164011399807708, 5.233569591524036, 5.451864193967082, 5.717918103525555, 5.930754216588168, 5.930754216588168, 5.930754216588168, 5.930754216588168, 5.930754216588168]
        self.steep = [0.44488321293429717, 0.3982524926040595, 0.4047478046095632, 0.4498647719298246, 0.5190990175438597, 0.5979461644306845, 0.6719018355693154, 0.7264616539387686, 0.74712124251806, 0.74712124251806, 0.74712124251806, 0.74712124251806, 0.74712124251806]
        self.peak_value = [0.5407899896800878, 4.267370315789474, 16.792480689370475, 36.306931859649126, 61.001534575851394, 89.0670995872033, 118.69443764293084, 148.07435949226007, 175.39767588441697, 175.39767588441697, 175.39767588441697, 175.39767588441697, 175.39767588441697]
        
        self.link_list_dict = {}                #initialize the dictionary to store links
        i = 0                                   #initialize link counter
        #creates links between neighboring apical compartments for 19 mesodermal cells(initialization)
        for cell in self.cell_list_by_type(self.APICAL):        
            if ((2*self.tot)+self.mid-10 <= cell.id <= (2*self.tot)+self.mid+9):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor:
                        if (neighbor.type == self.APICAL):
                            if neighbor.id > cell.id:
                                link = self.focalPointPlasticityPlugin.createFocalPointPlasticityLink(cell, neighbor, 20, self.target_distance, 200)
                                i += 1
                                self.link_list_dict[tuple([cell.id, neighbor.id])] = i
        
    def step(self,mcs):
        if (mcs <= 13*self.t_relax):
            index = int(mcs // self.t_relax)
            #the parameters for Myosin levels
            s = self.steep[index]
            w = self.width[index]
            peak_val = self.peak_value[index]
            
            #internal-links are created between neighboring apical compartments with myosin levels updated 
            #every 2000 MCS(relaxation time) for the 19 mesodermal cells
            for cell in self.cell_list_by_type(self.APICAL):
                if ((2*self.tot)+self.mid-10 <= cell.id < (2*self.tot)+self.mid):
                    for link in self.get_focal_point_plasticity_data_list(cell): 
                        initiated_cell = link.neighborAddress
                        initiator_cell = cell
                        if initiated_cell.id > cell.id:
                            x_val = self.link_list_dict[tuple([initiator_cell.id, initiated_cell.id])]    
                            #calculate the contact area between the respective apical compartment and Medium
                            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(initiator_cell):
                                if not neighbor:
                                    area_initiator = common_surface_area/2.0    
                            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(initiated_cell):
                                if not neighbor:
                                    area_initiated = common_surface_area/2.0
                            #total contact area        
                            area = area_initiator + area_initiated 
                            #strength of the intercellular force(time varying myosin levels)
                            lam_dis = self.lam_d + (self.lam_dist(11-x_val, s, w)*self.peak_val_multiplier*peak_val)   
                            #update the strength of apical-apical link parameters using myosin profile
                            self.focalPointPlasticityPlugin.setFocalPointPlasticityParameters(cell, link.neighborAddress, lam_dis/area, self.target_distance, 200.0)
                            
                
                elif ((2*self.tot)+self.mid <= cell.id <= (2*self.tot)+self.mid+9):
                    for link in self.get_focal_point_plasticity_data_list(cell): 
                        initiated_cell = link.neighborAddress
                        initiator_cell = cell  
                        if initiated_cell.id > cell.id:   
                            x_val = self.link_list_dict[tuple([initiator_cell.id, initiated_cell.id])]
                            #calculate the contact area between the respective apical compartment and Medium    
                            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(initiator_cell):
                                if not neighbor:
                                    area_initiator = common_surface_area/2.0
                            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(initiated_cell):
                                if not neighbor:
                                    area_initiated = common_surface_area/2.0
                            #total contact area  
                            area = area_initiator + area_initiated
                            #strength of the intercellular force(time varying myosin levels)
                            lam_dis = self.lam_d + (self.lam_dist(x_val-10, s, w)*self.peak_val_multiplier*peak_val)
                            #update the strength of apical-apical link parameters using myosin profile
                            self.focalPointPlasticityPlugin.setFocalPointPlasticityParameters(cell, link.neighborAddress, lam_dis/area, self.target_distance, 200.0)
                            
                cell.connectivityOn = True       #prevents fragmentation of apical compartments 
    
    #defines the function that describes the myosin level profiles
    def lam_dist(self,x,s,w):
        lam_func = (1+np.exp(-s*w))/(1+np.exp(s*(x-w)))
        return lam_func


#*******************************************************************************************************#
#METRIC-FURROW DEPTH
class metric(SteppableBasePy):
    """
    Returns furrow depth as the normalised distance between the apical-most coordinate of the central mesodermal 
    cell and the lowest point of the mesoderm from the initial condition in units of initial apical-basal distance.

    Arguments:
      _frequency  -- determines how often this Steppable is executed
    
    """    
    def __init__(self,_frequency):
        SteppableBasePy.__init__(self,_frequency)
        
    def start(self):
        self.mid = 118                                  #cell-id of the apical compartment at the ventral midline
        mid_cell = self.fetch_cell_by_id(self.mid)
        self.yCOM_initial = mid_cell.yCOM               #initial y-center of mass of the apical compartment at the ventral midline
        
        self.time_list = []
        self.depth_data_list = []   
        
        CELL1 =  self.fetch_cell_by_id(1)
        CELL3 =  self.fetch_cell_by_id(95)
        self.bas_ap_distance = self.distance_between_cells(CELL1, CELL3)     #initial distance between apical-basal compartments
        
    def step(self,mcs):
        mid_cell = self.fetch_cell_by_id(self.mid)
        self.yCOM_now = mid_cell.yCOM                   #y-center of mass of the cell at the ventral midline at each time step
        self.time_list.append(mcs)
        self.depth_data_list.append(self.yCOM_now-self.yCOM_initial)    #shift in y-COM at each time step(furrow depth) 

    def finish(self):
        output_file, path__ = self.open_file_in_simulation_output_folder("furrow_depth_data.txt", mode='w')
        for i in range(len(self.time_list)):
            output_file.write('{} {}\n'.format(self.time_list[i],self.depth_data_list[i]/(self.bas_ap_distance))) #furrow depth as a fraction of initial apical-basal distance saved into a .txt file
                   

#*******************************************************************************************************#
#GENERATES .PIFF FILE
class piff_generator(SteppableBasePy):
    """
    Generates Potts Initial File(a lattice initialization file) that lays out the lattice pixels to cells.
    Can be used to generate simulation output snapshots.
    
    Arguments:
    _frequency -- determines how often this Steppable is executed
    _Time -- Total time
    
    """ 
    def __init__(self,_frequency,_Time):
        SteppableBasePy.__init__(self,_frequency)
        self.Time = _Time
        
    def step(self,mcs):
        self.get_all_pixels()     #get pixels without using pixel tracker plugin
        self.save_piff(mcs)
                    
    def finish(self):
        self.save_piff(self.Time)
    
    def save_piff(self,mcs):
        out_folder = self.output_dir
        FileName = out_folder+"/PiffFile_"+str(mcs)+".piff"
        piffPath=Path(out_folder).joinpath(FileName)
        with open(piffPath, 'a') as fout:
            for i in self.L:
                pixel_list = self.L[i]
                cell = self.fetch_cell_by_id(i)   
                name = self.get_type_name_by_cell(cell)
                for pixel in pixel_list:
                    x = pixel[0]; y = pixel[1]; z = pixel[2]          
                    fout.write("%d %d %s %d %d %d %d %d %d \n" % (cell.clusterId, cell.id, name, x, x, y, y, z, z))
    
    def get_all_pixels(self):
        self.L = {}   
        for cell in self.cell_list:
            self.L[cell.id] = []
        for x, y, z in self.every_pixel():
            cell = self.cell_field[x, y, 0]
            if cell:
                self.L[cell.id].append([x,y,0])

#*******************************************************************************************************#
#*******************************************************************************************************#
#CHECK FOR ARTEFACTS IN SIMULATION      
class CheckArtefactsSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
   
    def checkApical(self,core):
        cluster_cell_list = self.get_cluster_cells(core.clusterId)
        for cell in cluster_cell_list:
            if cell.type == self.APICAL:
                return True, cell  # has apical compartment
        print(">>>>>>>>>>>>>>>>>>> NO APICAL COMPT")
        return False, 0     # does NOT have apical compartment
 
    def scanPixels(self,cell):
        dy=20
        if cell.type==self.APICAL: dy=10
        xx=[self.dim.x,0]; yy=[self.dim.y,0]
        for x in range(int(cell.xCOM)-10,int(cell.xCOM)+10):
            for y in range(int(cell.yCOM)-dy,int(cell.yCOM)+dy):
                if self.cell_field[x,y,0] and self.cell_field[x,y,0].id==cell.id:
                    xx[0]=min(xx[0],x); xx[1]=max(xx[1],x)
                    yy[0]=min(yy[0],y); yy[1]=max(yy[1],y)
        return xx, yy      
 
    def checkFingers(self,core,apical):
        S = 6
        if (self.dim.x/2-18)<core.xCOM<(self.dim.x/2+18):
            cx,cy = self.scanPixels(core)
            ax,ay = self.scanPixels(apical)
            if ay[1]>(cy[0]+S):
                print(">>>>>>>>>>>>>>>>>>>>>>> FINGER DETECTED")
                return False   # fingers detected
        return True   # no fingers detected
   
    def step(self, mcs):
        print(mcs,"> > > > > >, checking for artefacts")
        artefact = 0
        for core in self.cell_list_by_type(self.CORE,self.CORE_PRIME):
            ok, apical = self.checkApical(core)
            if ok:
                if not self.checkFingers(core,apical):
                    artefact = 2; break
            else:
                artefact = 1; break
       
        if artefact==1:
            print(mcs,"STOPPING DUE TO ARTEFACT DETECTION - NO APICAL")
            output_file, path__ = self.open_file_in_simulation_output_folder("furrow_depth_data_FAILED-APICAL.txt", mode='w')
            self.stop_simulation()  
        elif artefact==2:
            print(mcs,"STOPPING DUE TO ARTEFACT DETECTION - FINGER")
            output_file, path__ = self.open_file_in_simulation_output_folder("furrow_depth_data_FAILED-FINGER.txt", mode='w')
            self.stop_simulation()
