#!/usr/bin/env python
# coding: utf-8

# # Workspace Setup and Library Imports

# In[68]:


aprx = arcpy.mp.ArcGISProject('Current')
wkspc = arcpy.env.workspace
Scratch = aprx.homeFolder
gdb = aprx.defaultGeodatabase


# In[2]:


import arcpy
from gurobipy import *
import pandas as pd
from arcgis.features import GeoAccessor, GeoSeriesAccessor
import os.path
import matlab.engine
eng = matlab.engine.start_matlab()


# # Data Sourcing Fx

# In[32]:


#Check for Nodes_2 if on cycle 2 or more, else use Node feature input
def Data_Sourcing():
    Nodes_FP = r'C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\Nodes_2.csv'
    if os.path.exists(Nodes_FP):
        Nodes = pd.read_csv(Nodes_FP)
        Nodes.drop(Nodes.columns[[0]], axis = 1, inplace=True)
    else:
        Nodes = pd.DataFrame.spatial.from_featureclass(r'C:\Users\Penelope\Box\NSF_D-ISN\Data\DISN_Setup\DISN_Setup\DISN_Setup\DISN_Setup.gdb\NodeInfo')
        Nodes.drop(['OBJECTID'], axis = 1, inplace=True)
        #Nodes['Timestep'] = 1
    return Nodes


# Data to be optimized:
# -Node Info as a pandas dataframe. 
# -At each timestep new Node Info is to be generated and written as a dataframe including a field for timestep. 
# -At each timestep this table will be written to a text file 
# -As the NodeInfo solution dataframes are generated a composite solution dataframe should also be compiled.

# # MT-MCLP Fx

# In[33]:


def MTMCI_func(NodesDF, timestep):
    """Writes Interdicted Nodes as text file."""
    MTMCI = Model()
    MTMCI.setParam('OutputFlag', 0)    
    x = {}
    y = {}

    #IDarray = Nodes['ID'].values
    #IDset = set(IDarray)
    ID = NodesDF['ID'].squeeze()

    FPtype = set(range(1, 9))
    j = ID
    i = ID
    x = MTMCI.addVars(j, FPtype, vtype=GRB.BINARY, name = f'Node{j,FPtype}')
    y = MTMCI.addVars(i, FPtype, vtype=GRB.BINARY, name = f'Demand{i,FPtype}')

    ##Create a[i,t] dictionary 
    #make a dataframe with just the ait# and ID field
    Ait = NodesDF.filter(regex='^ait', axis=1)
    ID_Ait= Ait.join(ID, how='right')

    #Dataframe to MultiIndex TupleDict e.g. {(i,t):a, (i,t):a...}
    #stripping 'ait' string from list, converting to int
    IDs = ID_Ait.ID.tolist()
    others = list(ID_Ait.columns)
    others.remove("ID")
    others_t = list(map(lambda st: str.replace(st, 'ait', ''), others))
    t_int = list(map(lambda ele : int(ele), others_t))
    index_tuples = [(ID, other) for ID in IDs for other in t_int]
    multi_ix = pd.MultiIndex.from_tuples(index_tuples)
    df1 = pd.DataFrame(ID_Ait[others].values.ravel(), index=multi_ix, columns=["data"])
    a = df1.to_dict()["data"]

    ##Create U[j] lists for each node, list the Type#s that have 0 as a value for that ID# 
    #e.g. {5:[1,2,4,5,6,7,8]} ; I changed to {5: [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]}
    Uj = NodesDF.filter(regex='^Type',axis=1)
    Uj_ID= Uj.join(ID, how ='right')
    U = Uj_ID.set_index('ID').T.to_dict('list')
    ##The x[j,t] values are = 0 if type can't locate at that node (13)
    for j,t in x:
        if t in U[j]==0:
            UConstraint = MTMCI.addConstr(x[j,t], GRB.EQUAL, 0, f'FPtype {t} Constraint')

    ##Producer/Consumer Constraints
    consumers = [1,156,161,162,163]
    for j in consumers:
        for t in range (1,9): 
            MTMCI.addConstr(x[j,t], GRB.EQUAL, 0, f'Consumer{j,t}')

    ###Locate Number of Force Packages
    fpTdict ={1:2,2:2,3:1,4:1,5:1,6:1,7:1,8:1}
    for t in range(1,9):
        MTMCI.addConstr(quicksum(x[j,t] for j in range(1,164)), GRB.EQUAL, fpTdict[t], f'Type Force Packages{j,t}')


    ##Covering Constraints by Force Package Type  (14)
    for i,t in y:
        for j in range(1,164):
            MTMCI.addConstr(x[i,t],GRB.GREATER_EQUAL, y[i,t], f'Covering Constraint{i,t}')

    ##Set the Objective Function
    MTMCI.setObjective(quicksum(a[i,t]*y[i,t] for i,t in y), GRB.MAXIMIZE) 

    ##Update
    MTMCI.update()
    MTMCI.optimize()

    intSites = []
    for j, t in x:
        if x[j,t].X == 1:
            intSites.append(str(j))
            print(f'Node {j} | Type {t}') 

    obj = MTMCI.getObjective()
    print ("Objective Value: " + str(obj.getValue()))

    #write interdicted nodes to text file
    SolutionFP = r'C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\MTMCI_IntNodes\MTMCI_IntNodes.txt'
    with open (SolutionFP, 'w') as q:
        for site in intSites:
            q.write('%s\n' % site)
            
    print("MTMCI for timestep {} is complete.".format(timestep))


# # NarcoLogic Fx

# In[34]:


def NarcoLogic_Fxs(number):
    NLtimestep = float(number + 1)
#     import matlab.engine
#     eng = matlab.engine.start_matlab()

    if NLtimestep == 2.0:
        #Call NarcoLogic Initialization Fx to generate workspace variables
        """writes workspace.m file called in NarcoLogic_Dynamic Fx"""
        eng.addpath(r'C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic')
        eng.NarcoLogic_initialize_python(nargout=0)  

        eng.addpath(r'C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic')
        tflow = eng.NarcoLogic_dynamics_python(NLtimestep, nargout=1)

        #Pull up latest Matlab workspace to convert Tflow table to python variable
        eng.eval("NL_Tflow = table2struct(Tflow, 'ToScalar',true);", nargout=0)

        #Convert Tflow data to pandas database
        NarcoLogic_Tflow = eng.workspace["NL_Tflow"]
        Tflow = pd.DataFrame.from_dict(NarcoLogic_Tflow)
        print('NarcoLogic simulation for timestep {} complete!'.format(NLtimestep))    
        return Tflow
    
    else:
        eng.addpath(r'C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic')
        eng.NarcoLogic_dynamics_python(NLtimestep, nargout=1)

        #Pull up latest Matlab workspace to convert Tflow table to python variable
        eng.eval("NL_Tflow = table2struct(Tflow, 'ToScalar',true);", nargout=0)

        #Convert Tflow data to pandas database
        NarcoLogic_Tflow = eng.workspace["NL_Tflow"]
        Tflow = pd.DataFrame.from_dict(NarcoLogic_Tflow)

        print('NarcoLogic simulation for timestep {} complete!'.format(NLtimestep))    

        return Tflow


# # Data Processing
# 1) Merge Nodes and TFlow, retain all rows<br>
# 2) Create dicts/columns to allocate updated flow based on positive type per ID<br>
# 3) Perform the 25% rule (.75 * expected flow value (PMNode) + .25 * actual flow (InitFlow) = the new PMNode value<br>
# 4) Drop unnecessary columns (Start_Node, IntitFlow, DTO_y, SumType, PMNode_2, Update Flow, ObjectID)<br>
# 5) Consolidate rows by end node
# 6) Update Timestep by 1<br>
# 7) Concatinate updated flow data with the all previous data in a comprehensive dataframe, write to csv<br>
# 8) Write Nodes_2 as csv to be called in as Nodes in next cycle

# In[35]:


def Data_Processing(Tflow, Nodes, i):
    #cast the returned TFlow table and timestep to py variables
    tf2 = Tflow
    NLtimestep = i + 1

    #convert object type column values in TFlow to float by first converting to string and omitting special characters, drop DTO
    tf2['End_Node'] = tf2['End_Node'].astype(str)
    tf2['End_Node'] = tf2['End_Node'].str[1:4]
    tf2['IntitFlow'] = tf2['IntitFlow'].astype(str)
    tf2['IntitFlow'] = tf2['IntitFlow'].str[1:4]
    tf2['End_Node'] = tf2['End_Node'].astype(float)
    tf2['IntitFlow'] = tf2['IntitFlow'].astype(float)
    tf2.drop(['DTO'], axis=1, inplace=True)

    #combine Nodes with tf2 (tflow table) and write to new df
    Nodes_2 = Nodes.merge(tf2, on=['End_Node'], how = 'left')

    #Get the Type value by node stored as a dict by ID e.g. {1: [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], 2: [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]}
    ID = Nodes_2['ID']
    Types= Nodes_2.filter(regex ='^Type', axis = 1)
    ID_Type= Types.join(ID, how='right')
    ID_TypeDict = ID_Type.set_index('ID').T.to_dict('list')

    # Sum the positive type values by ID e.g. {2: 1.0, 3: 1.0, 4: 2.0,...} to divide the updated flow by per node.
    SumType = {k:sum(v) for k,v in ID_TypeDict.items()}

    #convert TypeSum dict to column
    Nodes_2['SumType'] = Nodes_2['ID'].map(SumType)

    #Perform the 25% rule (.75 * expected flow value (PMNode) + .25 * actual flow (InitFlow) =  PMNode2
    Nodes_2['PMNode2'] = (Nodes_2['PMNode'] *.75) + (Nodes_2['IntitFlow'] *.25)

    #Divide the updated flow value by number of nodes w/ Type = 1
    Nodes_2['UpdateFlow'] = Nodes_2['PMNode2']/Nodes_2['SumType']

    #allocate flow to positive ait/type values
    TypeNums = [1, 2, 3, 4, 5, 6, 7, 8]
    for num in TypeNums:
        Nodes_2['ait{}'.format(num)] = Nodes_2['Type{}'.format(num)] * Nodes_2['UpdateFlow']

    #Update Timestep
    Nodes_2['Timestep'] = NLtimestep

    #Remove Unnecessary Columns, drop duplicate rows, fill NaNs w 0
    Nodes_2['PMNode'] = Nodes_2['PMNode2']
    Nodes_2.drop(['Start_Node', 'IntitFlow', 'SumType', 'PMNode2', 'UpdateFlow'], axis=1, inplace=True)
    Nodes_2.drop_duplicates(subset=['ID'], inplace=True)
    Nodes_2.fillna(0, inplace=True)
    
    #Concatinate Nodes_2 table with initial Node table if TS 1->2, otherwise Concat CompDF with Nodes_2 and overwrite CompDF csv
    CompFP = r'C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\ComprehensiveNodesDF.csv'
    if os.path.exists(CompFP):
        Comp = pd.read_csv(CompFP)
        Comp.drop(Comp.columns[0], axis = 1, inplace=True)
        ComprehensiveNodesDF = pd.concat([Comp, Nodes_2])
        ComprehensiveNodesDF.to_csv(CompFP)
    else:
        ComprehensiveNodesDF = pd.concat([Nodes, Nodes_2])
        ComprehensiveNodesDF.to_csv(CompFP)
        
    #Write Nodes_2 to csv
    Nodes_FP = r'C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\Nodes_2.csv'
    Nodes_2.to_csv(Nodes_FP)
    
    print('Data Processing for timestep {} complete!'.format(NLtimestep))


# # Call Functions

# In[39]:


times = range(1, 180)
for i in times:

    MTMCI_func(Data_Sourcing(), i)

    Data_Processing(NarcoLogic_Fxs(i), Data_Sourcing(), i) 

