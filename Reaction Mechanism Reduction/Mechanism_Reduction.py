# Reaction Mechanism Reduction of GRI 3.0

"""
The following Steps ensue :
1. Extract Input State Variables from a Data File & Initialize an Ideal Gas Reactor Network.
2. Simulate a 0-D Reactor for the above & Calculate the Maximum Temperature & Ignition Delay
   for each of the cases, along with Sensitivities to Temperature.
3. Sort the Reactions based on absolute values of Sensitivities for the un-reduced mechanism in each case.
4. For each case, cumulatively add single Reactions to the Mechanism & recalculate 
   Maximum Temperature & Ignition Delay.
5. If Reactions are not occuring due to lesser number of reactions in the mechanism, Ignition Delay will
   not be calculated. This Exception is caught using NameError.
6. Check whether the calculated values fall within the Tolerance Limits of the referrence values. Only when
   all calculated values meet the Tolerance Criteria, store the number of reactions for the Reduced Mechanism.

The following short-hand notations have been used :
t_ig : Time of Ignition / Ignition Delay
rid : Reaction Index
rxn : Reaction Equation
sT : Sensitivity wrt. Temperature
mech_red : Reduced Mechanism

"""
import os
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

#######################################################################

# Class to import data from a text file
class File_Reader :
    def __init__(self, file_name) :
        self.file_name = file_name
        self.input_dict = {}

    def validate_input_file(self) :
        return os.path.isfile(self.file_name)

    def Read(self) :
        if self.validate_input_file() :
            for line in open(self.file_name, 'r') :
                self.input_dict[line.split()[0]] = [float(x) for x in line.split()[1:]]  

#######################################################################

# Function 1 : To perform 0-D Reactor Simulation
# ref is a Boolean used to switch between the if & else statements
def ZeroD_Reactor(gas, T, P, phi, t, ref) :
    gas.TP = T, P
    gas.set_equivalence_ratio(phi, "CH4:1", {'O2':1, 'N2':3.76})
    r = ct.IdealGasReactor(gas)
    sim = ct.ReactorNet([r])
    sim.rtol = 1e-6
    sim.atol = 1e-15
    T_max = 0
    t_ig_flag = True

    # Sensitivity Analysis on Reactor Network :
    # The If statement returns Max Temperature, Ignition 
    # Delay & Max Sensitivity of each Reaction
    if ref :
        for i in range(gas.n_reactions) :
            r.add_sensitivity_reaction(i)
        
        sim.rtol_sensitivity = 1e-6
        sim.atol_sensitivity = 1e-6

        sT_max = [0] * gas.n_reactions

        for time in t :
            sim.advance(time)
            
            for j in range(gas.n_reactions) :
                sT = sim.sensitivity('temperature', j)

                if np.abs(sT) >= np.abs(sT_max[j]) :
                    sT_max[j] = sT

            if gas.T > T_max :
                T_max = gas.T
                
            if gas.T >= T+400 and t_ig_flag :
                t_ig = time
                t_ig_flag = False
                
        return[T_max, t_ig, sT_max]

    # The Else statement returns only Max Temperature &
    # Max Sensitivity of each Reaction
    else :
        for time in t :
            sim.advance(time)

            if gas.T > T_max :
                T_max = gas.T

            if gas.T >= T+400 and t_ig_flag :
                    t_ig = time
                    t_ig_flag = False
                
        return[T_max, t_ig]

#######################################################################

# Function 2: To perform Sorting on 2-D Arrays
def group_sorting(rxn) :
    
    Rid_sorted = []
    for j in range(len(rxn[0])) :
        for i in range(len(rxn)) :
            if rxn[i][j] not in Rid_sorted:
                Rid_sorted.append(rxn[i][j])
    return Rid_sorted

#######################################################################

# Main body of the Program
# Reading Data from the text file
a = File_Reader('Parameters.txt')
a.Read()

# Placing the data in Arrays
T_in = [x for x in a.input_dict['Temperature']]
P_in = [1e5*x for x in a.input_dict['Pressure']]
phi_in = [x for x in a.input_dict['Equivalence_Ratio']]

# This dictionary will store the (T,P,phi) values for which
# the entire program has to be repeated
case_dict = {}
case = 1

# Updating the above Dictionary. Since there are 3 values each,
# there will be 3! = 27 cases
for T in T_in :
    for P in P_in : 
        for phi in phi_in :
            case_dict.update({case : [T, P, phi]})
            case += 1

gas = ct.Solution('gri30.xml')

# 2-D Array containing Reaction Indices of all Reactions, for each case
# [[All 325 Reac Eqns]_case1, ..., [All 325 Reac Eqns]_case27]
rxn = [gas.reaction_equations() for i in range(len(case_dict))]

# 2-D Array containing Reaction Indices of all Reactions, for each case
# [[All 325 Reac Indices]_case1, ..., [All 325 Reac Indices]_case27]
rid = [np.linspace(1, gas.n_reactions, gas.n_reactions).astype(int) for i in range(len(case_dict))]

# Creating a Species Object from existing definitions of Mechanism File
spec = ct.Species.listFromFile('gri30.xml')

# Creating & Updating a Dictionary to contain all the Reactions
rxn_dict = {}
for i in range(1, gas.n_reactions+1) :
    rxn_dict.update({i : gas.reaction_equation(i-1)})

# Claculating Referrence Values for each case along with Sensitivity Analysis
# Referrence Values are calculated with the original, un-reduced Mechanism
t = np.linspace(0, 5, 5001)
T_ref = []
t_ig_ref = []

# 2-D Array meant to contain Sensitivity Values of all Reactions, for each case
# initialized with all elements as 0
sT = [[0 for i in range(gas.n_reactions)] for j in range(len(case_dict))]

# Calculating & Storing 2 Referrence Values(Temp, delay) for each case
for case in case_dict :
    
    print('Calculating Referrence Values for Case : '+str(case)+'...', sep='', end='', flush=True)
    
    # Fetching the values of T,P,phi from the case_dict
    T, P, phi = case_dict[case]
    
    # Calculating Max Temperature, Ignition Delay & Max Sensitivity using Function 1
    T_max, t_ig, sT_max = ZeroD_Reactor(gas, T, P, phi, t, True)
    
    T_ref.append(T_max)
    t_ig_ref.append(t_ig)
    
    # Curating all three 2-D Arrays with Referrence Values
    sT[case-1] = sT_max
    rid[case-1] = [r for s, r in sorted(zip(np.abs(sT_max), rid[case-1]), reverse=True)]
    rxn[case-1] = [r for s, r in sorted(zip(np.abs(sT_max), rxn[case-1]), reverse=True)]
    
    print('Done')

# Sorting Reactions by Sensitivity for each case using Function 2
Rid_sorted = group_sorting(rid)
Reactions = [rxn_dict[i] for i in Rid_sorted]

# Calculating Maximum Temperature & Ignition Delay by cumulatively adding 
# reactions to the mechanism for each case
T_case = [[0 for i in range(gas.n_reactions)] for j in range(len(case_dict))]
t_ig_case = [[0 for i in range(gas.n_reactions)] for j in range(len(case_dict))]

T_tol = 0.05
t_ig_tol = 0.05

# We shall consider a minimum of 40 Reactions
min_red = 40
red_flag = True
mech_red = [0 for i in range(len(case_dict))]
ct_add = 1
rxns = []

for case in case_dict :
    
    for r in Rid_sorted :

        rxns.append(gas.reaction(r-1))
        gas1 = ct.Solution(thermo='IdealGas',kinetics='GasKinetics',species=spec,reactions=rxns)
        T, P, phi = case_dict[case]
        print('Case : '+str(case)+', Reactions : '+str(len(rxns)))

        # Catching Exceptions in case no Reaction occurs & Ignition Delay is not calculated
        try :
            T_max, t_ig = ZeroD_Reactor(gas1, T, P, phi, t, False)
            
            T_case[case-1][len(rxns)-1] = T_max
            t_ig_case[case-1][len(rxns)-1] = t_ig
            
            T_err = np.abs(T_max - T_ref[case-1])*100/T_ref[case-1]
            t_ig_err = np.abs(t_ig - t_ig_ref[case-1])*100/t_ig_ref[case-1]

            if T_err <= T_tol and t_ig_err <= t_ig_tol and ct_add >= min_red and red_flag :
                mech_red[case-1] = ct_add
                red_flag = False
        
        except NameError :
            pass

        ct_add += 1
    
    ct_add = 1
    rxns = []
    red_flag = True

# Plotting the Results
count = np.arange(40, 325, 1)

for case in case_dict :
    
    fig, ax = plt.subplots()

    plt.subplot(2, 1, 1)
    title = 'Case : '+str(case)+', T = '+str(case_dict[case][0])+' K, P = '+str(format('%.2f'%(1e-5*case_dict[case][1])))+' Bar, phi = '+str(case_dict[case][2])
    plt.title(title, fontweight='bold')
    plt.plot(count, t_ig_case[case-1][40:], color='blue', label='Ignition Delay')
    plt.xlabel('Number of Reactions', fontweight='bold')
    plt.ylabel('Ignition Delay [s]', fontweight='bold')
    plt.hlines(y=t_ig_ref[case-1], xmin=30, xmax=350,color='red', linestyle='--', label='Referrence Value')
    plt.legend()
    plt.grid('on')
    plt.tight_layout()

    plt.subplot(2, 1, 2)
    plt.plot(count, T_case[case-1][40:], color='black', label='Maximum T')
    plt.xlabel('Number of Reactions', fontweight='bold')
    plt.ylabel('Ignition Delay [s]', fontweight='bold')
    plt.hlines(y=T_ref[case-1],xmin=30, xmax=350, color='red',linestyle='--', label='Referrence Value')
    plt.legend()
    plt.grid('on')
    plt.tight_layout()

    caption = 'Case_'+str(case)+'.png'
    plt.savefig(caption, format='png', dpi=600)

plt.show()







