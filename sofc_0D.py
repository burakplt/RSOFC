''' Installation guide of Cantera packages:
https://cantera.org/install/conda-install.html
'''
__author__ = "Burak Polat"
__email__ = "burakplt@outlook.com"

import os
import cantera as ct
from math import log, asinh, exp, pi
from copy import copy, deepcopy
from matplotlib import pyplot as plt
from gas import *
from datetime import datetime
from pathlib import Path


R = 8.314       # J/mol.K
F = 96485.3329  # Faraday's constant


class SOFC_0D():

    """ Reversible SOFC Model, 0-Dimensional.

        Args:
            fuel_gas(GasStream): The gas stream is fed to "anode"
            oxygen_gas(GasStream): The gas stream is fed to "cathode"
            operation(str): Operation mode 'SOFC' or 'SOEC'
            fixed(str): Specify the fixed variables. 'FU,V' , 'FU,A' or 'V,A'
            
    """

    def __init__(self, fuel_gas, oxygen_gas, operation = "SOFC", fixed = "FU,A"):

        self.name = "RSOFC Reactor 0-D"

        ### Gas Feeds ###
        self.fuel_gas = fuel_gas
        self.oxygen_gas = oxygen_gas
        self.orig_oxy = {"T":copy(self.oxygen_gas.T), "P":copy(self.oxygen_gas.P), "X":copy(self.oxygen_gas.gas_stream.X),"feedrate":copy(self.oxygen_gas.feed_rate) }
        self.orig_fuel = {"T":copy(self.fuel_gas.T), "P":copy(self.fuel_gas.P), "X":copy(self.fuel_gas.gas_stream.X),"feedrate":copy(self.fuel_gas.feed_rate) }
        
        ### Operational Parameters ###
        self.operation = operation      # SOFC or SOEC
        self.fixed = fixed              # Varying "FU,V", "FU,A", or "V,A"
        self.fixed_oxygen_flow = True   # If False, Oxygen/Air feed rate will be adjusted according to Q_rxn, so variable air flow
        self.cell_area = 16             # cm2
        self.FU = 0.40                   # Fuel utilisation
        self.V_cell = 0.70              # Volt
        self.dT_max = 50               # Max temperature difference along fuel channel
        self.n_fuel = 2                 # number of electrons for fuel
        self.n_oxy = 4                  # number of electrons for oxygen
        self.P = self.fuel_gas.P        # Pressure [Pa]
        self.T_op = self.fuel_gas.T + self.dT_max/2    # Mean operation temperature of the cell [K]
        self.fuel_flowrate = self.fuel_gas.feed_rate   # mol/sec
        
        ### Activation parameters ###
        self.E_act_fuel = 1e+5       # J/mol
        self.E_act_oxy = 1.2e+5      # J/mol
        self.E_act_el = 85634        # J/mol
        self.gamma_fuel = 1.344e+6   # A/cm2
        self.gamma_oxy = 2.051e+5    # A/cm2

        ### Diffusion parameters ###
        self.porosity = 0.3
        self.tortuosity = 5
        self.pore_diameter = 1e-6       # m
        self.thickness_fuel = 3.125e-5    # m
        self.thickness_oxy = 1.75e-5    # m
        self.Vd = {"H2": 6.12, "N2": 18.5, "O2": 16.3, "H2O": 13.1} # Diffusion volumes
        
        ### Ohmic loss parameters ###
        self.sigma0_el = 333.3       # 1/ohm.cm2 
        self.r_ohmic_const = 0.057   # ohm.cm2
        self.thickness_el = 1.25e-5  # m
        

    def state(self):
        
        # TODO automatically change to SOEC mode if V_gibbs < V_op 
        if self.operation == "SOFC":
            self.SOFC_mode()
        else:
            self.SOEC_mode()

    
    def SOFC_mode(self):
        time = str(datetime.now().strftime("%d-%m-%Y %H_%M_%S")) + ".txt"
        file_name = "log.txt" # os.path.join( Path(__file__).parent.absolute() , time)
        #report = open(file_name, "a")
        #report.write("*"*10+ " Simulation Started "+ "*"*10+ "\n")
        print("Fixed: ", self.fixed, "\tFixed Airflow: ", self.fixed_oxygen_flow)

        # Molar flowrate of fuel and its components
        N_fuel = self.fuel_flowrate
        x_H2_in = copy(self.fuel_gas.mole_frac('H2'))
        x_H2O_in = copy(self.fuel_gas.mole_frac('H2O'))
        x_CO_in = copy(self.fuel_gas.mole_frac('CO'))
        x_CH4_in = copy(self.fuel_gas.mole_frac('CH4'))
        x_O2_in = copy(self.oxygen_gas.mole_frac('O2'))

        N_CH4 = copy(x_CH4_in*N_fuel)
        N_CO = copy(x_CO_in*N_fuel)
        N_H2 = copy(x_H2_in*N_fuel)
        N_H2O = copy(x_H2O_in*N_fuel)
    
        J_max = (2*N_H2 + 8*N_CH4 + 2*N_CO)*F   # Max current
        
        if self.fixed == "FU,A":
            ## Fixed Fuel Utilisation and fixed total cell area ##
            ## First, the maximum Gibbs voltage is calculated to determine SOFC or SOEC mode ##

            FU = self.FU
            fuel_copy = GasStream(fed_to="fuel_channel", T=self.orig_fuel["T"],P= self.orig_fuel["P"], X=self.orig_fuel["X"],feed_rate= self.orig_fuel["feedrate"] )
            air_copy = GasStream(fed_to="oxygen_channel", T=self.orig_oxy["T"],P= self.orig_oxy["P"],X= self.orig_oxy["X"],feed_rate= self.orig_oxy["feedrate"] )
            O2_anode = GasStream(fed_to="oxygen_channel", T=self.T_op, P=self.P, X="O2:1.0") # Oxygen flow at electrolyte side, i.e at TPB

            for i in range(0, 3):

                ### Temperature is equal to Fuel feed temperature ###
                fuel_copy.feed_rate = self.orig_fuel["feedrate"]
                fuel_copy.gas_stream.TPX = self.orig_fuel["T"], self.orig_fuel["P"], self.orig_fuel["X"]
                air_copy.gas_stream.TPX = self.orig_oxy["T"], self.orig_oxy["P"], self.orig_oxy["X"]
                O2_anode.gas_stream.TPX = self.T_op, self.P, "O2:1.0"

                j = FU*J_max/self.cell_area             # Current density
                J_total = J_max*FU                      # Total current 

                N_O2_an = J_total/4/F                   # Oxygen flow from cathode to anode
                O2_anode = GasStream(fed_to="oxygen_channel", T=self.T_op, P=self.P, X="O2:1.0", feed_rate= N_O2_an) # Oxygen flow at electrolyte side, i.e at TPB
                
                Cp_f1 = copy(fuel_copy.Cp_mole())   # Heat capacity fuel inlet J/mol.K
                T_fuel_in = copy(fuel_copy.gas_stream.T)
                T_oxy_in = copy(air_copy.gas.T)

                ### Temperature of fuel side is now at T_op = T_fuel_in + delta T_max/2 ###
                fuel_copy.gas_stream.TP = self.T_op, self.P

                H_fuel_in = copy(fuel_copy.enthalpy_stream())           # Enthalpy of Fuel cell inlet, before reforming
                Cp_f2 = copy(fuel_copy.Cp_mole())                       # Heat capacity fuel inlet at T_op J/mol.K
                
                fuel_copy.gas_stream.equilibrate('TP', solver="gibbs")  # Gas reforming
                G_fuel_in = copy(fuel_copy.gibbs_stream())
                
                fuel_bulk = GasMixture( [fuel_copy, O2_anode] )         # Mixing the fuel inlet stream and oxygen at TPB
                fuel_bulk.gas_stream.equilibrate('TP', solver= "gibbs")     # Gibbs minimization at mean cell temperature, electrochemical reactions
                Cp_f3 = copy(fuel_bulk.Cp_mole())                           # Heat capacity in the fuel cell, after electrochemical reactions j/mol.K
                H_fuel_out = copy(fuel_bulk.enthalpy_stream())              # Enthalpy of Fuel cell outlet, at T_op
                G_fuel_out = copy(fuel_bulk.gibbs_stream())                 # Gibbs energy of Fuel cell outlet, at T_op

                W_DC = self.V_cell* J_total
                Q_rxn = lambda H_O2_in, H_O2_out, W_DC: H_fuel_in + H_O2_in - H_fuel_out - H_O2_out - W_DC    # Reaction heat calculation [W]
                
                if self.fixed_oxygen_flow == False:
                    # Variable air flow
                    # First check the minimum oxygen-channel flowrate based on oxygen flow to fuel side
                    ### Temperature of oxygen side at T_in air stream ###
                    air_lambda_old = 1                                                    # Initial guess
                    N_air_stoic = N_O2_an / air_copy.mole_frac("O2")                      # Minimum required amount of O2, stoichiometric minimum mol
                    N_air_in = N_air_stoic * air_lambda_old                               # Initial guess for air inlet based on minimum O2 need mol
                    small_lambda = False

                    for i in range(0, 50):
                        
                        air_copy.gas_stream.TPX = self.orig_oxy["T"], self.orig_oxy["P"], self.orig_oxy["X"]
                        air_copy.feed_rate = N_air_in                                # Set oxygen/air mole flow as N_oxy_in mol
                        Cp_oxy1 = copy(air_copy.Cp_mole())                           # Heat capacity of air stream inlet J/mol.K
                        air_copy.gas_stream.TP = self.T_op, self.P                   # Increasing T of air stream to operation temperature
                        
                        ### Temperature of oxygen side at T_op ###
                        H_oxy_in = copy(air_copy.enthalpy_stream())
                        G_oxy_in = copy(air_copy.gibbs_stream())
                        Cp_oxy2 = copy(air_copy.Cp_mole())                  # Heat capacity of air stream internal J/mol.K
                        stream_operation( air_copy, -N_O2_an, {"O2":1.0} )             # Remove the oxygen stream which is transferred to anode
                        
                        Cp_oxy3 = copy(air_copy.Cp_mole())                  # Heat capacity of air stream outlet J/mol.K
                        H_oxy_out = copy(air_copy.enthalpy_stream())
                        G_oxy_out = copy(air_copy.gibbs_stream())
                        
                        ### Temperature is now T = T_op + delta T_max/2 ###
                        fuel_bulk.gas.TP = self.T_op + self.dT_max/2, self.P
                        air_copy.gas.TP = self.T_op + self.dT_max/2, self.P
                        Cp_f4 = copy(fuel_bulk.Cp_mole())
                        Cp_oxy4 = copy(air_copy.Cp_mole())

                        dH_stack = H_fuel_out + H_oxy_out - H_fuel_in - H_oxy_in
                        Q_anode = (Cp_f1 + Cp_f2 + Cp_f3 + Cp_f4)/4* self.dT_max * N_fuel
                        Q_Reaction = Q_rxn(H_oxy_in, H_oxy_out, W_DC)
                        Q_cathode = abs(Q_Reaction) - Q_anode
                        Cp_cat = (Cp_oxy1+ Cp_oxy2+ Cp_oxy3+ Cp_oxy4)/4

                        if not small_lambda:
                            N_air_calculated = Q_cathode / Cp_cat /(abs(T_fuel_in- T_oxy_in)+ self.dT_max)
                            air_lambda_new = N_air_calculated / N_air_stoic
                            N_air_in *= air_lambda_new/air_lambda_old
                            error = copy(abs(air_lambda_old - air_lambda_new)/air_lambda_new)
                            air_lambda_old = air_lambda_new
                        
                        else:
                            break

                        if error <= 1e-3:
                            if air_lambda_old < 1:
                                air_lambda_new = 1
                                N_air_in = N_air_stoic
                                small_lambda = True
                            else:
                                break    
                        
                    
                else:
                    # Temperatur difference between inlet-oulet must be calculated
                    air_copy.gas_stream.TP = self.T_op, self.P
                    air_copy.feed_rate = self.orig_oxy["feedrate"] 
                    G_oxy_in = copy(air_copy.gibbs_stream())
                    stream_operation(air_copy, -N_O2_an, {"O2":1.0} )             # Remove the oxygen stream which is transferred to anode
                    G_oxy_out = copy(air_copy.gibbs_stream())
                    
                delta_G = G_fuel_out - G_fuel_in + G_oxy_out - G_oxy_in   # delta Gibbs fuel cell
                V_gibbs = abs(delta_G)/J_total
                if V_gibbs < self.V_cell:
                    print("Switch to SOEC mode, because V_gibbs maximum is lower than cell voltage")
                    # TODO SOEC Mode
                    break
                
                else:
                    # Outlet molar fractions
                    x_H2_out = copy(fuel_bulk.mole_frac('H2'))
                    x_H2O_out = copy(fuel_bulk.mole_frac('H2O'))
                    x_CO_out = copy(fuel_bulk.mole_frac('CO'))
                    x_CH4_out = copy(fuel_bulk.mole_frac('CH4'))
                    x_O2_out = copy(air_copy.mole_frac('O2'))
                    x_H2_eq = (x_H2_in + x_H2_out)/2 + (x_CO_in + x_CO_out)/2 + 4*(x_CH4_in + x_CH4_out)/2
                    x_H2O = (x_H2O_in + x_H2O_out)/2
                    x_O2 = (x_O2_in + x_O2_out)/2
                    
                    ASR = ( self.activation_v(j) + self.diffusion_v(j, x_H2_eq, x_O2, x_H2O) + self.ohmic_v(j) )/j
                    V_op = V_gibbs - j*ASR 
                    self.V_cell = V_op
            
            print("\nASR: ", ASR)
            print("V_gibbs: ", V_gibbs, "\nV_op: ", V_op, "\nj[A/cm2]: ",j, "\nFU: ", FU, "\nJ: ",J_total, "\nair_lambda:", "\tN_airflow:", air_copy.feed_rate)
            return V_op, j, ASR
                    
                
        elif self.fixed == "V,A":
            ## Fixed Voltage and fixed total cell area ##
            ## First, the maximum Gibbs voltage is calculated to determine SOFC or SOEC mode  ##

            FU_max = 1.0001
            FU_min = 0.0001
            err_tolerance = 1e-6
            f_ASR = lambda j, x_H2_eq, x_O2, x_H2O: ( self.activation_v(j) + self.diffusion_v(j, x_H2_eq, x_O2, x_H2O) + self.ohmic_v(j) )/j

            fuel_copy = GasStream(fed_to="fuel_channel", T=self.orig_fuel["T"],P= self.orig_fuel["P"],X= self.orig_fuel["X"],feed_rate= self.orig_fuel["feedrate"] )
            air_copy = GasStream(fed_to="oxygen_channel", T=self.orig_oxy["T"], P=self.orig_oxy["P"], X=self.orig_oxy["X"],feed_rate= self.orig_oxy["feedrate"] )
            O2_anode = GasStream(fed_to="oxygen_channel", T=self.T_op, P=self.P, X="O2:1.0") # Oxygen flow at electrolyte side, i.e at TPB

        
            def V_cell_FU(FU, converged= False):
                
                fuel_copy.feed_rate = self.orig_fuel["feedrate"]
                fuel_copy.gas_stream.TPX = self.orig_fuel["T"], self.orig_fuel["P"], self.orig_fuel["X"]
                air_copy.gas_stream.TPX = self.orig_oxy["T"], self.orig_oxy["P"], self.orig_oxy["X"]
                O2_anode.gas_stream.TPX = self.T_op, self.P, "O2:1.0"
                
                ### Temperature is equal to Fuel feed temperature ###
                
                j = FU*J_max/self.cell_area             # Current density
                J_total = J_max*FU                      # Total current 

                N_O2_an = J_total/4/F                   # Oxygen flow from cathode to anode
                O2_anode.feed_rate = N_O2_an
                
                Cp_f1 = copy(fuel_copy.Cp_mole())       # Heat capacity fuel inlet J/mol.K
                T_fuel_in = copy(fuel_copy.gas_stream.T)
                T_oxy_in = copy(air_copy.gas.T)

                ### Temperature of fuel side is now at T_op = T_fuel_in + delta T_max/2 ###
                fuel_copy.gas_stream.TP = self.T_op, self.P

                H_fuel_in = copy(fuel_copy.enthalpy_stream())           # Enthalpy of Fuel cell inlet, before reforming
                Cp_f2 = copy(fuel_copy.Cp_mole())                       # Heat capacity fuel inlet at T_op J/mol.K
                
                fuel_copy.gas_stream.equilibrate('TP')  # Gas reforming
                G_fuel_in = copy(fuel_copy.gibbs_stream())
                
                fuel_bulk = GasMixture( [fuel_copy, O2_anode] )             # Mixing the fuel inlet stream and oxygen at TPB
                fuel_bulk.gas_stream.equilibrate('TP')     # Gibbs minimization at mean cell temperature, electrochemical reactions
                Cp_f3 = copy(fuel_bulk.Cp_mole())                           # Heat capacity in the fuel cell, after electrochemical reactions j/mol.K
                H_fuel_out = copy(fuel_bulk.enthalpy_stream())              # Enthalpy of Fuel cell outlet, at T_op
                G_fuel_out = copy(fuel_bulk.gibbs_stream())                 # Gibbs energy of Fuel cell outlet, at T_op
                
                W_DC = self.V_cell* J_total
                Q_rxn = lambda H_O2_in, H_O2_out, W_DC: H_fuel_in + H_O2_in - H_fuel_out - H_O2_out - W_DC    # Reaction heat calculation [W]
                
                
                if self.fixed_oxygen_flow == False:
                    # Variable air flow
                    # First check the minimum oxygen-channel flowrate based on oxygen flow to fuel side
                    ### Temperature of oxygen side at T_in air stream ###
                    air_lambda_old = 1                                                    # Initial guess
                    N_air_stoic = N_O2_an / air_copy.mole_frac("O2")             # Minimum required amount of O2, stoichiometric minimum mol
                    N_air_in = N_air_stoic * air_lambda_old                               # Initial guess for air inlet based on minimum O2 need mol
                    small_lambda = False
                    
                    for i2 in range(0, 50):

                        fuel_bulk.TP = self.T_op, self.P
                        air_copy.gas_stream.TPX = self.orig_oxy["T"], self.orig_oxy["P"], self.orig_oxy["X"]
                        air_copy.feed_rate = N_air_in                                # Set oxygen/air mole flow as N_oxy_in mol
                        Cp_oxy1 = copy(air_copy.Cp_mole())                           # Heat capacity of air stream inlet J/mol.K
                        air_copy.gas_stream.TP = self.T_op, self.P                   # Increasing T of air stream to operation temperature
                        
                        ### Temperature of oxygen side at T_op ###
                        H_oxy_in = copy(air_copy.enthalpy_stream())
                        G_oxy_in = copy(air_copy.gibbs_stream())
                        Cp_oxy2 = copy(air_copy.Cp_mole())                  # Heat capacity of air stream internal J/mol.K
                        stream_operation( air_copy, -N_O2_an, {"O2":1.0} )             # Remove the oxygen stream which is transferred to anode
                        
                        Cp_oxy3 = copy(air_copy.Cp_mole())                  # Heat capacity of air stream outlet J/mol.K
                        H_oxy_out = copy(air_copy.enthalpy_stream())
                        G_oxy_out = copy(air_copy.gibbs_stream())
                        
                        ### Temperature is now T = T_op + delta T_max/2 ###
                        fuel_bulk.gas_stream.TP = self.T_op + self.dT_max/2, self.P
                        air_copy.gas_stream.TP = self.T_op + self.dT_max/2, self.P
                        Cp_f4 = copy(fuel_bulk.Cp_mole())
                        Cp_oxy4 = copy(air_copy.Cp_mole())

                        dH_stack = H_fuel_out + H_oxy_out - H_fuel_in - H_oxy_in
                        Q_anode = (Cp_f1 + Cp_f2 + Cp_f3 + Cp_f4)/4* self.dT_max * N_fuel
                        Q_Reaction = Q_rxn(H_oxy_in, H_oxy_out, W_DC)
                        Q_cathode = abs(Q_Reaction) - Q_anode
                        Cp_cat = (Cp_oxy1+ Cp_oxy2+ Cp_oxy3+ Cp_oxy4)/4

                        if not small_lambda:
                            N_air_calculated = Q_cathode / Cp_cat /(abs(T_fuel_in- T_oxy_in)+ self.dT_max)
                            air_lambda_new = N_air_calculated / N_air_stoic
                            N_air_in *= air_lambda_new/air_lambda_old
                            error = copy(abs(air_lambda_old - air_lambda_new)/air_lambda_new)
                            air_lambda_old = air_lambda_new
                        
                        else:
                            break

                        if error <= 1e-3:
                            if air_lambda_old < 1:
                                air_lambda_new = 1
                                N_air_in = N_air_stoic
                                small_lambda = True
                            else:
                                break
                    
                else:
                    air_copy.gas_stream.TP = self.T_op, self.P
                    air_copy.feed_rate = self.orig_oxy["feedrate"] 
                    G_oxy_in = copy(air_copy.gibbs_stream())
                    N_air_stoic = N_O2_an / air_copy.mole_frac("O2")
                    air_lambda_new = air_copy.feed_rate / N_air_stoic
                    stream_operation(air_copy, -N_O2_an, {"O2":1.0} )             # Remove the oxygen stream which is transferred to anode
                    G_oxy_out = copy(air_copy.gibbs_stream())
                    
                delta_G = G_fuel_out - G_fuel_in + G_oxy_out - G_oxy_in   # delta Gibbs fuel cell
                V_gibbs = abs(delta_G)/J_total
                x_H2_out = copy(fuel_bulk.mole_frac('H2'))
                x_H2O_out = copy(fuel_bulk.mole_frac('H2O'))
                x_CO_out = copy(fuel_bulk.mole_frac('CO'))
                x_CH4_out = copy(fuel_bulk.mole_frac('CH4'))
                x_O2_out = copy(air_copy.mole_frac('O2'))
                x_H2_eq = (x_H2_in + x_H2_out)/2 + (x_CO_in + x_CO_out)/2 + 4*(x_CH4_in + x_CH4_out)/2
                x_H2O = (x_H2O_in + x_H2O_out)/2
                x_O2 = (x_O2_in + x_O2_out)/2
                
                ASR = f_ASR(j, x_H2_eq, x_O2, x_H2O)
                V_op = V_gibbs - j*ASR
                
                if converged:
                    print("\nASR: ", ASR)
                    print("V_gibbs: ", V_gibbs, "\nV_op: ", V_op, "\nj[A/cm2]: ",j, "\nFU: ", FU, "\nJ: ",J_total, "\nair_lambda:", air_lambda_new, "\tN_airflow:", air_copy.feed_rate)
                else:
                    return V_op - self.V_cell
            
            a = FU_min
            b = FU_max
            an = a
            bn = b
            
            for i in range(0, 100):
                
                mn = (an + bn)*0.5
                f_mn = V_cell_FU(mn)
                if V_cell_FU(an)*f_mn < 0:
                    bn = mn

                elif V_cell_FU(bn)*f_mn < 0:
                    an = mn

                if abs(f_mn) < err_tolerance:
                    FU = (an+bn)*0.5
                    V_cell_FU(FU, converged= True)
                    break


    def SOEC_mode(self):
        time = str(datetime.now().strftime("%d-%m-%Y %H_%M_%S")) + ".txt"
        file_name = "log.txt" # os.path.join( Path(__file__).parent.absolute() , time)
        #report = open(file_name, "a")
        #report.write("*"*10+ " Simulation Started "+ "*"*10+ "\n")

        self.dT_max = 1
        CR = ConversionReactor()    # Stoichiometric conversion Reactor
        
        if self.fixed == "FU,A":
            ## Fixed Fuel Utilisation and fixed total cell area ##
            ## First, the maximum Gibbs voltage is calculated to determine SOFC or SOEC mode ##

            FU = self.FU
            fuel_copy = GasStream(fed_to="fuel_channel", T=self.orig_fuel["T"],P= self.orig_fuel["P"],X= self.orig_fuel["X"],feed_rate= self.orig_fuel["feedrate"] )
            air_copy = GasStream(fed_to="oxygen_channel",T= self.orig_oxy["T"],P= self.orig_oxy["P"],X= self.orig_oxy["X"],feed_rate= self.orig_oxy["feedrate"] )
            N_H2O_rcl = 2.3477E-05   # Recycle stream of H2O after water-gas shift and reforming reactions
            
            for i in range(0, 1):
    
                ### Temperature is equal to Fuel feed temperature ###
                fuel_copy.gas_stream.TPX = self.orig_fuel["T"], self.orig_fuel["P"], self.orig_fuel["X"]
                fuel_copy.feed_rate = self.orig_fuel["feedrate"]
                air_copy.gas_stream.TPX = self.orig_oxy["T"], self.orig_oxy["P"], self.orig_oxy["X"]
                
                Cp_f1 = copy(fuel_copy.Cp_mole())           # Heat capacity fuel inlet J/mol.K
                T_fuel_in = copy(fuel_copy.gas_stream.T)
                T_oxy_in = copy(air_copy.gas.T)

                ### Temperature of fuel side is now at T_op = T_fuel_in + delta T_max/2 ###
                self.T_op = self.dT_max/2 + self.fuel_gas.T
                fuel_copy.gas_stream.TP = self.T_op, self.P

                H_fuel_in = copy(fuel_copy.enthalpy_stream())           # Enthalpy of Fuel cell inlet, before reforming
                Cp_f2 = copy(fuel_copy.Cp_mole())                       # Heat capacity fuel inlet at T_op J/mol.K
                
                fuel_copy.gas_stream.equilibrate('TP')                  # Gas reforming
                N_H2O = copy(fuel_copy.mole_flow("H2O"))
                # Molar flowrate of fuel and molarfractions inlet
                N_fuel = self.fuel_flowrate
                x_H2_in = copy(self.fuel_gas.mole_frac('H2'))
                x_H2O_in = copy(self.fuel_gas.mole_frac('H2O'))
                x_CO_in = copy(self.fuel_gas.mole_frac('CO'))
                x_CH4_in = copy(self.fuel_gas.mole_frac('CH4'))
                x_O2_in = copy(self.oxygen_gas.mole_frac('O2'))
                
                J_max = -(2*(N_H2O + N_H2O_rcl))*F                      # Max current
                j = FU*J_max/self.cell_area                             # Current density
                J_total = J_max*FU                                      # Total current
                Q1 = copy(-self.V_cell*J_total)
                
                stream_operation(fuel_copy, N_H2O_rcl, {"H2O":1})
                G_fuel_in = copy(fuel_copy.gibbs_stream())
                
                CR.react( [fuel_copy], {"H2O":-1, "H2":1, "O2":0.5}, "H2O", FU )    # Electrolysis reaction
                N_H2O_1 = copy(fuel_copy.mole_flow("H2O"))
                O2_produced = copy( fuel_copy.mole_flow("O2") )
                stream_operation(fuel_copy, - O2_produced, {"O2":1.0})              # O2 produced transfer to cathode
                
                Cp_f3 = copy(fuel_copy.Cp_mole())                           # Heat capacity in the fuel cell, after electrochemical reactions j/mol.K
                H_fuel_out = copy(fuel_copy.enthalpy_stream())              # Enthalpy of Fuel cell outlet, at T_op
                G_fuel_out = copy(fuel_copy.gibbs_stream())                 # Gibbs energy of Fuel cell outlet, at T_op
                
                fuel_copy.gas_stream.equilibrate("TP")                      # Water-gas shift reaction and other possible reforming reactions
                N_H2O_2 = copy(fuel_copy.mole_flow("H2O"))
                N_H2O_rcl = N_H2O_2 - N_H2O_1                               # Water recycle
                
                W_DC = self.V_cell* J_total
                Q_rxn = lambda H_O2_in, H_O2_out, W_DC: H_fuel_in + H_O2_in - H_fuel_out - H_O2_out - W_DC    # Reaction heat calculation [W]
                
                if self.fixed_oxygen_flow == False:
                    # Variable air flow
                    # First check the minimum oxygen-channel flowrate based on oxygen flow to fuel side
                    ### Temperature of oxygen side at T_in air stream ###
                    air_lambda_old = 1                                                    # Initial guess
                    N_air_stoic = N_O2_an / air_copy.mole_frac("O2")             # Minimum required amount of O2, stoichiometric minimum mol
                    N_air_in = N_air_stoic * air_lambda_old                               # Initial guess for air inlet based on minimum O2 need mol
                    small_lambda = False

                    for i in range(0, 50):
                        
                        air_copy.gas_stream.TPX = self.orig_oxy["T"], self.orig_oxy["P"], self.orig_oxy["X"]
                        air_copy.feed_rate = N_air_in                                # Set oxygen/air mole flow as N_oxy_in mol
                        Cp_oxy1 = copy(air_copy.Cp_mole())                           # Heat capacity of air stream inlet J/mol.K
                        air_copy.gas_stream.TP = self.T_op, self.P                   # Increasing T of air stream to operation temperature
                        
                        ### Temperature of oxygen side at T_op ###
                        H_oxy_in = copy(air_copy.enthalpy_stream())
                        G_oxy_in = copy(air_copy.gibbs_stream())
                        Cp_oxy2 = copy(air_copy.Cp_mole())                  # Heat capacity of air stream internal J/mol.K
                        stream_operation( air_copy, -N_O2_an, {"O2":1.0} )             # Remove the oxygen stream which is transferred to anode
                        
                        Cp_oxy3 = copy(air_copy.Cp_mole())                  # Heat capacity of air stream outlet J/mol.K
                        H_oxy_out = copy(air_copy.enthalpy_stream())
                        G_oxy_out = copy(air_copy.gibbs_stream())
                        
                        ### Temperature is now T = T_op + delta T_max/2 ###
                        fuel_bulk.gas.TP = self.T_op + self.dT_max/2, self.P
                        air_copy.gas.TP = self.T_op + self.dT_max/2, self.P
                        Cp_f4 = copy(fuel_bulk.Cp_mole())
                        Cp_oxy4 = copy(air_copy.Cp_mole())

                        dH_stack = H_fuel_out + H_oxy_out - H_fuel_in - H_oxy_in
                        Q_anode = (Cp_f1 + Cp_f2 + Cp_f3 + Cp_f4)/4* self.dT_max * N_fuel
                        Q_Reaction = Q_rxn(H_oxy_in, H_oxy_out, W_DC)
                        Q_cathode = abs(Q_Reaction) - Q_anode
                        Cp_cat = (Cp_oxy1+ Cp_oxy2+ Cp_oxy3+ Cp_oxy4)/4

                        if not small_lambda:
                            N_air_calculated = Q_cathode / Cp_cat /(abs(T_fuel_in- T_oxy_in)+ self.dT_max)
                            air_lambda_new = N_air_calculated / N_air_stoic
                            N_air_in *= air_lambda_new/air_lambda_old
                            error = copy(abs(air_lambda_old - air_lambda_new)/air_lambda_new)
                            air_lambda_old = air_lambda_new
                        
                        else:
                            break

                        if error <= 1e-3:
                            if air_lambda_old < 1:
                                air_lambda_new = 1
                                N_air_in = N_air_stoic
                                small_lambda = True
                            else:
                                break    
                        
                    
                else:
                    # Fixed air flow
                    air_copy.gas_stream.TP = self.T_op, self.P
                    air_copy.feed_rate = self.orig_oxy["feedrate"] 
                    G_oxy_in = copy(air_copy.gibbs_stream())
                    stream_operation(air_copy, O2_produced, {"O2":1.0} )             # Transfer the oxygen produced during electrolysis to air stream
                    G_oxy_out = copy(air_copy.gibbs_stream())
                    
                delta_G = G_fuel_out - G_fuel_in + G_oxy_out - G_oxy_in              # delta Gibbs fuel cell
                V_gibbs = abs(delta_G/J_total)
                
                if V_gibbs > self.V_cell:
                    print("Switch to SOFC mode, because V_gibbs maximum is higher than cell voltage")
                    # TODO SOFC Mode
                    break
                
                else:
                    # Outlet molar fractions
                    x_H2_out = copy(fuel_copy.mole_frac('H2'))
                    x_H2O_out = copy(fuel_copy.mole_frac('H2O'))
                    x_CO_out = copy(fuel_copy.mole_frac('CO'))
                    x_CH4_out = copy(fuel_copy.mole_frac('CH4'))
                    x_O2_out = copy(air_copy.mole_frac('O2'))
                    x_H2_eq = (x_H2_in + x_H2_out)/2 + (x_CO_in + x_CO_out)/2 + 4*(x_CH4_in + x_CH4_out)/2
                    x_H2O = (x_H2O_in + x_H2O_out)/2
                    x_O2 = (x_O2_in + x_O2_out)/2
                    
                    ASR = ( self.activation_v(j) + self.diffusion_v(j, x_H2_eq, x_O2, x_H2O) + self.ohmic_v(j) )/j
                    V_op = V_gibbs - j*ASR 
                    self.V_cell = V_op
            print("\nASR: ", ASR)
            print("V_gibbs: ", V_gibbs, "\nV_op: ", V_op, "\nj[A/cm2]: ",j, "\nFU: ", FU, "\nJ: ",J_total, "\nair_lambda:")
            
                    
                
        elif self.fixed == "V,A":
            ## Fixed Voltage and fixed total cell area ##
            ## First, the maximum Gibbs voltage is calculated to determine SOFC or SOEC mode  ##

            FU_max = 1.0001
            FU_min = 0.0001
            err_tolerance = 1e-6
            f_ASR = lambda j, x_H2_eq, x_O2, x_H2O: ( self.activation_v(j) + self.diffusion_v(j, x_H2_eq, x_O2, x_H2O) + self.ohmic_v(j) )/j

            fuel_copy = GasStream(fed_to="fuel_channel",T= self.orig_fuel["T"],P= self.orig_fuel["P"],X= self.orig_fuel["X"],feed_rate= self.orig_fuel["feedrate"] )
            air_copy = GasStream(fed_to="oxygen_channel",T= self.orig_oxy["T"],P= self.orig_oxy["P"],X= self.orig_oxy["X"],feed_rate= self.orig_oxy["feedrate"] )
            
            
            def V_cell_FU(FU, converged= False):
                N_H2O_rcl = 6.00039420820532E-05   # Recycle stream of H2O after water-gas shift and reforming reactions
                ### Temperature is equal to Fuel feed temperature ###
                fuel_copy.gas_stream.TPX = self.orig_fuel["T"], self.orig_fuel["P"], self.orig_fuel["X"]
                fuel_copy.feed_rate = self.orig_fuel["feedrate"]
                air_copy.gas_stream.TPX = self.orig_oxy["T"], self.orig_oxy["P"], self.orig_oxy["X"]
                
                Cp_f1 = copy(fuel_copy.Cp_mole())           # Heat capacity fuel inlet J/mol.K
                T_fuel_in = copy(fuel_copy.gas_stream.T)
                T_oxy_in = copy(air_copy.gas.T)

                ### Temperature of fuel side is now at T_op = T_fuel_in + delta T_max/2 ###
                self.T_op = self.dT_max/2 + self.fuel_gas.T
                fuel_copy.gas_stream.TP = self.T_op, self.P

                H_fuel_in = copy(fuel_copy.enthalpy_stream())           # Enthalpy of Fuel cell inlet, before reforming
                Cp_f2 = copy(fuel_copy.Cp_mole())                       # Heat capacity fuel inlet at T_op J/mol.K
                
                fuel_copy.gas_stream.equilibrate('TP')                  # Gas reforming
                N_H2O = copy(fuel_copy.mole_flow("H2O"))
                # Molar flowrate of fuel and molarfractions inlet
                N_fuel = self.fuel_flowrate
                x_H2_in = copy(self.fuel_gas.mole_frac('H2'))
                x_H2O_in = copy(self.fuel_gas.mole_frac('H2O'))
                x_CO_in = copy(self.fuel_gas.mole_frac('CO'))
                x_CH4_in = copy(self.fuel_gas.mole_frac('CH4'))
                x_O2_in = copy(self.oxygen_gas.mole_frac('O2'))
                
                J_max = -(2*(N_H2O + N_H2O_rcl))*F                      # Max current
                j = FU*J_max/self.cell_area                             # Current density
                J_total = J_max*FU                                      # Total current
                Q1 = copy(-self.V_cell*J_total)
                
                stream_operation(fuel_copy, N_H2O_rcl, {"H2O":1})
                G_fuel_in = copy(fuel_copy.gibbs_stream())
                
                CR.react( [fuel_copy], {"H2O":-1, "H2":1, "O2":0.5}, "H2O", FU )    # Electrolysis reaction
                N_H2O_1 = copy(fuel_copy.mole_flow("H2O"))
                O2_produced = copy( fuel_copy.mole_flow("O2") )
                stream_operation(fuel_copy, - O2_produced, {"O2":1.0})              # O2 produced transfer to cathode
                
                Cp_f3 = copy(fuel_copy.Cp_mole())                           # Heat capacity in the fuel cell, after electrochemical reactions j/mol.K
                H_fuel_out = copy(fuel_copy.enthalpy_stream())              # Enthalpy of Fuel cell outlet, at T_op
                G_fuel_out = copy(fuel_copy.gibbs_stream())                 # Gibbs energy of Fuel cell outlet, at T_op
                
                fuel_copy.gas_stream.equilibrate("TP")                      # Water-gas shift reaction and other possible reforming reactions
                N_H2O_2 = copy(fuel_copy.mole_flow("H2O"))
                N_H2O_rcl = N_H2O_2 - N_H2O_1                               # Water recycle
                
                W_DC = self.V_cell* J_total
                Q_rxn = lambda H_O2_in, H_O2_out, W_DC: H_fuel_in + H_O2_in - H_fuel_out - H_O2_out - W_DC    # Reaction heat calculation [W]
                
                if self.fixed_oxygen_flow == False:
                    # Variable air flow
                    # First check the minimum oxygen-channel flowrate based on oxygen flow to fuel side
                    ### Temperature of oxygen side at T_in air stream ###
                    air_lambda_old = 1                                                    # Initial guess
                    N_air_stoic = N_O2_an / air_copy.mole_frac("O2")             # Minimum required amount of O2, stoichiometric minimum mol
                    N_air_in = N_air_stoic * air_lambda_old                               # Initial guess for air inlet based on minimum O2 need mol
                    small_lambda = False

                    for i in range(0, 50):
                        
                        air_copy.gas_stream.TPX = self.orig_oxy["T"], self.orig_oxy["P"], self.orig_oxy["X"]
                        air_copy.feed_rate = N_air_in                                # Set oxygen/air mole flow as N_oxy_in mol
                        Cp_oxy1 = copy(air_copy.Cp_mole())                           # Heat capacity of air stream inlet J/mol.K
                        air_copy.gas_stream.TP = self.T_op, self.P                   # Increasing T of air stream to operation temperature
                        
                        ### Temperature of oxygen side at T_op ###
                        H_oxy_in = copy(air_copy.enthalpy_stream())
                        G_oxy_in = copy(air_copy.gibbs_stream())
                        Cp_oxy2 = copy(air_copy.Cp_mole())                  # Heat capacity of air stream internal J/mol.K
                        stream_operation( air_copy, -N_O2_an, {"O2":1.0} )             # Remove the oxygen stream which is transferred to anode
                        
                        Cp_oxy3 = copy(air_copy.Cp_mole())                  # Heat capacity of air stream outlet J/mol.K
                        H_oxy_out = copy(air_copy.enthalpy_stream())
                        G_oxy_out = copy(air_copy.gibbs_stream())
                        
                        ### Temperature is now T = T_op + delta T_max/2 ###
                        fuel_bulk.gas.TP = self.T_op + self.dT_max/2, self.P
                        air_copy.gas.TP = self.T_op + self.dT_max/2, self.P
                        Cp_f4 = copy(fuel_bulk.Cp_mole())
                        Cp_oxy4 = copy(air_copy.Cp_mole())

                        dH_stack = H_fuel_out + H_oxy_out - H_fuel_in - H_oxy_in
                        Q_anode = (Cp_f1 + Cp_f2 + Cp_f3 + Cp_f4)/4* self.dT_max * N_fuel
                        Q_Reaction = Q_rxn(H_oxy_in, H_oxy_out, W_DC)
                        Q_cathode = abs(Q_Reaction) - Q_anode
                        Cp_cat = (Cp_oxy1+ Cp_oxy2+ Cp_oxy3+ Cp_oxy4)/4

                        if not small_lambda:
                            N_air_calculated = Q_cathode / Cp_cat /(abs(T_fuel_in- T_oxy_in)+ self.dT_max)
                            air_lambda_new = N_air_calculated / N_air_stoic
                            N_air_in *= air_lambda_new/air_lambda_old
                            error = copy(abs(air_lambda_old - air_lambda_new)/air_lambda_new)
                            air_lambda_old = air_lambda_new
                        
                        else:
                            break

                        if error <= 1e-3:
                            if air_lambda_old < 1:
                                air_lambda_new = 1
                                N_air_in = N_air_stoic
                                small_lambda = True
                            else:
                                break    
                        
                    
                else:
                    # Fixed air flow
                    air_copy.gas_stream.TP = self.T_op, self.P
                    air_copy.feed_rate = self.orig_oxy["feedrate"] 
                    G_oxy_in = copy(air_copy.gibbs_stream())
                    stream_operation(air_copy, O2_produced, {"O2":1.0} )             # Transfer the oxygen produced during electrolysis to air stream
                    G_oxy_out = copy(air_copy.gibbs_stream())
                    
                delta_G = G_fuel_out - G_fuel_in + G_oxy_out - G_oxy_in              # delta Gibbs fuel cell
                V_gibbs = abs(delta_G/J_total)
                
                if V_gibbs > self.V_cell:
                    print("Switch to SOFC mode, because V_gibbs maximum is higher than cell voltage")
                    # TODO SOFC Mode
                
                else:
                    # Outlet molar fractions
                    x_H2_out = copy(fuel_copy.mole_frac('H2'))
                    x_H2O_out = copy(fuel_copy.mole_frac('H2O'))
                    x_CO_out = copy(fuel_copy.mole_frac('CO'))
                    x_CH4_out = copy(fuel_copy.mole_frac('CH4'))
                    x_O2_out = copy(air_copy.mole_frac('O2'))
                    x_H2_eq = (x_H2_in + x_H2_out)/2 + (x_CO_in + x_CO_out)/2 + 4*(x_CH4_in + x_CH4_out)/2
                    x_H2O = (x_H2O_in + x_H2O_out)/2
                    x_O2 = (x_O2_in + x_O2_out)/2
                    
                    ASR = ( self.activation_v(j) + self.diffusion_v(j, x_H2_eq, x_O2, x_H2O) + self.ohmic_v(j) )/j
                    V_op = V_gibbs - j*ASR
                
                if converged:
                    print("\nASR: ", ASR)
                    print("V_gibbs: ", V_gibbs, "\nV_op: ", V_op, "\nj[A/cm2]: ",j, "\nFU: ", FU, "\nJ: ",J_total, "\nair_lambda:")
                else:
                    return V_op - self.V_cell
            
            a = FU_min
            b = FU_max
            an = a
            bn = b
            
            for i in range(0, 100):
                
                mn = (an + bn)*0.5
                f_mn = V_cell_FU(mn)
                if V_cell_FU(an)*f_mn < 0:
                    bn = mn

                elif V_cell_FU(bn)*f_mn < 0:
                    an = mn

                if abs(f_mn) < err_tolerance:
                    FU = (an+bn)*0.5
                    V_cell_FU(FU, converged= True)
                    break                 


    def D_Knudsen(self, gas, i):
        """ 
            gas(GasStream): Which gas?
            i (str): Component name i
            >>> D_Knudsen(self.fuel_gas, "H2")
        """
        
        Mi = gas.MW(i)
        #print("dknudsen", i, self.pore_diameter/3* (8*R*self.T_op/pi/Mi)**0.5 )
        return self.pore_diameter/3* (8*R*self.T_op/pi/Mi)**0.5

    
    def D_binary(self, gas, i, k):
        
        Mi = gas.MW(i)
        Mk = gas.MW(k)
        Vd_i = self.Vd[ i ]
        Vd_k = self.Vd[ k ]
        P = self.P*1e-5         # convert from Pa to bar
        #print("dbin", i,k, 1.43e-7* self.T_op**1.75 / P / ( 2/(1/Mi+1/Mk) )**0.5 / (Vd_i**(1/3) + Vd_k**(1/3))**2 )
        return 1.43e-7* self.T_op**1.75 / P / ( 2/(1/Mi+1/Mk) )**0.5 / (Vd_i**(1/3) + Vd_k**(1/3))**2


    def D_eff(self, gas, i, k):

        D_Knudsen = self.D_Knudsen(gas, i)
        D_binary = self.D_binary(gas, i, k)
        #print("deff", i,k, self.porosity/self.tortuosity* D_Knudsen* D_binary/(D_Knudsen + D_binary) )
        return self.porosity/self.tortuosity* D_Knudsen* D_binary/(D_Knudsen + D_binary)


    def tpb_H2(self, j, x_H2):

        P = self.P 
        P_bulk = x_H2 * P
        D_eff = self.D_eff(self.fuel_gas, "H2", "H2O")
        
        return P_bulk - R*self.T_op*self.thickness_fuel* j*1e+4 /self.n_fuel/F/D_eff


    def tpb_H2O(self, j, x_H2O):

        P = self.P 
        P_bulk = x_H2O * P
        D_eff = self.D_eff(self.fuel_gas, "H2O", "H2")
        
        return P_bulk + R*self.T_op*self.thickness_fuel* j*1e+4 /(self.n_fuel+2)/F/D_eff


    def tpb_O2(self, j, x_O2):

        P = self.P 
        P_bulk = x_O2 * P
        D_eff = self.D_eff(self.oxygen_gas, "O2", "N2")
        
        return P_bulk - R*self.T_op*self.thickness_oxy* (j*1e+4) /self.n_oxy / F / D_eff


    def diffusion_v(self, j, x_H2, x_O2, x_H2O):

        P_H2_bulk = x_H2*self.P
        P_H2O_bulk = x_H2O*self.P
        P_H2_tpb = self.tpb_H2(j, x_H2)
        P_H2O_tpb = self.tpb_H2O(j, x_H2O)
        P_O2_bulk = x_O2* self.P
        P_O2_tpb = self.tpb_O2(j, x_O2)

        if self.operation == "SOFC":
            ovp_fuel =  R*self.T_op* log( P_H2O_tpb* P_H2_bulk/ P_H2O_bulk/ P_H2_tpb )/(2*F)
        
        else:
            ovp_fuel =  R*self.T_op* log( P_H2O_bulk* P_H2_tpb / P_H2O_tpb / P_H2_bulk )/(2*F)
        ovp_oxy =   R*self.T_op* log( (P_O2_bulk/P_O2_tpb)**0.5 )/(2*F)
        #print("Loss diffusion_fuel: ", ovp_fuel)
        #print("Loss diffusion_oxy: ", ovp_oxy)
        return ovp_fuel+ovp_oxy

    
    def activation_v(self, j, alpha= 0.5):
        
        j0_fuel = self.gamma_fuel* exp( -self.E_act_fuel/ R/ self.T_op )
        j0_oxy = self.gamma_oxy* exp( -self.E_act_oxy/ R/ self.T_op )
        
        act_fuel = R*self.T_op/(F* self.n_fuel*alpha)* asinh( j/ (2 * j0_fuel) )
        act_oxy = R*self.T_op/(F* self.n_oxy*alpha)* asinh( j/ (2 *j0_oxy) )
        #print("Loss act_fuel: ", act_fuel)
        #print("Loss act_oxy: ", act_oxy)
        return act_fuel+act_oxy

    
    def ohmic_v(self,j):

        sigma_el = self.sigma0_el*exp(-self.E_act_el/R/self.T_op)
        r_ohmic_el = self.thickness_el*100/sigma_el
        ovp = j*(r_ohmic_el + self.r_ohmic_const)
        #print("Loss ohmic_v: ", ovp)
        return ovp

    def V_J_profile(self):

        kazempoor_5050_850_j = [-1.2820741135490246,-1.1631944183952885,-0.9961187979475132,-0.8322385789362716,-0.6458332632700867,-0.4690762387289791,-0.3083748646701441,
        -0.15732592085741892,-0.04160851417196221,0.05802426818399953,0.21873805960593007,0.5594622246999461,0.8455210183272519,1.0801677896188262,1.2730177167314918,1.4593816411873592]
        kazempoor_5050_850_V = [1.2550751444488566,1.219230356313647,1.1693060892937754,1.1258116881017086,1.0862145280561266,1.0466018463066753,1.014685342725689,0.9814663094950143,
        0.9623474507541608,0.9419157145610398,0.9138602348175101,0.857780318738189,0.8041864622285805,0.7646669107023443,0.7289411222967983,0.6764738827930272]
        kazempoor_25_750_j = [0.048342864, 0.141444113, 0.224913628, 0.321197862,0.414216329]
        kazempoor_25_750_V = [0.931604129, 0.880273855, 0.83407609, 0.772454926, 0.695384492]
        kazempoor_25_H2_850_j = [0.09822710055304751,0.22009608751951393,0.36848604906402627,0.5433493230959134,0.7155703299761969,0.9196087612604023,1.107704223214976,1.2692787105915035,1.3962922244698295]
        kazempoor_25_H2_850_V = [0.9391302040812527,0.9130970212124231,0.8897679172690367,0.8547608563909558,0.822449682516484,0.7910463225252715,0.7488519673179553,0.7039535331851962,0.6302813459699933]

        if self.operation =="SOFC" and self.fixed == "FU,A":
            j = []
            V = []
            for FU in range(1, 99, 3):
                self.FU = FU / 100
                #try:
                Vi, ji, ASRi = self.SOFC_mode()
                j.append(ji)
                V.append(Vi)
                #except:
                    #print("Error: FU = ", FU, " couldn't be calculated.")
        
            
            plt.plot(j, V)
            plt.plot(kazempoor_25_H2_850_j, kazempoor_25_H2_850_V, "*")
            plt.show()

# Gas stream definition for anode and cathode. See gas.py for details. "H2:0.21, H2O:0.50, N2:0.125, CO:0.07, CH4:0.025 CO2:0.07"
gas_file = 'gri30.yaml'
anode_gas = GasStream (gas_file, "fuel_channel", T = 850+273, feed_rate= 0.01859/60, X= ":CO2:0.25, H2O:0.25, CO:0.25, H2:0.25") # H2:0.12, H2O:0.57, CH4: 0.28, CO:0.03
cathode_gas = GasStream (gas_file, "oxygen_channel", T = 850+273, feed_rate= 0.02/60, X= "O2:1.0")
anode_gas.change_feed_unit("mol")   # Call this if the feedrate as "kmol/s" or "kg/s" given. 
cathode_gas.change_feed_unit("mol")
rsofc = SOFC_0D(anode_gas, cathode_gas, operation="SOEC", fixed="FU,A")
rsofc.fixed_oxygen_flow = True
rsofc.V_J_profile()



