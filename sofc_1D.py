__author__ = "Burak Polat"
__email__ = "burakplt@outlook.com"

import os
import sys
import time
import cantera as ct
from math import log, asinh, exp, pi
from copy import copy, deepcopy
from matplotlib import pyplot as plt
from gas import *
from scipy.optimize import least_squares, fsolve, curve_fit, minimize
import numpy as np


R = 8.3145
F = 96485.3329

class RSOFC_1D():

    """ Reversible SOFC Model, 1-Dimensional.

        Args:
            fuel_gas(GasStream): The gas stream is fed to "anode"
            oxygen_gas(GasStream): The gas stream is fed to "cathode"
            operation(str): Operation mode 'SOFC' or 'SOEC'
            segment_number(int): How many computational segments will be created?
            
    """

    def __init__(self, gas_file, fuel_gas, air, segment_N, operation = "SOFC"):

        self.name = "RSOFC Reactor 1-D"

        ### Gas Feeds ###
        self.gas_file = gas_file
        self.fuel_gas = fuel_gas
        self.air = air
        self.orig_oxy = {"T":copy(self.air.T), "P":copy(self.air.P), "X":copy(self.air.gas_stream.X),"feedrate":copy(self.air.feed_rate) }
        self.orig_fuel = {"T":copy(self.fuel_gas.T), "P":copy(self.fuel_gas.P), "X":copy(self.fuel_gas.gas_stream.X),"feedrate":copy(self.fuel_gas.feed_rate) }
        
        ### Operational Parameters ###
        self.operation = operation      # SOFC or SOEC
        #self.fixed = fixed              # Varying "FU,V", "FU,A", or "V,A"
        self.flow = "co"                # "co" or "counter"
        self.N = segment_N
        self.cell_area = 100             # cm2
        self.tpb_area_vol = 7e+5         # Active Surface area per catalyst volume at TPB 
        self.w_ch = 0.3                 # Channel width cm
        self.l_ch = 10                  # Channel length cm
        self.h_ch = 0.1                  # Channel height cm
        self.channel_number = 33        # Number of channels
        self.FU = 0.80                   # Fuel utilisation
        self.V_cell = 0.85               # Volt
        self.ne_fuel = 2                 # number of electrons for fuel
        self.ne_air = 4                  # number of electrons for air
        self.P = self.fuel_gas.P        # Pressure [Pa]
        
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
        self.thickness_fuel = 3.2e-5    # m
        self.thickness_oxy = 1.75e-5    # m
        self.k_solid = 2                # W/mK
        self.Vd = {"H2": 6.12, "N2": 18.5, "O2": 16.3, "H2O": 13.1} # Diffusion volumes
        
        ### Ohmic loss parameters ###
        self.sigma0_el = 333.3       # 1/ohm.cm2 
        self.r_ohmic_const = 0.057   # ohm.cm2
        self.thickness_el = 1.25e-5  # m


    def run_simulation(self):
        
        sys.stdout = open("log.txt", 'w')
        self.set_initial_conditions()
        print("Initial conditions set.")

        #solution = least_squares (fun= self.residual_function,x0= self.initial_guess, bounds=(self.lower_bounds, self.upper_bounds), xtol=1e-12, ftol=1e-5, diff_step=1e-2)
        #print(solution)
        self.residual_function(self.initial_guess)

    def set_initial_conditions (self):
        
        self.As = self.w_ch* self.l_ch / self.N                                              # Heat transfer area of a segment
        self.Ai = self.tpb_area_vol *self.thickness_fuel* self.w_ch* self.l_ch/ self.N       # Electrochemically active area of a segment cm2
        self.vol_segment = self.As*self.h_ch
        self.A_lumped = self.w_ch*(self.thickness_el+self.thickness_fuel+self.thickness_oxy)*1.2    # TODO interconnect usw.
        B = self.h_ch / self.w_ch
        Dh = 2*self.w_ch*self.h_ch / (self.h_ch + self.w_ch)
        self.Nu_Dh = (7.541 -2.61*B + 4.97*B**2 -5.119*B**3 +2.702*B**4 +0.548*B**5)/Dh

        j_avg = 0.4     # TODO get the value by calculating with 0D Model

        Tf_func = lambda x : -2083.3*x**5 + 6796.3*x**4 - 8084.9*x**3 + 4048.6*x**2 - 567.14*x + self.orig_fuel["T"]
        Ts_func = lambda x : -1362.2*x**5 + 4731.2*x**4 - 5898.5*x**3 + 3004.6*x**2 - 350.7*x + self.orig_fuel["T"]-15
        Ta_func = lambda x : 448.72*x**5 - 882.87*x**4 + 232.23*x**3 + 455.33*x**2 - 208.84*x + self.orig_oxy["T"]
        j_func = lambda x : 11.538*x**5 - 34.441*x**4 + 38.234*x**3 - 19.096*x**2 + 3.7606*x + j_avg - 0.1
        
        J = []
        Tf = []
        Ts = []
        Ta = []
        low = {"J":[], "Tf":[],"Ts":[],"Ta":[]}
        up = {"J":[], "Tf":[],"Ts":[],"Ta":[]}

        for i in range (0, self.N):
            J.append( j_func(i/(self.N-1)) )
            Tf.append( Tf_func(i/(self.N-1)) )
            Ts.append( Ts_func(i/(self.N-1)) )
            Ta.append( Ta_func(i/(self.N-1)) )
            low["J"].append(0)
            low["Tf"].append(Tf[i]-100)
            low["Ta"].append(Ta[i]-100)
            low["Ts"].append(Ts[i]-100)
            up["J"].append(1.2)
            up["Tf"].append(Tf[i]+100)
            up["Ts"].append(Ts[i]+100)
            up["Ta"].append(Ta[i]+100)
            
        self.lower_bounds = low["J"] #+ low["Tf"] + low["Ts"] + low["Ta"]
        self.upper_bounds = up["J"] #+ up["Tf"] + up["Ts"] + up["Ta"]
        self.initial_guess = J #+Ts+Tf+Ta
        self.temps = J+Ts+Tf+Ta
        #self.initial_guess.append( self.V_cell )

    def residual_function (self, vars):

        # vars =  J J..Ts Ts...Tf Tf...Ta Ta...V
        J = vars #[0:self.N]
        
        Ts = self.temps[self.N:2*self.N]
        Tf = self.temps[2*self.N:3*self.N]
        Ta = self.temps[3*self.N:4*self.N]
        V = self.V_cell 

        
        # Reset fuel and air gas streams
        fg = self.fuel_gas
        fuel = self.fuel_gas.gas_stream
        fuel.TPX = self.orig_fuel["T"], self.P, self.orig_fuel["X"]
        fg.feed_rate = self.orig_fuel["feedrate"]
        air = self.air.gas
        airg = self.air.gas_stream
        airg.TPX = self.orig_oxy["T"], self.P, self.orig_oxy["X"]
        self.air.feed_rate = self.orig_oxy["feedrate"]
        #surface = self.fuel_gas.surface

        # Inlet flows cantera.Reservoir
        fuel_inlet = ct.Reservoir(fg.gas)
        air_inlet = ct.Reservoir(air)
        fuel_outlet = ct.Reservoir(fg.gas)
        
        # cantera.Reactor objects to simulate anode and cathode channels
        anode = ct.IdealGasReactor(fg.gas, energy='off')
        anode.volume = self.vol_segment
        #r_surface = ct.ReactorSurface(surface, anode, A = self.Ai)
        
        # Mass flowrates of channels cantera.MassFlowController
        mfc_fuel = ct.MassFlowController(fuel_inlet, anode)
        #mfc_O2 = ct.MassFlowController(air_inlet, anode)                # O2 flow to anode reaction zone
        mfc_fuel.mass_flow_rate = self.fuel_gas.mass_flow()
        
        N_O2_an = J[0]*self.Ai/4/F  # mol/s
        O2_an = GasStream(self.gas_file, "oxygen_channel", self.orig_oxy["T"], self.P, X ="O2:1.0", feed_rate= N_O2_an )     # Required for energy calculations
        #mfc_O2.mass_flow_rate = O2_an.mass_flow()

        # Outlet Valve (required for cantera.Reactor)
        v1 = ct.Valve(anode, fuel_outlet, K=1e-5)
        
        # Cantera reactor network, represents a complete SOFC cell 
        cell = ct.ReactorNet([anode])
        cell.max_err_test_fails = 12
        cell.rtol = 1.0e-8
        cell.atol = 1.0e-8
        #cell.reinitialize()

        # Residual vectors
        res_ele = []
        res_en_fuel = []
        res_en_air = []
        res_en_tpb = []
        
        thermal_conductivity_f = thermal_conductivity({ "H2":fg.mole_frac("H2"), "CH4":fg.mole_frac("CH4"), "H2O":fg.mole_frac("H2O"), "CO":fg.mole_frac("CO"), "CO2":fg.mole_frac("CO2") }, Tf[0])
        thermal_conductivity_a = thermal_conductivity({ "O2": self.air.mole_frac("O2"), "N2":self.air.mole_frac("N2")}, Ta[0]  )
        #I_max = (2*self.fuel_gas.mass_frac("H2") + 8*self.fuel_gas.mass_frac("CH4") + 2*self.fuel_gas.mass_frac("CO"))*F*self.orig_fuel["feedrate"]
        #I_SOFC = I_max * self.FU
        
        # Iteration through discrete cell segments 
        for i in range(0, self.N):
            
            O2_an.feed_rate = J[i]*self.Ai/8/F                                              # O2 flow to TPB mol/s
            H_O2_tpb = copy( O2_an.enthalpy_stream() )                                      # Enthalpy of O2 stream flowing into TPB Watt
            fuel.TP = Tf[i], self.P
            H_fuel_in = copy( fg.enthalpy_stream() )                                        # Enthalpy of fuel inlet at T = Tf
            x_H2_eq_in = copy( fg.mole_frac("H2") + fg.mole_frac("CO") + 4*fg.mole_frac("CH4") )
            x_H2_in = copy( fg.mole_frac("H2") )
            G_fuel_in = copy( fg.gibbs_stream() )
            fuel.TDY = Ts[i], anode.thermo.density, anode.thermo.Y                          # Temperature should be "Tf", but we need to set surface temperature as "Ts". Surface temp. is always same as the bulk gas temperature. 
            airg.TP = Ta[i], self.P
            H_air_in = copy( self.air.enthalpy_stream() )                                   # Enthalpy of air inlet at T = Ta
            G_air_in = copy( self.air.gibbs_stream() )
            x_O2_in = copy( self.air.mole_frac("O2") )

            stream_operation(fg, O2_an.feed_rate, {"O2":1.0} )
            x_O2_an = copy( fg.mole_frac("O2") )
            fg.gas()
            # Cantera Reactor calculations
            fuel_inlet.syncState()
            air_inlet.syncState()
            anode.syncState()
            cell.reinitialize()
            cell.advance_to_steady_state()
            
            fuel.TPX = Tf[i], self.P, anode.thermo.X
            fg.feed_rate += O2_an.feed_rate
            H_fuel_out = copy( fg.enthalpy_stream() )                                # Enthalpy of fuel outlet at T = Tf with outlet mass flowrate
            G_fuel_out = copy( fg.gibbs_stream() )
            stream_operation(self.air, -O2_an.feed_rate, {"O2":1.0} )                # Update cathode air stream by removing O2_to_anode
            H_air_out = copy( self.air.enthalpy_stream() )                           # Enthalpy of air outlet at T = Ta with outlet mass flowrate
            G_air_out = copy( self.air.gibbs_stream() )
  
            # Electrochemical balance
            x_H2_eq_out = fg.mole_frac("H2") + fg.mole_frac("CO") + 4*fg.mole_frac("CH4")
            delta_G = G_fuel_out - G_fuel_in + G_air_out - G_air_in
            ASR = self.calculate_ASR( J[i], (x_H2_in +fg.mole_frac("H2"))/2, (self.air.mole_frac("O2")+x_O2_in)/2, fg.mole_frac("H2O"), (x_H2_eq_in+x_H2_eq_out)/2, Ts[i] )
            res_ele.append( abs(delta_G) / self.Ai/J[i] - ASR*J[i] - V )
            print({ "H2":fg.mole_frac("H2"), "CH4":fg.mole_frac("CH4"), "H2O":fg.mole_frac("H2O"), "CO":fg.mole_frac("CO"), "CO2":fg.mole_frac("CO2") })
            # Energy balance
            #E_react = surface.heat_relase_rate*anode.volume         # Watt 
            """
            #x_f = { "H2":fg.mole_frac("H2"), "CH4":fg.mole_frac("CH4"), "H2O":fg.mole_frac("H2O"), "CO":fg.mole_frac("CO"), "CO2":fg.mole_frac("CO2") }
            hf_f = self.Nu_Dh * thermal_conductivity_f

            Q_conv_sf = hf_f * self.As * (Ts[i] - Tf[i])
            res_en_fuel.append( H_fuel_in + Q_conv_sf - H_fuel_out )
            #x_a = { "O2": (self.air.mole_frac("O2")+x_O2_in)/2, "N2":self.air.mole_frac("N2") } 
            hf_a = self.Nu_Dh * thermal_conductivity_a
            
            Q_conv_sa = hf_a * self.As * (Ts[i] - Ta[i])
            res_en_air.append( H_air_in - H_O2_tpb + Q_conv_sa - H_air_out )
            
            W_DC = J[i]*self.Ai*V
            
            if i == 0:
                Q_cond = self.k_solid * self.l_ch/self.N * self.A_lumped*( Ts[i]-Ts[i+1] )
            elif i == self.N-1:
                Q_cond = self.k_solid * self.l_ch/self.N * self.A_lumped*( Ts[i-1]-Ts[i] )
            else:
                Q_cond = self.k_solid * self.l_ch/self.N * self.A_lumped*( Ts[i-1]-Ts[i+1] )
            
            res_en_tpb.append( Q_cond + H_O2_tpb - Q_conv_sa - Q_conv_sf - W_DC)
            
            
        residuals = res_ele + res_en_fuel + res_en_air + res_en_tpb
        """
        #print("\nparams\n",J, Tf, Ts, Ta)
        #print("\nresiduals\n", residuals)
        #residuals.append( I_SOFC/self.cell_area - sum(J)/self.N )
        
        residuals = res_ele
        #print(vars, residuals)
        return residuals
    

    def calculate_ASR (self, j, x_H2, x_O2, x_H2O, x_H2_eq, T_s ):

        def D_Knudsen(gas, i):
            """
                i (str): Component name i
            """
            
            Mi = gas.MW(i)
            #print("dknudsen", i, self.pore_diameter/3* (8*R*T_s/pi/Mi)**0.5 )
            return self.pore_diameter/3* (8*R*T_s/pi/Mi)**0.5

        
        def D_binary(gas, i, k):
            
            Mi = gas.MW(i)
            Mk = gas.MW(k)
            Vd_i = self.Vd[ i ]
            Vd_k = self.Vd[ k ]
            P = self.P*1e-5         # convert from Pa to bar
            #print("dbin", i,k, 1.43e-7* T_s**1.75 / P / ( 2/(1/Mi+1/Mk) )**0.5 / (Vd_i**(1/3) + Vd_k**(1/3))**2 )
            return 1.43e-7* T_s**1.75 / P / ( 2/(1/Mi+1/Mk) )**0.5 / (Vd_i**(1/3) + Vd_k**(1/3))**2


        def D_eff(gas, i, k):

            DKnudsen = D_Knudsen(gas, i)
            Dbinary = D_binary(gas, i, k)
            #print("deff", i,k, self.porosity/self.tortuosity* D_Knudsen* D_binary/(D_Knudsen + D_binary) )
            return self.porosity/self.tortuosity* DKnudsen* Dbinary/(DKnudsen + Dbinary)


        def tpb_H2(j, x_H2):

            P = self.P 
            P_bulk = x_H2 * P
            Deff = D_eff(self.fuel_gas, "H2", "H2O")
            tpb = P_bulk - R*T_s*self.thickness_fuel* j*1e+4 /self.ne_fuel/F/Deff
            if tpb<0:
                tpb = 1e-5
            return tpb


        def tpb_H2O(j, x_H2O):

            P = self.P 
            P_bulk = x_H2O * P
            Deff = D_eff(self.fuel_gas, "H2O", "H2")
            
            return P_bulk + R*T_s*self.thickness_fuel* j*1e+4 /self.ne_fuel/F/Deff


        def tpb_O2(j, x_O2):

            P = self.P 
            P_bulk = x_O2 * P
            Deff = D_eff(self.air, "O2", "N2")
            
            return P_bulk - R*T_s*self.thickness_oxy* (j*1e+4) /self.ne_air / F / Deff


        def diffusion_v(j, x_H2, x_O2, x_H2O):

            P_H2_bulk = x_H2*self.P
            P_H2O_bulk = x_H2O*self.P
            P_H2_tpb = tpb_H2(j, x_H2)
            P_H2O_tpb = tpb_H2O(j, x_H2O)
            P_O2_bulk = x_O2* self.P
            P_O2_tpb = tpb_O2(j, x_O2)
            
            if self.operation == "SOFC":
                ovp_fuel =  R*T_s* log( P_H2O_tpb* P_H2_bulk/ P_H2O_bulk/ P_H2_tpb )/(2*F)
            
            else:
                ovp_fuel =  R*T_s* log( P_H2O_bulk* P_H2_tpb / P_H2O_tpb / P_H2_bulk )/(2*F)
            ovp_oxy =   R*T_s* log( (P_O2_bulk/P_O2_tpb)**0.5 )/(2*F)
            #print("Loss diffusion_fuel: ", ovp_fuel)
            #print("Loss diffusion_oxy: ", ovp_oxy)
            return ovp_fuel+ovp_oxy

        
        def activation_v(j, alpha= 0.5):
            
            j0_fuel = self.gamma_fuel* exp( -self.E_act_fuel/ R/ T_s )
            j0_oxy = self.gamma_oxy* exp( -self.E_act_oxy/ R/ T_s )
            
            act_fuel = R*T_s/(F* self.ne_fuel*alpha)* asinh( j/ (2 * j0_fuel) )
            act_oxy = R*T_s/(F* self.ne_air*alpha)* asinh( j/ (2 *j0_oxy) )
            #print("Loss act_fuel: ", act_fuel)
            #print("Loss act_oxy: ", act_oxy)
            return act_fuel+act_oxy

        
        def ohmic_v(j):

            sigma_el = self.sigma0_el*exp(-self.E_act_el/R/T_s)
            r_ohmic_el = self.thickness_el*100/sigma_el
            ovp = j*(r_ohmic_el + self.r_ohmic_const)
            #print("Loss ohmic_v: ", ovp)
            return ovp


        ASR = ( activation_v(j) + diffusion_v(j, x_H2_eq, x_O2, x_H2O) + ohmic_v(j) )/j
        return ASR


gas_file =  'surface1-gas.yaml'
anode_gas = GasStream(gas_file, "fuel_channel", T = 1073, feed_rate= 1.32e-7, feed_unit="kg", X="H2:0.12, H2O:0.57, CH4: 0.28, CO:0.03") # H2:0.12, H2O:0.57, CH4: 0.28, CO:0.03
cathode_gas = GasStream("air.yaml", "oxygen_channel", T = 1073, feed_rate= 4.5e-6, feed_unit="kg", X="O2:0.21, N2:0.79")
anode_gas.add_surface(gas_file, 'YSZ', [anode_gas.gas])

anode_gas.change_feed_unit("mol")   # Call this if the feedrate as "kmol/s" or "kg/s" given. 
cathode_gas.change_feed_unit("mol")
SOFC = RSOFC_1D(gas_file, anode_gas, cathode_gas, 5)
SOFC.run_simulation()