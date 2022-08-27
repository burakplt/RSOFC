__author__ = "Burak Polat"
__email__ = "burakplt@outlook.com"

import cantera as ct
from copy import copy
from math import log
from matplotlib import pyplot as plt

ct.add_module_directory()

GAS_DEFAULT = {"T": 1000, "P": 101325, "fuel_channel": "H2:0.50, H2O:0.50", "oxygen_channel": "O2:0.21, N2:0.79", "feed_rate": 1e-4}


def stream_operation(stream, amount, composition):
    """ Removing desired part of a stream (Quantity object). Cantera allows mixing two strams,
        but it is not possible to subtract or add a portion of a stream by default.
        With this function these operations can be done. 
        If amount is negative, two streams will be seperated, Otherwise will be combined.
        
        Args:
            stream(GasStream): The gas stream
            amount(float): Molar flowrate [mole/s] of the part that will be removed from/added to stream
            composition(dict): Molar composition of the part to be removed.
        Example:
        >>> stream_operation(fuel_gas, -0.05, {"O2":0.9, "N2":0.1})
    """

    if amount != 0:
        unit = copy(stream.feed_unit)
        stream.change_feed_unit("mol")
        flowrate = stream.feed_rate                         # Current flowrate [mole/s]
        x_mole = stream.gas_stream.mole_fraction_dict()     # Current mole fractions of components in stream

        flowrates = {}
        for component, frac in x_mole.items():
            flowrates[component] = frac*flowrate

        for component, frac in composition.items():
            try:
                current_mole = x_mole[component] * flowrate     # Current mole of the component
            except KeyError:
                current_mole = 0
                flowrates[component] = 0
            
            change = frac*amount                            # Amount of this component to be removed/added
            
            if change < 0 and abs(change) > current_mole:
                change = -current_mole                       # If the amount to remove greater than current mole, remove all of it
            
            flowrates[component] += change

        new_fractions = {}
        new_flowrate = flowrate + amount

        for component, N in flowrates.items():
            new_fractions[component] = N/new_flowrate

        stream.gas_stream.X = new_fractions
        stream.feed_rate = new_flowrate
        stream.change_feed_unit(unit)
    

class GasStream ():
    
    """ Gas object based on Cantera gri30 gas model. 
        Args:
            fed_to(str): This gas stream is fed to "fuel_channel" or "oxygen_channel"
            T(float): Temperature of gas stream in Kelvin. If empty, default value 1000 K.
            P(float): Pressure of gas stream in Pascal. If empty, default 101325 Pa.
            X(str): Molar composition of stream. 
                    Default "H2:0.50, H2O:0.50" for anode, "O2:0.21, N2:0.79" for cathode
            feed_rate(float): Molar flow rate of stream. If empty, default 0.0001 moles/second.
            feed_unit(str): Unit of the flowrate. "mol", "kmol" or "kg" per second.

    """
    
    def __init__(self, gas_file="gri30.yaml", fed_to="fuel_channel", T = GAS_DEFAULT["T"], P = GAS_DEFAULT["P"], X = None, feed_rate = GAS_DEFAULT["feed_rate"], feed_unit="mol", flash_type = 'TP'):

        self.gas = ct.Solution(gas_file)
        self.fed_to = fed_to  # Fuel or air
        self.T = T
        self.P = P
        self.feed_unit = feed_unit
        self.flash_type = flash_type

        if X is not None:
            self.gas.TPX = T, P, X
        else:
            if fed_to == "fuel_channel":
                self.gas.TPX = T, P, GAS_DEFAULT["fuel_channel"]
            elif fed_to == "oxygen_channel":
                self.gas.TPX = T, P, GAS_DEFAULT["oxygen_channel"]
        
        if feed_unit == "kg":
            self.gas_stream = ct.Quantity(self.gas, moles = feed_rate*self.gas.mean_molecular_weight, constant= flash_type)
            self.feed_rate = feed_rate    # kg/s

        elif feed_unit == "kmol":
            self.gas_stream = ct.Quantity(self.gas, moles = feed_rate, constant= flash_type)
            self.feed_rate = feed_rate # kmol/s
        
        else:
            self.gas_stream = ct.Quantity(self.gas, moles = feed_rate/1000, constant= flash_type)
            self.feed_rate = feed_rate  # mol/s
        
        
        # Gas stream with a defined molar flow rate in mole/s. Cantera.Quantity object
        # Gas stream is equivalent to a material stream in Aspen. It has a flowrate and flash type like 'TP' or 'HS' etc. 
        # In Cantera the mole flow of a stream is represented in [kmol], mass flow in [kg]. But we define a mole flow in fuel cell in [moles/s]
        # ... that's why the feed_rate is divided to 1000 when defining a gas_stream.


    @property
    def feed_rate(self):
        
        if self.feed_unit == "kg":
            return self.gas_stream.mass

        elif self.feed_unit == "kmol":
            return self.gas_stream.moles
        
        else:
            return self.gas_stream.moles*1000
    
    @feed_rate.setter
    def feed_rate(self, new_flowrate):

        if self.feed_unit == "kg":
            self._feed_rate = new_flowrate
            self.gas_stream.mass = self._feed_rate
            
        elif self.feed_unit == "kmol":
            self._feed_rate = new_flowrate
            self.gas_stream.moles = self._feed_rate
        
        else:
            self._feed_rate = new_flowrate
            self.gas_stream.moles = self._feed_rate/1000
    
    
    def change_feed_unit(self, new_unit):
        if new_unit != self.feed_unit:
            feed_rate = self.gas_stream.moles       #kmol
            self.feed_unit = new_unit
            
            if self.feed_unit == "kg":
                self.feed_rate = feed_rate * self.gas.mean_molecular_weight
                
            elif self.feed_unit == "kmol":
                self.feed_rate = feed_rate
            
            else:
                self.feed_rate = feed_rate*1000
        
    
    def add_surface(self, surf_file, phase_name, gases ):
        self.surface = ct.Interface(surf_file, phase_name, gases)

    def enthalpy_stream(self):
        return self.gas_stream.enthalpy   #Total enthalpy in J

    def entropy_stream(self):
        return self.gas_stream.entropy  # Total entropy in J/K

    def gibbs_stream(self):
        return self.gas_stream.gibbs    # Gibbs free energy in J

    def Cp_mole(self):
        return self.gas_stream.cp_mole/1000    # Heat capacity J/mol.K


    def MW (self, component):
        """Returns the Molecular Weight of a given component
        Args:
            component(str): Formula of the component, i.e "H2" for hydrogen
        """
        
        index = self.gas.species_index(component)
        return self.gas.molecular_weights[index]


    def mole_frac(self, component):
        """Returns mole fraction of a given component
        Args:
            component(str): Formula of the component, i.e "H2" for hydrogen
        """

        index = self.gas_stream.species_index(component)
        return self.gas_stream.X[index]
    

    def mass_frac(self, component):
        """Returns mass fraction of a given component
        Args:
            component(str): Formula of the component, i.e "H2" for hydrogen
        """

        index = self.gas_stream.species_index(component)
        return self.gas_stream.Y[index]

    
    def mole_flow(self, component=None):
        """Returns molar flowrate of a given component
        Args:
            component(str): Formula of the component, i.e "H2" for hydrogen
        """
        if component != None:
            return self.mole_frac(component) * self.gas_stream.moles
        else:
            self.gas_stream.moles

    
    def mass_flow(self, component=None):
        """Returns mass flowrate of a given component
        Args:
            component(str): Formula of the component, i.e "H2" for hydrogen
        """
        if component != None:
            return self.mass_frac(component) * self.gas_stream.mass
        else:
            return self.gas_stream.mass


class GasMixture():

    """ Gas stream mixture that consists of Cantera.Quantity objects i.e GasStream objects.
        Args:
            gas_streams(list): Gas streams to be combined in a mixture.
            >>> GasMixture ( [Gas1, Gas2, ...GasN] )
    """

    def __init__(self, gas_streams):

        self.gas_stream = None
        self.gas = None
        self.create_mixture(gas_streams)
        self.feed_rate = self.gas_stream.moles*1000
    
    def create_mixture(self, gas_streams):
        
        if len(gas_streams) > 2:
            mix = gas_streams[0].gas_stream + gas_streams[1].gas_stream
            for i in range(2, len(gas_streams)):
                mix += gas_streams[i].gas_stream
            self.gas_stream = mix
            self.gas = self.gas_stream._phase

        elif len(gas_streams) == 2:
            mix = gas_streams[0].gas_stream + gas_streams[1].gas_stream
            self.gas_stream = mix
            self.gas = self.gas_stream._phase

        else:
            #print("Gas mixture must consist of at least 2 different gas streams!")
            self.gas_stream = gas_streams[0].gas_stream
            self.gas = gas_streams[0].gas
    

    @property
    def feed_rate(self):
        return self.gas_stream.moles*1000
    
    @feed_rate.setter
    def feed_rate(self, new_flowrate):
        self.gas_stream.moles = new_flowrate/1000
        self._feed_rate = new_flowrate
    
    def enthalpy_stream(self):
        return self.gas_stream.enthalpy   #Total enthalpy in J

    def entropy_stream(self):
        return self.gas_stream.entropy  # Total entropy in J/K

    def gibbs_stream(self):
        return self.gas_stream.gibbs    # Gibbs free energy in J

    def Cp_mole(self):
        return self.gas_stream.cp_mole/1000    # Heat capacity J/mol.K

    def MW (self, component):
        """Returns the Molecular Weight of a given component
        Args:
            component(str): Formula of the component, i.e "H2" for hydrogen
        """
        
        index = self.gas.species_index(component)
        return self.gas.molecular_weights[index]


    def mole_frac(self, component):
        """Returns mole fraction of a given component
        Args:
            component(str): Formula of the component, i.e "H2" for hydrogen
        """

        index = self.gas_stream.species_index(component)
        return self.gas_stream.X[index]
    

    def mass_frac(self, component):
        """Returns mass fraction of a given component
        Args:
            component(str): Formula of the component, i.e "H2" for hydrogen
        """

        index = self.gas_stream.species_index(component)
        return self.gas_stream.Y[index]
    

    def mole_flow(self, component):
        """Returns molar flowrate of a given component
        Args:
            component(str): Formula of the component, i.e "H2" for hydrogen
        """

        return self.mole_frac(component) * self.feed_rate

    
    def mass_flow(self, component):
        """Returns mass flowrate of a given component
        Args:
            component(str): Formula of the component, i.e "H2" for hydrogen
        """

        return self.mass_frac(component) * self.gas_stream.mass
    

class ConversionReactor():

    def __init__(self):
        self.inlets = None
        self.reaction = None
        self.key_reactant = None
        self.conversion = None
        self.heat_inlet = None

        self.thermo_inlet = {}
        self.thermo_outlet = {}
        self.delta = {}
        self.H_rxn = 0
        self.S_rxn = 0
        self.G_rxn = 0

        
    def react (self, inlets, reaction, key_reactant, conversion, heat_inlet= None):
        """ Reaction is carried out in conversional reactor. 
            Args:
                inlets(list): [gas1, gas2 ...] the gas streams to react
                reaction(dict): {"H2O":-1, "H2":1, "O2":0.5} denotes for example H2O --> H2 + 1/2 O2
                key_reactant(str): for example "H2O" for the reaction above
                converison(float): Conversion rate of key reactant
                heat_inlet(float): Heat stream inlet to the reactor in [W]. Default None
        """
        self.inlets = inlets
        self.reaction = reaction
        self.key_reactant = key_reactant
        self.conversion = conversion
        self.heat_inlet = heat_inlet
        
        # Reaction input ex. {"H2O":-1, "H2":1, "O2": 0.5}
        content = GasMixture(self.inlets)           # Combine inlet streams to get reactor mixture
        reactants = []                              # Reactants
        products = []                               # Products
        N_change = {}                               # Molar flowrates of components after the reaction
        N_initial = {}
        coeff = {}
        max_conversion = 0
        
        self.thermo_inlet = {"H":copy(content.enthalpy_stream()), "S":copy(content.entropy_stream()), "G":copy(content.gibbs_stream()) }

        for component, stoic in self.reaction.items():
            N_initial[component] = copy( content.mole_flow(component) )
            coeff[component] = stoic
            if stoic < 0:
                reactants.append(component)
            else:
                products.append(component)
            
        reaction_extend = N_initial[self.key_reactant] * self.conversion

        for component in self.reaction.keys():
            N_change[component] = - reaction_extend * coeff[component]/coeff[self.key_reactant]

        for product in products:
            stream_operation( content, N_change[product], {product: 1.0}  )
        
        for reactant in reactants:
            stream_operation( content, N_change[reactant], {reactant: 1.0}  )

        
        self.thermo_outlet = {"H":copy(content.enthalpy_stream()), "S":copy(content.entropy_stream()), "G":copy(content.gibbs_stream()) }
        self.delta = {"H": self.thermo_outlet["H"] - self.thermo_inlet["H"], "S": self.thermo_outlet["S"] - self.thermo_inlet["S"], "G": self.thermo_outlet["G"] - self.thermo_inlet["G"] }
        self.H_rxn = self.delta["H"]
        self.S_rxn = self.delta["S"]
        self.G_rxn = self.delta["G"]



def h_molar (comp, T):
    gas = ct.Solution("sofc_gas.yaml")
    gas.TPX = T + 273.15, 101325, {comp:1}
    h = gas.enthalpy_mole/1000
    return h

def thermal_conductivity(species, T):
    # Species : {dict}
    gas = ct.Solution("sofc_gas.yaml")
    gas.TPX = T, 101325, species
    return gas.thermal_conductivity



"""
g = GasStream("fuel_channel", T = 600+273.15, feed_rate= 1, X= "H2O:1.0")
r = ConversionReactor() 
r.react([g], {"H2O":-1, "H2":1, "O2":0.5}, "H2O",  0.5)
print(r.H_rxn)
"""
