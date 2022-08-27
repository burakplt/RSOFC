# 0-D and 1-D Models of a Solid Oxide Fuel/Electrolysis Cell

A Python interface provided by developers of Cantera gives all necessary functionality of the source code. However, the Python interface itself does not include any functional source code, that causes to check the documentation frequently. In order to have an intuitive and maintainable Cantera interface, Python module "gas.py" was developed, which contains all the required functionality related to Cantera.

A gas object defined with Cantera ’Solution’ class does not have a mass
related to it. Instead, a Cantera ’Quantity’ object must also be defined, in order to represent
an amount of a gas. Therefore "gas.py" module contains a Python class ’GasStream’, that
provides the required properties and functions to model a fluid stream with a flow rate.
The ’GasStream’ object includes helper functions, which ensures the consistent and easier code
development. For instance, the molar fraction of a component in a gas is accessed in Cantera
as follows:

**gas_object["component"].X[0]**

where *’gas_object’* is a *’Solution’* object and *’component’* is the chemical formula of the
desired component. The same can be achieved with a helper-function *’mole_frac’* defined in
’GasStream’:

**gas_object.mole_frac("component")**

which is more intuitive and enhances code readability. Similarly, the molecular weight of a
desired component can be achieved:
with Cantera default functions:

**gas_object.molecular_weights[species_index("component")]**

with GasStream "MW" function:

**gas_object.MW("component")**

Moreover, the class **’ConversionReactor’** was added to "gas.py" module, to be able to simulate a
conversional/stoichiometric reactor, as Cantera supports by default only kinetic-based reactions.
This type of reactor is useful with water-splitting reaction, when the SOFC cell used as
electrolysis cell. For a **’ConversionReactor’** the inlet streams, stoichiometric reaction, heat
supply/removal rate and conversion rate of the key species can be specified. As a result the
properties such as change in enthalpy, entropy and Gibbs free energy are calculated.
Although Cantera allows mixing two steams, it is not possible to subtract or add a portion of a
stream by default. With the *’stream_operation’* function in "gas.py" module a desired amount
of a stream with a given composition can be added to or removed from the bulk phase.

A 0-D SOFC model is represented as a Python object **’SOFC_0D’** in **’sofc_0D.py´** module. By
default holds a ’SOFC_0D’ object the fuel and oxygen gas streams, operational parameters,
and material properties as Python class variables. That enables the user to access and alter the
parameters even during the run-time. In addition, multiple instances of a ’SOFC_0D’ object
can be created and executed simultaneously, which could be useful by comparing multiple cell
configurations or making sensitivity analysis. In Table below all class variables of a ’SOFC_0D’
model is listed with their explanations.

![image](https://user-images.githubusercontent.com/22001926/187039490-4cb49c48-6ce4-4363-a799-7c8bc226aa28.png)

# Operational Properties of 0-D Model
’SOFC_0D’ model calculates the steady-state results of a fuel cell operation, if the fuel and
oxygen gases are defined properly. A typical simulation with the 0-D model can be run as
follows:

![image](https://user-images.githubusercontent.com/22001926/187039609-68c1cb9f-b399-4ee1-a3f6-8fb5bbb36934.png)

***0-D SOFC Model schematic representation of main components and their object types.***
![image](https://user-images.githubusercontent.com/22001926/187039671-c45c7103-4050-4c44-86f3-5c32bc4f6861.png)

# Validation of 0-D SOFC Model
![image](https://user-images.githubusercontent.com/22001926/187039710-a40dac18-963c-4b55-904a-a686a5d3c79a.png)
![image](https://user-images.githubusercontent.com/22001926/187039736-30d1c4c4-ae2c-4005-b8c0-151b1a003b42.png)

# 1-D Model 
Convergence problems arise if the initial guesses are badly chosen. No guarantee it will work properly!  
