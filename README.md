# well_stirred_reactor
Simulating  will-stirred reactor

Well-Stirred Reactor Model
This code implements a well-stirred reactor model using Cantera, a popular thermochemical modeling library, to simulate the behavior of a single well-stirred reactor. The well-stirred reactor is modeled using a system of ordinary differential equations (ODEs) that describe the time evolution of temperature and species mass fractions in the reactor. The reactor is assumed to be well-mixed, meaning that the reactants are thoroughly mixed and their concentrations are uniform throughout the reactor.

How to Use
To use this code, you need to have Cantera and scipy libraries installed in your Python environment. Cantera is used for handling chemical kinetics and thermodynamics, while scipy is used for integrating the ODE system.

The WellStirredReactor class is defined in the code, and you can create an instance of this class to simulate the behavior of a well-stirred reactor. The class takes several input parameters, including the reactor volume, residence time, inlet gas temperature, reactor temperature, heat loss, pressure, inlet gas composition, and mechanism file for gas-phase reactions.
