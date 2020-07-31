import os

from c2e2 import *

# Loading a DryVR Black Box File
thermo_filename = os.path.dirname(os.path.abspath(__file__)) + "/Thermostats_ODE.py"

# FileHandler
automaton = FileHandler.load_model(thermo_filename)

TEST_INFO = """
=================== TEST_FILE: {filename} =======================
===== This example will load Thermostats.py and simulate    =====
===== by passing automaton constructed by the blackbox and  =====
===== some initial conditions, the length of the trace will =====
=====  be printed in the end.                               =====
=================================================================
""".format(filename=os.path.basename(__file__))
print(TEST_INFO)

# Set Up Simulation Scenario
complete = Complete_Automata(automata=automaton,
                             initial_state=[65],
                             initial_mode="On",
                             time_bound=10,
                             name="DryVR_Thermo-TEST")

# simulation trace
trace = complete.simulate()

# Display trace
print("Length of the simulation trace: {}".format(len(trace)))