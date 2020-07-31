import os
from c2e2 import *

# Load From TotalMotionV2 File
linThermo = os.path.dirname(os.path.abspath(__file__)) + "/linThermo.hyxml"

# Load into Automaton & Properties
automata, properties = FileHandler.load_model(linThermo)

# Get the first property to simulate
prop = properties[0]

TEST_INFO = """
=================== TEST_FILE: {filename} =======================
===== This example will load linThermo.hyxml and simulate   =====
===== the first property of the hyxml file. The length of   =====
===== the trace will be printed in the end.                 =====
=================================================================
""".format(filename=os.path.basename(__file__))
print(TEST_INFO)


# Construct scenario for simulate
sim = Complete_Automata(automata=automata,
                        prop=prop,
                        time_bound=10,
                        unsafe_check=False,
                        name="linThermo-TEST")

# Simulate
trace = sim.simulate()

print("Length of the simulation trace: {}".format(len(trace)))

