import os
from c2e2 import FileHandler

# Test load model
TOTALV2 = os.path.dirname(os.path.abspath(__file__)) + "/../sample_hyxml/TotalMotionV2.hyxml"

print("=== TEST Begin ===\n")
print("--- Testing FileHandler.load_model(filename) ---")

automata, properties = FileHandler.load_model(TOTALV2)
print(automata)
print(properties)


print("=== End of TEST ===")