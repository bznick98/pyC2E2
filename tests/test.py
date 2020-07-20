from FileHandler import *

# Test load model
TOTALV2 = "../TotalMotionV2.hyxml"

print("=== TEST Begin ===\n")
print("--- Testing FileHandler.load_model(filename) ---")

models = FileHandler.load_model(TOTALV2)
print(models)


print("=== End of TEST ===")