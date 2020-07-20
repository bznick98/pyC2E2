# pyC2E2

# Current Functionalities
- [x] Be able to load model from .hyxml file.
- [ ] Using C2E2 C++ Backend to do Simualtion/Verification

# Install
```zsh
# Clone this repo.
git clone https://github.com/bznick98/pyC2E2.git
# Go to the project directory
cd pyC2E2
# install into global environment
pip install .
```

# Use Case
```python
from pyC2E2 import FileHandler

# Using FileHandler to load a hyxml file
automata, properties = FileHandler.load_model("TotalMotionV2.hyxml")

# All the automata will be loaded into 'automata'
# All the properties will be loaded in to 'properties'
# You can display info of an automaton by directly printing automaton's name
print(automata)
print(properties)
```