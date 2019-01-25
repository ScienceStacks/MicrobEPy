"""Initializations for examples."""

import os
import sys

# Add the microbepy code to the path
project_path = os.path.dirname(os.path.abspath(__file__))
project_path = os.path.dirname(project_path)
code_path = os.path.join(project_path, "microbepy")
code_path = project_path
sys.path.insert(0, code_path)
