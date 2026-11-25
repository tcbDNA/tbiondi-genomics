"""Project root and work folder, from this file's own place. The manuscript is an argument and is not in the repository."""
import os
ROOT=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..",".."))+"/"
WORK=ROOT+"results/work/preprint_trace/"
os.makedirs(WORK,exist_ok=True)
