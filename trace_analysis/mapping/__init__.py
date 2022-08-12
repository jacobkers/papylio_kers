import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parents[0]))

#from mapping.mapping import Mapping2
#Or rename mapping.py file into main.py file.
# Since the mapping folder is now added to the path, all imports starting with mapping are now recognized as the folder instead of the file.