# First round of model testing
kept for archiving

# how to run
Should be possible to run from this folder without adjustments -  
a best effort for refactoring after files had been moved around without reexecuting.

The paths in the notebook for importing src modules should be like this:

```python
import os
notebookdir = Path.cwd().parents[2]
sys.path.append(str(notebookdir)) # this way we can import src modules even in different subfolders
# and later on:
working_dir = Path.cwd().parents[2]
```

similarly for the py scripts.

Also ensure that you use `from src.db_utils import get_all_data` to `from src.db_utils import get_selected_data` for loading the data.
