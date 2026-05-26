# First round of model testing
kept for archiving

# how to run
Will not run from this folder without adjustments, e.g. to file paths such as src/, data loading, etc.

Either move one level up or add before the `src`related imports these lines:

```python
import os
notebookdir = Path.cwd().parent.parent
sys.path.append(str(notebookdir)) # this way we can import src modules even in different subfolders
# and later on:
working_dir = Path.cwd().parent.parent
```

also ensure that you use `from src.db_utils import get_all_data` to `from src.db_utils import get_selected_data` for loading the data.
