"""Database utility functions for retrieving and updating data using SQLAlchemy
20205-10-28 Alexander Minidis
Update: 10/29: descriptor read/write, get all data with descriptors
Update: 29/04/26: new get_selected_data function to exclude specified columns.
Easiest to include the point-estimate-SVR model with the interval-based model
without having to do major refactoring.
"""

from typing import Optional, List

import pandas as pd

DEFAULT_EXCLUDE_COLUMNS = [
    "T_half_class_days",
    "T_half_lower_bound",
    "T_half_upper_bound",
    "T_half_log10_lower",
    "T_half_log10_upper",
]


def _get_table_map():
    """Return a mapping of table identifiers to their corresponding SQLAlchemy ORM classes.
    This module is only loaded once, and subsequent imports are just lookups, so there’s minimal overhead.

    """
    from src.db_schema import AirData, SedimentData, SoilData, WaterData

    return {
        "air": AirData,
        "soil": SoilData,
        "water": WaterData,
        "sediment": SedimentData,
    }


def get_basic_data(table_id: str, Session) -> pd.DataFrame:
    """Retrieve basic data (reference, T_half_days, Canonical_smiles) from the specified table.
    table_id should be one of: "air", "soil", "water", "sediment"
    """
    session = Session()
    _table_map = _get_table_map()

    if table_id not in _table_map:
        raise ValueError("Invalid table_id. Must be one of: 'air', 'soil', 'water', 'sediment'.")
    table = _table_map[table_id]

    df = pd.DataFrame(
        session.query(table.reference, table.T_half_days, table.Canonical_smiles).all(),
        columns=["reference", "T_half_days", "Canonical_smiles"],
    )
    session.close()
    return df


def get_all_data(table_id: str, Session) -> pd.DataFrame:
    """Retrieve all data from the specified table as a pandas DataFrame.
    table_id should be one of: "air", "soil", "water", "sediment"

    Retrieves all columns, including descriptors and T_half bounds, which may be useful for the interval-based SVR model or future models that utilize the additional information.
    """
    session = Session()
    _table_map = _get_table_map()

    if table_id not in _table_map:
        raise ValueError("Invalid table_id. Must be one of: 'air', 'soil', 'water', 'sediment'.")
    table = _table_map[table_id]

    df = pd.read_sql(session.query(table).statement, session.bind)
    session.close()
    return df


def get_selected_data(table_id: str, Session, exclude_columns: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Retrieve all data from the specified table except the columns in exclude_columns.

    table_id should be one of: "air", "soil", "water", "sediment"

    Currently set to exclude the additional T_half upper/lower bound columns for use in the legacy point-based SVR model which only uses  T_half_days column as  target variable.

    """

    if exclude_columns is None:
        exclude_columns = DEFAULT_EXCLUDE_COLUMNS

    session = Session()
    _table_map = _get_table_map()

    if table_id not in _table_map:
        raise ValueError("Invalid table_id. Must be one of: 'air', 'soil', 'water', 'sediment'.")
    table = _table_map[table_id]

    # Get all column names from the ORM class, excluding those in exclude_columns
    columns = [col for col in table.__table__.columns if col.name not in exclude_columns]

    # Query only the selected columns
    df = pd.DataFrame(
        session.query(*columns).all(),
        columns=[col.name for col in columns],
    )
    session.close()
    return df


def get_row_from_smiles(table_id: str, canonical_smiles: str, Session) -> Optional[pd.Series]:
    """
    Retrieve a single row from the specified table by Canonical_smiles.
    Returns a pandas Series of only the descriptor columns or None if not found.
    """
    session = Session()
    _table_map = _get_table_map()

    if table_id not in _table_map:
        raise ValueError("Invalid table_id. Must be one of: 'air', 'soil', 'water', 'sediment'.")
    table = _table_map[table_id]

    row = session.query(table).filter(table.Canonical_smiles == canonical_smiles).first()
    session.close()
    if row:
        return pd.Series(row.__dict__).drop("_sa_instance_state", errors="ignore")
    else:
        return None


def write_descriptors(table_id: str, reference: str, descriptors: dict, Session) -> bool:
    """
    Update a single row's descriptor columns by reference using _table_map.
    Returns True if update was successful, False otherwise.
    """
    session = Session()
    _table_map = _get_table_map()

    if table_id not in _table_map:
        raise ValueError("Invalid table_id. Must be one of: 'air', 'soil', 'water', 'sediment'.")
    table = _table_map[table_id]

    row = session.query(table).filter(table.reference == reference).first()
    if not row:
        session.close()
        return False

    for key, value in descriptors.items():
        if hasattr(row, key):
            setattr(row, key, value)
    session.commit()
    session.close()
    return True
