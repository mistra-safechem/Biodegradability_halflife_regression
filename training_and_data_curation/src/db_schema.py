"""Database schema definitions using SQLAlchemy
20205-10-28 Alexander Minidis
Updates: 1/11 schema updated, descriptors added
Updates: 04/21 interval mapping columns added, descriptors/fingerprints made non-nullable
"""

from typing import Optional

from sqlalchemy import Float, Integer, String
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column
from src.rdkit_tools import DESCRIPTOR_NAMES


class Base(DeclarativeBase):
    pass


class AirData(Base):
    __tablename__ = "air_data"
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    reference: Mapped[str] = mapped_column(String)
    T_half_days: Mapped[float] = mapped_column(Float)
    # interval mapping columns (NULL for unmapped/missing rows)
    T_half_class_days: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_lower_bound: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_upper_bound: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_log10_lower: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_log10_upper: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    Canonical_smiles: Mapped[str] = mapped_column(String)
    # rdkit descriptors
    for name in DESCRIPTOR_NAMES:
        vars()[name] = mapped_column(Float)
    # MACCS fingerprints
    for i in range(1, 167):
        vars()[f"MACCS_{i:03d}"] = mapped_column(Float)


class SoilData(Base):
    __tablename__ = "soil_data"
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    reference: Mapped[str] = mapped_column(String)
    T_half_days: Mapped[float] = mapped_column(Float)
    # interval mapping columns (NULL for unmapped/missing rows)
    T_half_class_days: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_lower_bound: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_upper_bound: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_log10_lower: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_log10_upper: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    Canonical_smiles: Mapped[str] = mapped_column(String)
    # rdkit descriptors
    for name in DESCRIPTOR_NAMES:
        vars()[name] = mapped_column(Float)
    # MACCS fingerprints
    for i in range(1, 167):
        vars()[f"MACCS_{i:03d}"] = mapped_column(Float)


class WaterData(Base):
    __tablename__ = "water_data"
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    reference: Mapped[str] = mapped_column(String)
    T_half_days: Mapped[float] = mapped_column(Float)
    # interval mapping columns (NULL for unmapped/missing rows)
    T_half_class_days: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_lower_bound: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_upper_bound: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_log10_lower: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_log10_upper: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    Canonical_smiles: Mapped[str] = mapped_column(String)
    # rdkit descriptors
    for name in DESCRIPTOR_NAMES:
        vars()[name] = mapped_column(Float)
    # MACCS fingerprints
    for i in range(1, 167):
        vars()[f"MACCS_{i:03d}"] = mapped_column(Float)


class SedimentData(Base):
    __tablename__ = "sediment_data"
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    reference: Mapped[str] = mapped_column(String)
    T_half_days: Mapped[float] = mapped_column(Float)
    # interval mapping columns (NULL for unmapped/missing rows)
    T_half_class_days: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_lower_bound: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_upper_bound: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_log10_lower: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    T_half_log10_upper: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    Canonical_smiles: Mapped[str] = mapped_column(String)
    # rdkit descriptors
    for name in DESCRIPTOR_NAMES:
        vars()[name] = mapped_column(Float)
    # MACCS fingerprints
    for i in range(1, 167):
        vars()[f"MACCS_{i:03d}"] = mapped_column(Float)
