import pytest
import altair as alt
import pandas as pd
from rdkit import Chem
from altmol.main import (
    get_rdkit_object,
    mol_to_png_data_url,
    mol_to_svg_data_url,
    encode_molecules,
    mol_plot,
)


def test_get_rdkit_object_molecule():
    smiles = "CCO"  # Ethanol
    result = get_rdkit_object(smiles)
    assert isinstance(result, Chem.Mol)
    assert Chem.MolToSmiles(result) == "CCO"


def test_get_rdkit_object_reaction():
    reaction_smiles = "[O:1]=[C:2]([OH:3])-[H].[H]-[O:4]-[H]>>[O:1]=[C:2]([O:4])-[OH:3]"
    result = get_rdkit_object(reaction_smiles)
    assert isinstance(result, Chem.rdChemReactions.ChemicalReaction)


def test_get_rdkit_object_invalid():
    smiles = "invalid_smiles"
    result = get_rdkit_object(smiles)
    assert result is None


def test_mol_to_png_data_url():
    smiles = "CCO"
    mol = get_rdkit_object(smiles)
    result = mol_to_png_data_url(mol)
    assert result.startswith("data:image/png;base64,")


def test_mol_to_svg_data_url():
    smiles = "CCO"
    mol = get_rdkit_object(smiles)
    result = mol_to_svg_data_url(mol)
    assert result.startswith("data:image/svg+xml;base64,")


def test_encode_molecules():
    data = {"smiles": ["CCO", "C1=CC=CC=C1"]}  # Ethanol and Benzene
    df = pd.DataFrame(data)
    encoded_df = encode_molecules(df)
    assert "image" in encoded_df.columns
    assert "MolecularWeight" in encoded_df.columns
    assert encoded_df["MolecularWeight"][0] > 0  # Check for a positive molecular weight


def test_mol_plot():
    data = {"smiles": ["CCO"], "MolecularWeight": [46.07]}
    df = encode_molecules(pd.DataFrame(data))
    chart = mol_plot(df, x_axis="MolecularWeight", y_axis="smiles", tooltip=["smiles"])
    assert isinstance(chart, alt.Chart)
    assert chart.encoding.x.shorthand == "MolecularWeight"
    assert chart.encoding.y.shorthand == "smiles"
