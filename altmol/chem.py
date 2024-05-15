import base64
from io import BytesIO

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdChemReactions import ReactionFromSmarts


def get_rdkit_object(smiles: str):
    """
    Get RDKit object from a SMILES string.

    Parameters
    ----------
    smiles : str
        The SMILES string representing the molecule or reaction.

    Returns
    -------
    object
        RDKit molecule or reaction object.

    Notes
    -----
    This function distinguishes between molecule and reaction SMILES strings.

    Examples
    --------
    >>> mol = get_rdkit_object('CCO')
    >>> type(mol)
    <class 'rdkit.Chem.rdchem.Mol'>
    >>> rxn = get_rdkit_object('CCO>>COC')
    >>> type(rxn)
    <class 'rdkit.Chem.rdChemReactions.ChemicalReaction'>
    """
    if ">" in smiles:
        if ";" in smiles:
            return AllChem.ReactionFromSmarts(smiles)
        else:
            try:
                return ReactionFromSmarts(smiles, useSmiles=True)
            except:
                return AllChem.ReactionFromSmarts(smiles)
    else:
        return Chem.MolFromSmiles(smiles)


def mol_to_png_data_url(mol, width=200, height=200):
    """
    Convert a molecule or reaction object to a PNG data URL.

    Parameters
    ----------
    mol : RDKit object
        RDKit molecule or reaction object.
    width : int, optional
        Width of the image, by default 200.
    height : int, optional
        Height of the image, by default 200.

    Returns
    -------
    str
        PNG image data URL.

    Notes
    -----
    This function generates a PNG image of the given molecule or reaction.

    Examples
    --------
    >>> mol = get_rdkit_object('CCO')
    >>> png_url = mol_to_png_data_url(mol, width=100, height=100)
    >>> print(png_url)
    data:image/png;base64,...
    """
    if mol is None:
        return None
    if isinstance(mol, Chem.Mol):
        img = Draw.MolToImage(mol, size=(width, height))
    else:
        img = Draw.ReactionToImage(mol, size=(width, height))
    with BytesIO() as buffer:
        img.save(buffer, format="PNG")
        data = base64.b64encode(buffer.getvalue()).decode("utf8")
    return f"data:image/png;base64,{data}"


def mol_to_svg_data_url(mol, width=200, height=200):
    """
    Convert a molecule or reaction object to an SVG data URL.

    Parameters
    ----------
    mol : RDKit object
        RDKit molecule or reaction object.
    width : int, optional
        Width of the image, by default 200.
    height : int, optional
        Height of the image, by default 200.

    Returns
    -------
    str
        SVG image data URL.

    Notes
    -----
    This function generates an SVG image of the given molecule or reaction.

    Examples
    --------
    >>> mol = get_rdkit_object('CCO')
    >>> svg_url = mol_to_svg_data_url(mol, width=100, height=100)
    >>> print(svg_url)
    data:image/svg+xml;base64,...
    """
    d2d = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = d2d.drawOptions()
    opts.clearBackground = False
    if isinstance(mol, Chem.Mol):
        d2d.DrawMolecule(mol)
    else:
        d2d.DrawReaction(mol)
    d2d.FinishDrawing()
    img_str = d2d.GetDrawingText()
    with BytesIO() as buffer:
        buffer.write(str.encode(img_str))
        data = base64.b64encode(buffer.getvalue()).decode("utf8")
    return f"data:image/svg+xml;base64,{data}"


def encode_molecules(
    df: pd.DataFrame,
    smiles_col_name="smiles",
    img_format="svg",
    img_width=200,
    img_height=200,
):
    """
    Encode molecules in a DataFrame and generate image URLs.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing molecules.
    smiles_col_name : str, optional
        Column name containing SMILES strings, by default "smiles".
    img_format : str, optional
        Image format ('png' or 'svg'), by default 'svg'.
    img_width : int, optional
        Width of the image, by default 200.
    img_height : int, optional
        Height of the image, by default 200.

    Returns
    -------
    pd.DataFrame
        DataFrame with added image URLs.

    Notes
    -----
    This function adds an 'image' column with data URLs for each molecule.

    Examples
    --------
    >>> df = pd.DataFrame({'smiles': ['CCO', 'CCC']})
    >>> df_encoded = encode_molecules(df, smiles_col_name='smiles', img_format='svg', img_width=100, img_height=100)
    >>> df_encoded.head()
       smiles                                             image
    0    CCO  data:image/svg+xml;base64,...
    1    CCC  data:image/svg+xml;base64,...
    """
    df["Mol"] = df[smiles_col_name].apply(lambda x: get_rdkit_object(x))
    if isinstance(df["Mol"][0], Chem.Mol):
        df["MolecularWeight"] = df["Mol"].apply(lambda x: Descriptors.MolWt(x))
    if img_format == "png":
        df["image"] = df["Mol"].apply(mol_to_png_data_url, args=(img_width, img_height))
    elif img_format == "svg":
        df["image"] = df["Mol"].apply(mol_to_svg_data_url)
    return df.drop(columns=["Mol"])
