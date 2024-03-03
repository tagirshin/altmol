import base64
from io import BytesIO

import altair as alt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdChemReactions import ReactionFromSmarts


def get_rdkit_object(smiles: str):
    if ">" in smiles:
        return ReactionFromSmarts(smiles, useSmiles=True)
    else:
        return Chem.MolFromSmiles(smiles)


def mol_to_png_data_url(mol, width=200, height=200):
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


def encode_molecules(df: pd.DataFrame, smiles_col_name="smiles", img_format="svg"):
    df["Mol"] = df[smiles_col_name].apply(lambda x: get_rdkit_object(x))
    df["MolecularWeight"] = df["Mol"].apply(lambda x: Descriptors.MolWt(x))
    if img_format == "png":
        df["image"] = df["Mol"].apply(mol_to_png_data_url)
    elif img_format == "svg":
        df["image"] = df["Mol"].apply(mol_to_svg_data_url)
    return df.drop(columns=["Mol"])


def mol_scatter(
    source,
    x_axis,
    y_axis,
    tooltip: list = None,
    width=600,
    height=400,
    title="",
    mark_size=20,
    mark_opacity=1.0,
    mark_fill=False,
    mark_shape="circle",
):
    if tooltip is None:
        tooltip = []
    tooltip.append("image")
    # Altair chart with tooltips displaying PNG images, ensuring the column is named 'image'
    chart = (
        alt.Chart(source)
        .mark_point(
            size=mark_size, opacity=mark_opacity, filled=mark_fill, shape=mark_shape
        )
        .encode(
            x=x_axis,
            y=y_axis,
            tooltip=tooltip,  # Referencing 'image' column for tooltips
        )
        .properties(
            width=width,
            height=height,
            title=title,
        )
        .interactive()
    )

    return chart
