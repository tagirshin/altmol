import base64
from io import BytesIO

import altair as alt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdChemReactions import ReactionFromSmarts


def get_rdkit_object(smiles: str):
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


def encode_molecules(
    df: pd.DataFrame,
    smiles_col_name="smiles",
    img_format="svg",
    img_width=200,
    img_height=200,
):
    df["Mol"] = df[smiles_col_name].apply(lambda x: get_rdkit_object(x))
    if isinstance(df["Mol"][0], Chem.Mol):
        df["MolecularWeight"] = df["Mol"].apply(lambda x: Descriptors.MolWt(x))
    if img_format == "png":
        df["image"] = df["Mol"].apply(mol_to_png_data_url, args=(img_width, img_height))
    elif img_format == "svg":
        df["image"] = df["Mol"].apply(mol_to_svg_data_url)
    return df.drop(columns=["Mol"])


def mol_plot(
    source,
    x_axis,
    y_axis,
    tooltip: list = None,
    color=None,
    selector=False,
    interactive=True,
    width=600,
    height=400,
    title="",
    chart_type="scatter",
    mark_size=20,
    mark_opacity=1.0,
    mark_fill=False,
    mark_shape="circle",
    mark_color="#1f77b4",
):
    if tooltip is None:
        tooltip = []
    if isinstance(selector, bool) and selector:
        selector = alt.selection_point()

    if selector and color:
        color = alt.condition(selector, color.legend(None), alt.value("lightgray"))
    elif selector:
        color = alt.condition(selector, alt.value(mark_color), alt.value("lightgray"))

    if color is None:
        color = alt.Color()

    tooltip.append("image")
    # Altair chart with tooltips displaying PNG images, ensuring the column is named 'image'
    if chart_type == "scatter":
        chart = alt.Chart(source).mark_point(
            size=mark_size,
            opacity=mark_opacity,
            filled=mark_fill,
            shape=mark_shape,
            color=mark_color,
        )
    elif chart_type == "bar":
        chart = alt.Chart(source).mark_bar(
            size=mark_size,
            opacity=mark_opacity,
            filled=mark_fill,
            shape=mark_shape,
            color=mark_color,
        )
    else:
        raise NotImplementedError(f"Unsupported chart type: {chart_type}")

    chart = chart.encode(
        x=x_axis,
        y=y_axis,
        color=color,
        tooltip=tooltip,  # Referencing 'image' column for tooltips
    ).properties(
        width=width,
        height=height,
        title=title,
    )
    if selector:
        chart = chart.add_params(selector)
    if interactive:
        chart = chart.interactive()
    return chart
