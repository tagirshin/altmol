# altmol

`altmol` is a Python package designed to enrich Altair visualizations with interactive 2D molecular structures, 
inspired by the capabilities of [molplotly](https://github.com/wjm41/molplotly). 
The current version uses RDKit for molecular rendering.

## Installation

`altmol` can be easily installed using `pip` or `poetry`, accommodating both traditional and modern Python workflows.

### Using pip

```sh
pip install altmol
```

### Using Poetry

If you're managing your project with Poetry, add `altmol` to your project using:

```sh
poetry add altmol
```

Alternatively, you can manually add it to your `pyproject.toml` and then install it using:

```sh
poetry install
```

## Quick Start

To get started with `altmol`, ensure you have Altair and RDKit installed in your environment. 
Here's a simple example to illustrate how to create an interactive scatter plot with molecule visualizations:

```python
import altair as alt
import pandas as pd

from altmol import encode_molecules, mol_scatter


# URL of the ESOL dataset
url_esol = 'https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/delaney-processed.csv'

# Load your dataset
df = pd.read_csv(url_esol)

# Prepare your data (ensure you have a 'smiles' column for molecular structures)
df = encode_molecules(df, smiles_col_name='smiles', img_format='svg')

x_col_name = "measured log solubility in mols per litre"
y_col_name = "ESOL predicted log solubility in mols per litre"

# Generate and display the interactive plot
chart = mol_scatter(
    df,
    x_axis=alt.X(f"{x_col_name}:Q", title="True Log solubility"),
    y_axis=alt.Y(f"{y_col_name}:Q", title="Pred Log solubility"),
    tooltip=[x_col_name, y_col_name],
    title="ESOL Regression"
)
chart.display()
```

## Documentation

For detailed documentation, including a full list of features, installation instructions, 
and advanced usage examples, please visit altmol documentation.

## Contributing

Contributions to `altmol` are welcome! Whether it's bug reports, feature requests, 
or contributions to the code, we value your input. Please refer to our Contributing Guidelines 
for more information on how you can contribute.

## License

`altmol` is released under the [MIT License](LICENSE.txt). See the LICENSE file for more details.

## Acknowledgements

`altmol` is built upon the powerful capabilities of [RDKit](https://www.rdkit.org/) 
and [Altair](https://altair-viz.github.io/), and we are grateful to the developers and contributors of these projects.

 