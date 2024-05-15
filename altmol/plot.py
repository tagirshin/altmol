import altair as alt
import pandas as pd


def chem_plot(
    source: pd.DataFrame,
    x_axis: alt.X,
    y_axis: alt.Y,
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
) -> alt.Chart:
    """
    Generate an Altair plot for molecular data.

    Parameters
    ----------
    source : pd.DataFrame
        Data source for the plot.
    x_axis : alt.X
        Altair X encoding.
    y_axis : alt.Y
        Altair Y encoding.
    tooltip : list, optional
        List of columns to show as tooltips, by default None.
    color : alt.Color, optional
        Color encoding, by default None.
    selector : bool or alt.selection, optional
        Selection for interactivity, by default False.
    interactive : bool, optional
        Enable interactivity, by default True.
    width : int, optional
        Width of the chart, by default 600.
    height : int, optional
        Height of the chart, by default 400.
    title : str, optional
        Title of the chart, by default "".
    chart_type : str, optional
        Type of chart ('scatter' or 'bar'), by default 'scatter'.
    mark_size : int, optional
        Size of the marks, by default 20.
    mark_opacity : float, optional
        Opacity of the marks, by default 1.0.
    mark_fill : bool, optional
        Fill the marks, by default False.
    mark_shape : str, optional
        Shape of the marks, by default "circle".
    mark_color : str, optional
        Color of the marks, by default "#1f77b4".

    Returns
    -------
    alt.Chart
        Altair chart object.

    Notes
    -----
    This function generates a scatter or bar plot with molecular images as tooltips.

    Examples
    --------
    >>> import pandas as pd
    >>> import altair as alt
    >>> from altmol.plot import chem_plot
    >>> data = {
    ...     'x': [1, 2, 3],
    ...     'y': [4, 5, 6],
    ...     'image': [
    ...         'data:image/svg+xml;base64,...',
    ...         'data:image/svg+xml;base64,...',
    ...         'data:image/svg+xml;base64,...'
    ...     ]
    ... }
    >>> df = pd.DataFrame(data)
    >>> x_axis = alt.X('x:Q', title='X Axis')
    >>> y_axis = alt.Y('y:Q', title='Y Axis')
    >>> points = chem_plot(
    ...     df,
    ...     x_axis=x_axis,
    ...     y_axis=y_axis,
    ...     tooltip=['x', 'y'],
    ...     title='Example Plot'
    ... )
    >>> points.show()
    """
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
