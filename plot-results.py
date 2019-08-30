from plotly.offline import plot
from plotly.io import write_image
import plotly.graph_objs as go

import pandas as pd

data = pd.read_excel("results/results_other-methods.xlsx")

print(data.values[1, :])

trace1 = go.Bar(
    x=data.columns.values[[1, 2, 3, 4, 5]],
    y=data.values[0, [1, 2, 3, 4, 5]],
    name='pcaReduce'
)
trace2 = go.Bar(
    x=data.columns.values[[1, 2, 3, 4, 5]],
    y=data.values[1, [1, 2, 3, 4, 5]],
    name='SEURAT'
)

trace3 = go.Bar(
    x=data.columns.values[[1, 2, 3, 4, 5]],
    y=data.values[2, [1, 2, 3, 4, 5]],
    name='SC3'
)

trace4 = go.Bar(
    x=data.columns.values[1:],
    y=data.values[3, 1:],
    name='SIMLR'
)


trace5 = go.Bar(
    x=data.columns.values[1:],
    y=data.values[4, 1:],
    name='SINCERA'
)

trace6 = go.Bar(
    x=data.columns.values[1:],
    y=data.values[5, 1:],
    name='TSCAN'
)

trace7 = go.Bar(
    x=data.columns.values[[1, 2, 3, 4, 5]],
    y=data.values[6, [1, 2, 3, 4, 5]],
    name='FeatClust'
)

data = [trace7, trace1, trace2, trace3, trace4, trace5, trace6]
layout = go.Layout(
    barmode='group',
    xaxis = dict(title = 'Datasets'),
    yaxis = dict(title = 'Adjusted Rand Index'))


fig = go.Figure(data=data, layout=layout)
plot(fig, filename='ari-barplot.html')

write_image(fig, 'ari-barplot.pdf')