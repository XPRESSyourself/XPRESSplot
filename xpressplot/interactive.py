"""
XPRESSplot
A toolkit for navigating and analyzing gene expression datasets
alias: xpressplot

Copyright (C) 2019  Jordan A. Berg
jordan <dot> berg <at> biochem <dot> utah <dot> edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import print_function

"""IMPORT DEPENDENCIES"""
import plotly
import plotly.plotly as py
import plotly.graph_objs as go

"""Create 3D scatterplot for PCA with plotly"""
def interactive_3D(
    df_pca, size, palette,
    plotly_login, save_fig):

    plotly.tools.set_credentials_file(username=str(plotly_login[0]), api_key=str(plotly_login[1]))

    df_pca.columns = ['PCa', 'PCb', 'PCc', 'label']
    unique_labels = df_pca['label'].unique() # Gather unique labels

    pca0 = df_pca.loc[df_pca['label'] == unique_labels[0]]
    pca1 = df_pca.loc[df_pca['label'] == unique_labels[1]]
    pca2 = df_pca.loc[df_pca['label'] == unique_labels[2]]

    x0 = pca0.PCa.values
    y0 = pca0.PCb.values
    z0 = pca0.PCc.values
    trace0 = go.Scatter3d(
        x=x0,
        y=y0,
        z=z0,
        name=unique_labels[0],
        mode='markers',
        marker=dict(
            size=size,
            color=palette[unique_labels[0]],
            opacity=0.8
        )
    )

    x1 = pca1.PCa.values
    y1 = pca1.PCb.values
    z1 = pca1.PCc.values
    trace1 = go.Scatter3d(
        x=x1,
        y=y1,
        z=z1,
        name=unique_labels[1],
        mode='markers',
        marker=dict(
            size=size,
            color=palette[unique_labels[1]],
            opacity=0.8
        )
    )

    x2 = pca2.PCa.values
    y2 = pca2.PCb.values
    z2 = pca2.PCc.values
    trace2 = go.Scatter3d(
        x=x2,
        y=y2,
        z=z2,
        name=unique_labels[2],
        mode='markers',
        marker=dict(
            size=size,
            color=palette[unique_labels[2]],
            opacity=0.8
        )
    )

    data = [trace0, trace1, trace2]
    layout = go.Layout(
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0
        ),
        hovermode='closest'
    )
    fig = go.Figure(data=data, layout=layout)

    if save_fig != None:
        py.offline.plot(fig, filename=str(save_fig))
    else:
        py.iplot(fig, filename='3d_pca')

"""Create interactive scatter plot that displays relevant sample/gene name and coordinates"""
def interactive_scatter(
    data, info, x, y,
    plotly_login, palette,
    highlight='sample', save_fig=None):

    plotly.tools.set_credentials_file(username=plotly_login[0], api_key=plotly_login[1])
    data_c = data.copy()

    if highlight == 'sample':
        # Add labels if not there
        if 'label' not in data_c.index:
            labels = pd.Series(info[1].values,index=info[0]).to_dict()
            data_c.loc['label'] = data_c.columns.map(labels.get)

        data_c = data_c.T
        df_plot = data_c[[x, y,'label']]

        traces = []
        for index, row in df_plot.iterrows():

            _x = [df_plot.loc[index][x]]
            _y = [df_plot.loc[index][y]]

            traces.append(go.Scatter(
                                x = _x,
                                y = _y,
                                name = index,
                                mode = 'markers',
                                marker = dict(
                                    size = 10,
                                    color = palette[df_plot.loc[index]['label']],
                                    line = dict(
                                        width = 2,
                                        color = 'rgb(0, 0, 0)'
                                    )
                                )
                            )
                         )

    elif highlight == 'gene':
        df_plot = data_c[[x, y]]

        traces = []
        for index, row in df_plot.iterrows():

            _x = [df_plot.loc[index][x]]
            _y = [df_plot.loc[index][y]]

            traces.append(go.Scatter(
                                x = _x,
                                y = _y,
                                name = index,
                                mode = 'markers',
                                marker = dict(
                                    size = 2,
                                    line = dict(
                                        width = 2
                                    )
                                )
                            )
                         )

    else:
        print('Error: Invalid input for highlight')
        return

    data = traces
    layout = dict(yaxis = dict(zeroline = False),
                  xaxis = dict(zeroline = False),
                  showlegend=False,
                  hovermode='closest'
                 )

    fig = dict(data=data, layout=layout)
    if save_fig != None:
        py.offline.plot(fig, filename=str(save_fig))
    else:
        py.iplot(fig, filename='2d_scatter')
