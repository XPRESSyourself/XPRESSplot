"""
XPRESStools
A toolkit for navigating and analyzing gene expression datasets
alias: xpresstools

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

"""
IMPORT DEPENDENCIES
"""
import plotly
import plotly.plotly as py
import plotly.graph_objs as go

"""
DESCRIPTION: Create interactive scatter plot that displays relevant sample/gene name and coordinates

VARIABLES:
data= MICARtools formatted expression matrix
"""
def interactive_scatter(data, info, x, y, plotly_login, file_name, gene_list=None, palette=None):

    plotly.tools.set_credentials_file(username=plotly_login[0], api_key=plotly_login[1])
    data_c = data.copy()

    if highlight == 'sample':
        #Add labels if not there
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
    py.offline.plot(fig, filename=str(file_name))
