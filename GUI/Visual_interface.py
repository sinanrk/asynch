# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 13:56:29 2022

@author: srasiyakoya
"""

import pandas as pd
import plotly.express as px  # (version 4.7.0 or higher)
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output
import dash_bootstrap_components as dbc

#import dash_core_components as dcc
#import dash_html_components as html

import os

#%% FUNCTIONS
def read_dat2df(qsim, tfreq, sdate, edate): #sdate and edate are dt.datetime
    link2usgs = pd.DataFrame({'Link':['17225', '29742', '34424', '51645', '50301', '72675', '65483' ], #'52025'- 'NA', '65483'-'06799445', '62576' check the numbers again
                          'STAID':['06800000','06799500', '06799350', '06799315', '06799000', '06799100', '06797500']})
    # Up todowstream -> 'STAID':['06797500', '06799100', '06799000', '06799315', '06799445', '06799350', '06799500', '06800000']
    nlinks = int(qsim[0][0])
    nstates = int(qsim[0][1])
    pind = 2
    date_list = pd.date_range(sdate,edate,freq=tfreq)
    op = pd.DataFrame(index=date_list)
    for n in np.arange(0,nlinks) :
        linkid = qsim[0][pind].split(" ")[0]
        npoints = int(qsim[0][pind].split(" ")[1])
        curr_mat = qsim[0][pind+1:pind+npoints+1]
        df_now = pd.DataFrame()
        for s in np.arange(0,nstates):
            ts_now = np.array(curr_mat.str.split(' ', expand=True)[s].astype('float32'))          
            op[linkid+"_op_"+str(s)] = ts_now
        pind = pind+npoints+1
    return op

def usgs_txt2df(usgs_sno):
    usgs_q = pd.read_table('E:\\DATA\\USGS\\Elkhorn\\'+str(int(usgs_sno)).zfill(8)+'.txt', sep="\t", skiprows=37, header=None)
    print(usgs_q)
    usgs_q = usgs_q[[2,6]] #change columns according to input txt
    usgs_q.columns = ['date', 'Q'] # {'Q', 'date'} #
    usgs_q['date'] = pd.to_datetime(usgs_q['date'])
    usgs_q = usgs_q.loc[((usgs_q['date'] >= '2018-01-01') & (usgs_q['date'] < '2020-01-01'))]
    usgs_q['Q'] = usgs_q['Q']*0.028316847 # cfs to cumec
    #usgs_q = usgs_q.fillna(method='bfill')
    usgs_q = usgs_q.dropna()
    usgs_q = usgs_q.set_index('date')
    usgs_q = usgs_q.astype({'Q': 'float64'})
    usgs_q = usgs_q.resample('H').mean()
    return usgs_q

#%% LOAD THE DATA

fname = r'C:\Users\srasiyakoya\OneDrive - University of Nebraska-Lincoln\MATC_MyWork\Elkhorn\Elkhorn_GUI\DATA\dat_files\m608new_2018_2019.csv'
bfnm = os.path.basename(fname)
qsim = pd.read_csv(fname)
qsim['Time'] = pd.to_datetime(qsim['Unnamed: 0'])
qsim = qsim.set_index('Time')
Qdf = qsim.loc[:,qsim.columns.str.endswith('op_0')] # Subsetting Discharge only
#Qdf.plot(subplots =True) # Data check

DF_xy = pd.read_csv(r'.\Elkhorn_GUI\Components\Link_Centroids.csv')
#DF_xy= DF_xy.iloc[:,0:4] # For CLICK ON USGS STATIONS 

mapbox_access_token = open(r'.\Elkhorn_GUI\mapbox_token.txt', 'r')
mapbox_access_token = mapbox_access_token.read()

#%% MAIN

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# CREATE MAP GRAPH OBJECT
df_sub = DF_xy
    
# Create figure
locations=[go.Scattermapbox(
                lon = df_sub['xCentroid'],
                lat = df_sub['yCentroid'],
                mode='markers',
                marker={'color' :'blue', 'size':5},
                unselected={'marker' : {'opacity':1}},
                selected={'marker' : {'opacity':0.5, 'size':25}},
                hoverinfo='text',
                #hovertext=df_sub['STANAME'],
                customdata=df_sub['LINKNO']
)]

# LAYOUT
app.layout = dbc.Container([
    dbc.Row([
       dbc.Col([
           html.H1("Elkhorn Discharge Level", style={'text-align': 'center'})    
       ], width=8)
    ], justify="center"),

    dbc.Row([
        dbc.Col([
            html.Label('Station Details'),
            html.Br(),
           html.Div(id='output_container', children=[])
        ], width=3, align="center"),
        dbc.Col([
            dcc.Graph(id='Elk_graph', figure={'data': locations,
                                                   'layout': go.Layout( uirevision= 'foo', #preserves state of figure/map after callback activated
                                                                        clickmode= 'event+select',
                                                                        hovermode= 'closest',
                                                                        hoverdistance=2,
                                                                        title=dict(text="Click anywhere on the streams",font=dict(size=20, color='green')),
                                                                        mapbox=dict(accesstoken=mapbox_access_token,
                                                                                    #bearing=25,
                                                                                    #style='light',
                                                                                    center=dict(lat=41.9952778,lon=-97.0580556),
                                                                                    #pitch=40,
                                                                                    zoom=6),
                                                                     )
                                                   },
                        config={'displayModeBar': False, 'scrollZoom': True})#,
                        #style={'background':'#00FC87','padding-bottom':'1px','padding-left':'1px','height':'70vh'}),
        ], width=9)
    ], align="center"),

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='sim_hydrograph', figure={})
                ], width=12)
            ]),
])

# CALLBACKS

@app.callback(
    [Output(component_id='output_container', component_property='children'),
     Output(component_id='sim_hydrograph', component_property='figure')],
    [Input('Elk_graph', 'clickData')])

def display_click_data(clickData):
    
    if clickData is None:
        print('----------------ClickData NONE--------------------')
        print(clickData)
        #print(clickData['points'][0]['customdata'])
        option_slctd_link = '17225'
        
    else:   
        print('----------------ClickData THERE IS SOMETHING--------------------')
        print(clickData)
        option_slctd_link=clickData['points'][0]['customdata']
        option_slctd_link=str(int(option_slctd_link))
        print('HLM Link -------------'+option_slctd_link)
        
    # usgs_sno = link2usgs.loc[link2usgs['Link'] == lnk]['STAID']  
    ch_link = option_slctd_link
    ch_X = str(DF_xy.loc[DF_xy['LINKNO']==int(option_slctd_link)].xCentroid.values[0])
    
    ch_Y = str(DF_xy.loc[DF_xy['LINKNO']==int(option_slctd_link)].yCentroid.values[0])
    ch_strOrd = str(DF_xy.loc[DF_xy['LINKNO']==int(option_slctd_link)].strmOrder.values[0])
    
    container = html.Ul([
                            html.Li("Selected HLM Link: {}".format(ch_link), 
                                    className='circle', 
                                    style={'background': 'white','color':'black', 'list-style':'none','text-indent': '10px'}),
                            html.Li("Selected Latitude: {}".format(ch_Y), 
                                    className='circle', 
                                    style={'background': 'white','color':'black', 'list-style':'none','text-indent': '10px'}),
                            html.Li("Selected Longitude: {}".format(ch_X), 
                                    className='circle', 
                                    style={'background': 'white','color':'black', 'list-style':'none','text-indent': '10px'}),
                            html.Li("Order of Selected Stream: {}".format(ch_strOrd), 
                                    className='circle', 
                                    style={'background': 'white','color':'black', 'list-style':'none','text-indent': '10px'})
                            
                            
                            # html.Li("Electronics", className='circle', style={'background': '#0000ff','color':'black',
                            #     'list-style':'none','text-indent': '17px','white-space':'nowrap'}),
                            # html.Li("Hazardous_waste", className='circle', style={'background': '#FF0000','color':'black',
                            #     'list-style':'none','text-indent': '17px'}),
                            # html.Li("Plastic_bags", className='circle', style={'background': '#00ff00','color':'black',
                            #     'list-style':'none','text-indent': '17px'}),
                            # html.Li("Recycling_bins", className='circle',  style={'background': '#824100','color':'black',
                            #     'list-style':'none','text-indent': '17px'}),
                            
                        ]#, style={'border-bottom': 'solid 3px', 'border-color':'#00FC87','padding-top': '6px'}
                        )
    
    dff = Qdf[option_slctd_link+'_op_0'] # option_slctd_link = "17225"

    fig = px.line(x=dff.index,
              y=dff.values,
              labels = dict(x="Time", y="Discharge (m3/s)"),
              title = 'Discharge at HLM '+option_slctd_link)
    return container, fig


        
if __name__ == '__main__':
    app.run_server(debug=True, use_reloader=False)



#%% CLICK ON USGS STATIONS  

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# LOAD THE DATA
fname = r'C:\Users\srasiyakoya\OneDrive - University of Nebraska-Lincoln\MATC_MyWork\Elkhorn\HLM608new\final\SIM\m608new_2018_2019.csv'
bfnm = os.path.basename(fname)
qsim = pd.read_csv(fname)
qsim['Time'] = pd.to_datetime(qsim['Unnamed: 0'])
qsim = qsim.set_index('Time')
Qdf = qsim.loc[:,qsim.columns.str.endswith('op_0')] # Subsetting Discharge only
Qdf.plot(subplots =True) # Data check

DF_xy = pd.read_csv(r'.\Elkhorn_GUI\Components\USGS_XY.csv')
DF_xy= DF_xy.iloc[:,0:4]

mapbox_access_token = open(r'.\Elkhorn_GUI\mapbox_token.txt', 'r')
mapbox_access_token = mapbox_access_token.read()


# CREATE MAP GRAPH OBJECT
df_sub = DF_xy
    
# Create figure
locations=[go.Scattermapbox(
                lon = df_sub['X'],
                lat = df_sub['Y'],
                mode='markers',
                marker={'color' :'red', 'size':5},
                unselected={'marker' : {'opacity':1}},
                selected={'marker' : {'opacity':0.5, 'size':25}},
                hoverinfo='text',
                hovertext=df_sub['STANAME'],
                customdata=df_sub['STAID']
)]

# LAYOUT
app.layout = dbc.Container([
    dbc.Row([
       dbc.Col([
           html.H1("Elkhorn Discharge Level", style={'text-align': 'center'})    
       ], width=8)
    ], justify="center"),

    dbc.Row([
        dbc.Col([
            html.Label('USGS Station'),
           html.Div(id='output_container', children=[])
        ], width=3),
        dbc.Col([
            dcc.Graph(id='Elk_graph', figure={'data': locations,
                                                   'layout': go.Layout( uirevision= 'foo', #preserves state of figure/map after callback activated
                                                                        clickmode= 'event+select',
                                                                        hovermode= 'closest',
                                                                        hoverdistance=2,
                                                                        title=dict(text="USGS Stations at Elkhorn",font=dict(size=20, color='green')),
                                                                        mapbox=dict(accesstoken=mapbox_access_token,
                                                                                    #bearing=25,
                                                                                    #style='light',
                                                                                    center=dict(lat=41.9952778,lon=-97.0580556),
                                                                                    #pitch=40,
                                                                                    zoom=6),
                                                                     )
                                                   },
                        config={'displayModeBar': False, 'scrollZoom': True})#,
                        #style={'background':'#00FC87','padding-bottom':'1px','padding-left':'1px','height':'70vh'}),
        ], width=9)
    ]),

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='sim_hydrograph', figure={})
                ], width=12)
            ]),
])

# CALLBACKS

@app.callback(
    [Output(component_id='output_container', component_property='children'),
     Output(component_id='sim_hydrograph', component_property='figure')],
    [Input('Elk_graph', 'clickData')])

def display_click_data(clickData):
    
    link2usgs = pd.DataFrame({'Link':['17225', '29742', '34424', '51645', '50301', '72675', '65483' ], #'52025'- 'NA', '65483'-'06799445', '62576' check the numbers again
                      'STAID':['06800000','06799500', '06799350', '06799315', '06799000', '06799100', '06797500']})
    
    if clickData is None:
        print('----------------ClickData NONE--------------------')
        print(clickData)
        #print(clickData['points'][0]['customdata'])
        option_slctd_link = '17225'
        
    else:   
        print('----------------ClickData THERE IS SOMETHING--------------------')
        print(clickData)
        option_slctd_USGS=clickData['points'][0]['customdata']
        option_slctd_USGS=str(int(option_slctd_USGS)).zfill(8)
        print('USGS---------'+option_slctd_USGS)
        option_slctd_link = link2usgs.loc[link2usgs['STAID']==option_slctd_USGS].Link.values[0]
        print('HLM Link -------------'+option_slctd_link)
        
    # usgs_sno = link2usgs.loc[link2usgs['Link'] == lnk]['STAID']  
    container = "The link chosen by user is: {}".format(option_slctd_link)

    dff = Qdf.copy()
    dff = Qdf[option_slctd_link+'_op_0'] # option_slctd_link = "17225"

    fig = px.line(x=dff.index,
              y=dff.values,
              labels = dict(x="Time", y="Discharge (m3/s)"),
              title = 'Discharge at HLM '+option_slctd_link)
    return container, fig
        
if __name__ == '__main__':
    app.run_server(debug=True, use_reloader=False)

#%% DRPDOWN USGS STATIONS

app = Dash(__name__)

# LOAD THE DATA
fname = r'C:\Users\srasiyakoya\OneDrive - University of Nebraska-Lincoln\MATC_MyWork\Elkhorn\HLM608new\final\SIM\m608new_2018_2019.csv'
bfnm = os.path.basename(fname)
qsim = pd.read_csv(fname)
qsim['Time'] = pd.to_datetime(qsim['Unnamed: 0'])
qsim = qsim.set_index('Time')
Qdf = qsim.loc[:,qsim.columns.str.endswith('op_0')] # Subsetting Discharge only
Qdf.plot(subplots =True) # Data check

DF_xy = pd.read_csv(r'.\Elkhorn_GUI\Components\USGS_XY.csv')
DF_xy= DF_xy.iloc[:,0:4]

mapbox_access_token = open(r'.\Elkhorn_GUI\mapbox_token.txt', 'r')
mapbox_access_token = mapbox_access_token.read()

# APP LAYOUT
app.layout = html.Div([

    html.H1("Elkhorn Discharge Level", style={'text-align': 'center'}),

    html.Div(children = [dcc.Dropdown(id="slct_link",
                 options=[
                     {"label": "06800000", "value": "17225"},
                     {"label": "06799500", "value": "29742"},
                     {"label": "06799350", "value": "34424"},
                     {"label": "06799315", "value": "51645"},
                     {"label": "06799000", "value": "50301"},
                     {"label": "06799100", "value": "72675"},
                     {"label": "06797500", "value": "65483"}],
                 multi=False,
                 value="17225",
                 style={'width': "40%"}
                 )],
             className='left_menu'),

    html.Div(id='output_container', children=[]),
    html.Br(),
    
    html.Div([dcc.Graph(id='Elk_graph', 
                        config={'displayModeBar': False, 'scrollZoom': True},
                        style={'background':'#00FC87','padding-bottom':'2px','padding-left':'2px','height':'70vh'})], 
             className='right_content'),
    
    dcc.Graph(id='sim_hydrograph', figure={})

])

# Connect the Plotly graphs with Dash Components
@app.callback(
    [Output(component_id='output_container', component_property='children'),
     Output(component_id='sim_hydrograph', component_property='figure')],
    [Input(component_id='slct_link', component_property='value')]
)
def update_graph(option_slctd_USGS):
    print(option_slctd_USGS)
    print(type(option_slctd_USGS))
    
    

    dff = Qdf.copy()
    dff = Qdf[option_slctd_link+'_op_0'] # option_slctd_link = "17225"
    #dff = dff[dff["Year"] == option_slctd_link]
    #dff = dff[dff["Affected by"] == "Varroa_mites"]
    
    # fig = dff.plot(title = 'Discharge at USGS '+option_slctd_link,
    #                ylabel = 'Discharge (m3/s)')
    
    fig = px.line(x=dff.index,
                  y=dff.values,
                  labels = dict(x="Time", y="Discharge (m3/s)"),
                  title = 'Discharge at HLM '+option_slctd_link)
    
    # Plotly Express
    # fig = px.choropleth(
    #     data_frame=dff,
    #     locationmode='USA-states',
    #     locations='state_code',
    #     scope="usa",
    #     color='Pct of Colonies Impacted',
    #     hover_data=['State', 'Pct of Colonies Impacted'],
    #     color_continuous_scale=px.colors.sequential.YlOrRd,
    #     labels={'Pct of Colonies Impacted': '% of Bee Colonies'},
    #     template='plotly_dark'
    # )

    # Plotly Graph Objects (GO)
    # fig = go.Figure(
    #     data=[go.Choropleth(
    #         locationmode='USA-states',
    #         locations=dff['state_code'],
    #         z=dff["Pct of Colonies Impacted"].astype(float),
    #         colorscale='Reds',
    #     )]
    # )
    #
    # fig.update_layout(
    #     title_text="Bees Affected by Mites in the USA",
    #     title_xanchor="center",
    #     title_font=dict(size=24),
    #     title_x=0.5,
    #     geo=dict(scope='usa'),
    # )

    return container, fig

# Output of Graph
@app.callback(Output('Elk_graph', 'figure'),
              [Input(component_id='slct_link', component_property='value')])

def update_figure(selected_link):
    df_sub = DF_xy
    
    # Create figure
    locations=[go.Scattermapbox(
                    lon = df_sub['X'],
                    lat = df_sub['Y'],
                    mode='markers',
                    marker={'color' :'red'},
                    unselected={'marker' : {'opacity':1}},
                    selected={'marker' : {'opacity':0.5, 'size':25}},
                    hoverinfo='text',
                    hovertext=df_sub['STANAME']
                    #customdata=df_sub['website']
    )]

    # Return figure
    return {
        'data': locations,
        'layout': go.Layout(
            uirevision= 'foo', #preserves state of figure/map after callback activated
            clickmode= 'event+select',
            hovermode='closest',
            hoverdistance=2,
            title=dict(text="USGS Stations at Elkhorn",font=dict(size=20, color='green')),
            mapbox=dict(
                accesstoken=mapbox_access_token,
                #bearing=25,
                #style='light',
                center=dict(
                    lat=41.9952778,
                    lon=-97.0580556
                ),
                #pitch=40,
                zoom=6
            ),
        )
    }


if __name__ == '__main__':
    app.run_server(debug=True, use_reloader=False)


#%% EXTRAS


# Connect the Plotly graphs with Dash Components
# @app.callback(
#     [Output(component_id='output_container', component_property='children'),
#      Output(component_id='sim_hydrograph', component_property='figure')],
#     [Input(component_id='slct_link', component_property='value')]
# )
# def update_graph(option_slctd_link):
#     print(option_slctd_link)
#     print(type(option_slctd_link))

#     container = "The link chosen by user is: {}".format(option_slctd_link)

#     dff = Qdf.copy()
#     dff = Qdf[option_slctd_link+'_op_0'] # option_slctd_link = "17225"
    
#     fig = px.line(x=dff.index,
#                   y=dff.values,
#                   labels = dict(x="Time", y="Discharge (m3/s)"),
#                   title = 'Discharge at HLM '+option_slctd_link)

#     return container, fig

# Output of Graph
# @app.callback(Output('Elk_graph', 'figure'))#,
#               # [Input(component_id='slct_link', component_property='value')])

# def update_figure(selected_link):
#     df_sub = DF_xy
    
#     # Create figure
#     locations=[go.Scattermapbox(
#                     lon = df_sub['X'],
#                     lat = df_sub['Y'],
#                     mode='markers',
#                     marker={'color' :'red'},
#                     unselected={'marker' : {'opacity':1}},
#                     selected={'marker' : {'opacity':0.5, 'size':25}},
#                     hoverinfo='text',
#                     hovertext=df_sub['STANAME'],
#                     customdata=df_sub['STAID']
#     )]

#     # Return figure
#     return {
#         'data': locations,
#         'layout': go.Layout(
#             uirevision= 'foo', #preserves state of figure/map after callback activated
#             clickmode= 'event+select',
#             hovermode='closest',
#             hoverdistance=2,
#             title=dict(text="USGS Stations at Elkhorn",font=dict(size=20, color='green')),
#             mapbox=dict(
#                 accesstoken=mapbox_access_token,
#                 #bearing=25,
#                 #style='light',
#                 center=dict(
#                     lat=41.9952778,
#                     lon=-97.0580556
#                 ),
#                 #pitch=40,
#                 zoom=6
#             ),
#         )
#     }



