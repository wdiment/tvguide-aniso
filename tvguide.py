"""
Author - William Diment - CIRES IT 


"""

import dash
from dash import dcc
from dash import html
from dash import dash_table
from dash.dependencies import Input, Output, State
from sys import stdout 
import datetime
import appSelectionOptions
import os
import csv
import subprocess
import pandas as pd
import numpy
import json
import urllib
import json
import vector_functions
import plotly.express as px

#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
ciresitMapboxToken = os.environ.get('MAPBOX_TOKEN', '')

#app.css.append_css({"external_url": "static/stylesheet.css"})

app = dash.Dash(__name__)
app.title = 'TVGuide Aniso'
app.scripts.config.server_locally = True
server = app.server
app.config['suppress_callback_exceptions'] = True

config = {"modeBarButtonsToAdd":[['drawline']] }
tab_height = '25px'
tab_style = {'height': tab_height, 'padding': '0'}
tablet_style = {'line-height': tab_height, 'padding': '0'}


degree_selector = [{'label': '1', 'value':'1'}, {'label': '2', 'value':'2'}, {'label': '3', 'value':'3'}]

phase_selector = [{'label': 'Whole', 'value':'0'}, {'label': 'Quartz', 'value':'1'}, {'label': 'Feldspar', 'value':'2'}, {'label': 'Biotite', 'value':'3'}, {'label': 'Muscovite', 'value':'4'}, {'label': 'Chlorite', 'value':'5'}, {'label': 'Hornblende', 'value':'6'}, {'label': 'Garnet', 'value':'7'}, {'label': 'Pyrite', 'value':'8'}]

sampleN_selector = [{'label': '1', 'value':'1'}, {'label': '2', 'value':'2'}, {'label': '3', 'value':'3'}, {'label': '7', 'value':'7'}, {'label': '8', 'value':'8'}, {'label': '9', 'value':'9'}, {'label': '11', 'value':'11'}, {'label': '13', 'value':'13'}, {'label': '14', 'value':'14'} , {'label': '15', 'value':'15'} , {'label': '21', 'value':'21'}, {'label': '22', 'value':'22'} , {'label': '27', 'value':'27'}]

tensor_selector = [{'label': i, 'value': i} for i in range(1,97)]

average_selector = [{'label': i, 'value': i} for i in range(2,6)]



dataButtonStyle={'padding':'5x', 'text-align': 'center', 'line-height':'5'}

velocity_tab_labels = [
  ('Vp', 'tab-2'),
  ('Vs1', 'tab-3'),
  ('Vs2', 'tab-4'),
  ('VP/VS1', 'tab-13'),
  ('Vs1 Polarization & splitting time', 'tab-1'),
  ('3D Quiver Plot/VS1', 'tab-5'),
  ('VP Fig 3D Plot', 'tab-6'),
  ('VS1 3D Plot', 'tab-7'),
  ('VS2 3D Plot', 'tab-8'),
  ('VP/VS1 3D Plot', 'tab-9'),
  ('Back Azimuthal Plot', 'tab-10'),
  ('Radial Plots', 'tab-11'),
]

fold_tab_labels = [
  ('Vp', 'fold-model-tab-2'),
  ('Vs1', 'fold-model-tab-3'),
  ('Vs2', 'fold-model-tab-4'),
  ('VpVs1', 'fold-model-tab-5'),
  ('Vs1 Polarization & splitting time', 'fold-model-tab-1'),
  ('3D Quiver Plot/VS1', 'fold-model-tab-6'),
  ('VP Fig 3D Plot', 'fold-model-tab-7'),
  ('VS1 3D Plot', 'fold-model-tab-8'),
  ('VS2 3D Plot', 'fold-model-tab-9'),
  ('VP/VS1 3D Plot', 'fold-model-tab-10'),
]


navbar = html.Nav(
          className="top-bar fixed",
          children=[
            html.Div([
              html.A(
                  href="https://cires.colorado.edu",
                  children=html.Img(src='assets/cireslogo.png', className='navbar-logo')
                  ),
              html.A(
                  href="https://www.noaa.gov",
                  children=html.Img(src='assets/noaalogo.png', className='navbar-logo')
                  )],
              style={'background-color': 'white',
                     'width': '600px',
                     'text-align': 'center',
                     'float': 'right',
                     'margin-right': '-3px',
                     'margin-top': '-5px',
                     'border': '3px solid rgb(175, 221, 246)',
                     'border-radius': '5px'},
              className='row'),
             # End Sponser Logos


        html.Button(
          children="Database",
          type='button',
          n_clicks=0,
          title='Click to show database. Click again to hide database',
          id="databaseButton",
          style={'height': '45px',
                 'padding': '9px',
                 'background-color': 'rgb(175, 221, 246)',
                 'border-radius': '1px',
                 'font-family': 'Arial',
                 'font-size': '10px',
                 'margin-top': '-5px',
                 'float': 'left',
                 'margin-left': '-5px'}),
        html.Button(
          children="Enter Your Own Tensor",
          type='button',
          n_clicks=0,
          title='Click to enter your own tensor. Click again to hide input area.',
          id="addRowButton",
          style={'height': '45px',
                 'padding': '9px',
                 'background-color': 'rgb(175, 221, 246)',
                 'border-radius': '1px',
                 'font-family': 'Arial',
                 'font-size': '10px',
                 'margin-top': '-5px',
                 'float': 'left',
                 'margin-left': '-5px'}),               
                 ],
          style={'position': 'fixed','top': '0px', 'left': '0px',
                 'background-color': 'rgb(0, 98, 174)', 'height': '50px',
                 'width': '100%', 'zIndex': '9999'
                 },
)
                 # End Toggle Buttons
body = html.Div([
  html.Hr(),
      # Title
    html.Div([
        html.H3('TV Guide - UNDER CONSTRUCTION'),
        html.H4('Interactive elastic tensor visualization tool and database'),
        html.Div(id='desc_div', children=[
          dcc.Markdown(id='description',
          children=[],
          style={'text-align': 'center',
                'width': '70%',
                'margin': '0px auto'}
        )],
        style={'text-align': 'center',
                'margin': '0 auto',
                'width' :'100%'}
        ),
          html.Hr()],
            style={'font-weight': 'bolder',
                   'text-align': 'center',
                   'font-size': '50px',
                   'font-family': 'Times New Roman',
                    'margin-bottom': '50',
                    'margin-top': '100'
                    }),
                   # End Title         
     # Options
    html.Div(id='unit_selection',
    children=[
              dcc.Dropdown(
              id='unit_selection_dropdown',
              searchable=False,
              placeholder='Select Units',
              options=[{'label':'Pascal', 'value':'1'}, {'label':'Gigapascal', 'value':'2'}],
              style={'width':'213px'}
              )
    ]),
    html.Div(id='dataOptions',
              children=[
                html.Button(
                  id='calculate_velocities_button',
                  children="Visualize Velocities",
                  title = "Click here to visualize the velocities",
                  type='button',
                  n_clicks=0,
                ),
                html.Button(
                  id='average_tensors_button',
                  children="Average Tensors",
                  title = "Click here to average tensors",
                  type='button',
                  n_clicks=0,
                ),
                html.Button(
                  id='calculate_fold_model_button',
                  children="Calculate Fold Model",
                  title = "Click here to calculate the fold model for the tensors",
                  type='button',
                  n_clicks=0,
                ),
                html.Button(
                  id='calculate_decomp_button',
                  children="Decompose Tensors",
                  title = "Click here to decompose the tensors",
                  type='button',
                  n_clicks=0,
                ),


                html.Div(id='test-div')],
               style={'text-align': 'center','position':'relative','float':'left','display':'inline'}                 
               ), 
    html.Br(),
    html.Br(),
    html.Div(id='tensor_average_options', children=[

              html.Div(id='number_of_tensors_to_average_div', children=[
                dcc.Dropdown(
                  id='number_of_tensors_to_average_dropdown',
                searchable=False,
                clearable=False,
                  placeholder = 'Select number of tensors to average',
                  options=average_selector,
                  value = 2,
                  style={'height': '30px', 'width': '165px'}
                ),
              ]
              
              ),
              dcc.Dropdown(
                id='select_average_tensor1',
                searchable=False,
                clearable=False,
                options=tensor_selector,
                value=1,
                placeholder='Select Tensor 1...',
                style={'height': '30px', 'width': '165px'}
                ),
              dcc.Dropdown(
                id='select_average_tensor2',
                searchable=False,
                clearable=False,
                options=tensor_selector,
                placeholder='Select Tensor 2...',
                value=1,
                style={'height': '30px', 'width': '165px'}
                ),
              dcc.Dropdown(
                id='select_average_tensor3',
                searchable=False,
                clearable=False,
                options=tensor_selector,
                placeholder='Select Tensor 3...',
                value=1,
                style={'height': '30px', 'width': '165px', 'display':'none'}
                ),
              dcc.Dropdown(
                id='select_average_tensor4',
                searchable=False,
                clearable=False,
                options=tensor_selector,
                placeholder='Select Tensor 4...',
                value=1,
                style={'height': '30px', 'width': '165px', 'display':'none'}
                ),
              dcc.Dropdown(
                id='select_average_tensor5',
                searchable=False,
                clearable=False,
                options=tensor_selector,
                placeholder='Select Tensor 5...',
                value=1,
                style={'height': '30px', 'width': '165px', 'display':'none'}
                ),
                html.Button(id='submit_average_tensors_button',
                n_clicks=0,
                title=('Submit the selected tensor options'),
                children='Submit Average Tensors Options',
                type='button',
                style={'height': '37px', 'width': '325px'}
                ),

    ], style={'display':'none'}),

    html.Div(id='tensor_decomp_options', children=[
      html.Br(),
              dcc.Dropdown(
                id='select_decomp_tensor1',
                searchable=False,
                clearable=False,
                options=tensor_selector,
                placeholder='Select Tensor...',
                value = 1,
                style={'height': '30px', 'width': '165px'}
                ),
                html.Button(id='submit_decomp_tensors_button',
                n_clicks=0,
                title=('Submit the selected tensor options'),
                children='Submit Decomposition',
                type='button',
                style={'height': '37px', 'width': '325px'}
                )
    ], style={'display':'none'}),
    html.Br(),
    html.Div(id='decomp-breakdown-displays',children=[
                html.Button(
                children="Decomposition Breakdown",
                type='button',
                n_clicks=0,
                id="decomp_breakdown"),   
                html.Button(
                children="Isotropic Tensors",
                type='button',
                n_clicks=0,
                id="decomp_iso"),
                html.Button(
                children="Hexagonal Tensors",
                type='button',
                n_clicks=0,
                id="decomp_hexag"),
                html.Button(
                children="Orthorhombic Tensors",
                type='button',
                n_clicks=0,
                id="decomp_ortho"),   
    ], style={'display':'none'}),

    html.Div(id='averaged_tensors_div', children=[

      html.Div(id='weight_selection_div', children=[
        dcc.Input(
          id='weight_selection_1',
          type='number',
          min=1,
          value=1,
          placeholder='Weight 1',
          style={'height': '35px', 'width': '165px', 'display':'inline-block'}
          ),
        dcc.Input(
          id='weight_selection_2',
          type='number',
          value=1,
          min=1,
          placeholder='Weight 2',
          style={'height': '35px', 'width': '165px', 'display':'inline-block'}
          ),
        dcc.Input(
          id='weight_selection_3',
          type='number',
          min=1,
          value=1,
          placeholder='Weight 3',
          style={'height': '35px', 'width': '165px', 'display':'inline-block'}
          ),
        dcc.Input(
          id='weight_selection_4',
          type='number',
          min=1,
          placeholder='Weight 4',
          value=1,
          style={'height': '35px', 'width': '165px', 'display':'inline-block'}
          ),
        dcc.Input(
          id='weight_selection_5',
          type='number',
          min=1,
          value=1,
          placeholder='Weight 5',
          style={'height': '35px', 'width': '165px', 'display':'inline-block'}
          )], style={'margin-left':'165px'}),   
      dcc.Input(
        id='v_ave_1',
        type='text',
        placeholder='V_ave 1',
        value='',
        style={'width':'75%', 'margin-top':'15px'}
      ),
      dcc.Textarea(id='text_area_averaged_tensor',
      value=vector_functions.averagingTensorDescription(),
      style={'width':'1009px', 'resize':'none'},
      disabled=True,
      readOnly=True,
      draggable=False),
      ], style={'display':'none'}),
      

    

    html.Br(),
    html.Div(id='tensor_selection_options',
             children=[
                 dcc.Dropdown(
                  id='select_tensor',
                searchable=False,
                clearable=False,
                  options=tensor_selector,
                  placeholder='Select Tensor...',
                  value = 1,
                  style={'height': '30px', 'width': '165px'}
                  ),
                html.Button(id='submit_button',
                n_clicks=0,
                title=('Submit the selected tensor options'),
                children='Submit Tensor',
                type='button',
                style={'height': '37px', 'width': '165px'}
                ),
            html.Div(id='calculate_velocities_description',
            children=[
            dcc.Textarea(id='text_area_calculate_velocities',
            value=vector_functions.calculateVelocitiesDescription(),
            style={'width':'330px', 'height':'110px', 'resize':'none'},
            disabled=True,
            readOnly=True,
            draggable=False),
            ],
            style={'margin-left':'-330px', 'margin-top':'50px'}
            ),                       
    ],
     style={'margin-top':'50px', 'display':'flex'}
    ),
    html.Br(),
    html.Div(id='tensors_not_filled', children=['You have not correctly filled in the number of tensors - you need 22. This will result in a run time error - please amend your entry.'],
    style={'font-style': 'oblique', 'color':'black', 'text-decoration':'underline', 'display':'none'}
  
    ),
    html.Div(id='tensor_input_table', children=[
      html.Br(),
      dcc.Input(
        id='custom_tensor_list',
        type='text',
        placeholder='Enter your tensors here!',
        value='',
        style={'width':'800px'}
      ),
      html.Button(id='addTensorToDatabase',
        n_clicks=0,
        title=('Submit your own input tensor'),
        children='Submit Input Tensor',
        type='button',
        ),
      html.Br(),
      dcc.Textarea(id='text_area_input_tensor',
      value='Input your own tensors here by copy and pasting into the field above as comma separated elements. Click on Submit Input Tensor when done with the tensor. You can input multiple tensors (one by one), once you are done click on Enter Your Own Tensor at the top to remove these displays',
      style={'width':'1009px', 'resize':'none'},
      draggable=False,
      disabled=True),
    ], style={'display':'none'}),
    html.H4(id='database-table-title', children='Database', style={'display':'none'}),
    html.Div(id='database_table', children=vector_functions.print_csv_table()[0], style={'display':'none'}),
    html.Div(id='database_dataframe', children=vector_functions.print_csv_table()[1], style={'display':'none'}),
    html.Div(id='database_dataframe_user_input', children=vector_functions.print_csv_table()[1], style={'display':'none'}),
    html.Div(id='database_dataframe_user_averaging_input', children=vector_functions.print_csv_table()[1], style={'display':'none'}),
    html.Div(id='decomp_breakdown_table', children=vector_functions.print_breakdown(), style={'display':'none'}),
    html.Div(id='decomp_iso_table', children='hello', style={'display':'none'}),
    html.Div(id='decomp_hexag_table', children='hello', style={'display':'none'}),
    html.Div(id='decomp_ortho_table', children='hello', style={'display':'none'}),
    html.Div(id='decomp_component_export', children=[dcc.Textarea(
      id='decomp_component_export_text',
      value='',
      placeholder='Components'
    )], style={'display':'none'}),
    ###put holder here - 7/27/2023

    html.Br(),
        html.Div([
          html.H4(id='tool-title-display', children='Current Tool: Visualize Velocities'),
          dcc.Tabs(id='tabs-div', value='tab-2', style={'display':'none'},   children=[
            dcc.Tab(label=label, value=value, className='custom-tab', selected_className='custom-tab--selected')
            for label, value in velocity_tab_labels
          ]),
          dcc.Tabs(id='fold-model-tabs-div', value='fold-model-tab-1', style={'display':'none'},   children=[
            dcc.Tab(label=label, value=value, className='custom-tab', selected_className='custom-tab--selected')
            for label, value in fold_tab_labels
          ]),
       html.Div(id='json_hidden_table', children=[
           dcc.Loading(
              id='loading-2',
              type='default',
              children=html.Div(id='loading-output-2')
            )
       ], 
       title='', style={'margin-top':'10px'}),
          html.Div(id='tab-content-display', children=[
            dcc.Loading(
              id='loading-1',
              type='default',
              children=html.Div(id='loading-output-1')
            ),
            dcc.Graph(id='graph-display', figure={})
          ], style={'margin-top':'50px','overflow':'scroll'}),
          html.Div(id='fold-model-display', children=[], style={'margin-top':'50px','overflow':'scroll'})
          
          ], style={'text-align': 'center', 'top':'15px', 'max-height':'1000px', 'min-width':'1600px'}),
        
       ],
       
       className='ten columns offset-by-one', style={'min-width':'1200px'})
       
       # End Static Elements


@app.callback(
  Output(component_id='v_ave_1', component_property='value'),
  [Input(component_id='submit_average_tensors_button', component_property='n_clicks'),
  Input(component_id='number_of_tensors_to_average_dropdown', component_property='value')],
  [State(component_id='select_average_tensor1', component_property='value'),
  State(component_id='select_average_tensor2', component_property='value'),
  State(component_id='select_average_tensor3', component_property='value'),
  State(component_id='select_average_tensor4', component_property='value'),
  State(component_id='select_average_tensor5', component_property='value'),

  State(component_id='weight_selection_1', component_property='value'),
  State(component_id='weight_selection_2', component_property='value'),
  State(component_id='weight_selection_3', component_property='value'),
  State(component_id='weight_selection_4', component_property='value'),
  State(component_id='weight_selection_5', component_property='value'),
  State(component_id='database_dataframe_user_averaging_input', component_property='children')],
)
def average_tensors(n_clicks, numberOfTensorsToAverage, tensor1, tensor2, tensor3, tensor4, tensor5, weight1, weight2, weight3, weight4, weight5, averagingDatabaseFrame):
  averagingDatabaseFrame = averagingDatabaseFrame

  display_tensor_style = {'display':'block'}
  if n_clicks == 0:
    return 'Select tensors to average!'
  if numberOfTensorsToAverage == 2:
    weight = [weight1, weight2]   
    if tensor1 < 96 and tensor2 < 96:
     V_ave, R_ave = vector_functions.tv_averaging(tensor1, tensor2, None, None, None, False, None, weight)
    if tensor1 >= 96 or tensor2 >= 96:
      V_ave, R_ave = vector_functions.tv_averaging(tensor1, tensor2, None, None, None, True, averagingDatabaseFrame, weight)

  if numberOfTensorsToAverage == 3:
    weight = [weight1, weight2, weight3]   
    if tensor1 < 96 and tensor2 < 96 and tensor3 < 96:
     V_ave, R_ave = vector_functions.tv_averaging(tensor1, tensor2, tensor3, None, None, False, None, weight)
    if tensor1 >= 96 or tensor2 >= 96 or tensor3 >= 96:
      V_ave, R_ave = vector_functions.tv_averaging(tensor1, tensor2, tensor3, None, None, True, averagingDatabaseFrame, weight)    

  if numberOfTensorsToAverage == 4:
    weight = [weight1, weight2, weight3, weight4]  
    if tensor1 < 96 and tensor2 < 96 and tensor3 < 96 and tensor4 < 96:
     V_ave, R_ave = vector_functions.tv_averaging(tensor1, tensor2, tensor3, tensor4, None, False, None, weight)
    if tensor1 >= 96 or tensor2 >= 96 or tensor3 >= 96 or tensor4 >= 96:
      V_ave, R_ave = vector_functions.tv_averaging(tensor1, tensor2, tensor3, tensor4, None, True, averagingDatabaseFrame, weight)   

  if numberOfTensorsToAverage == 5:
    weight = [weight1, weight2, weight3, weight4, weight5]
    if tensor1 < 96 and tensor2 < 96 and tensor3 < 96 and tensor4 < 96 and tensor5 < 96:
     V_ave, R_ave = vector_functions.tv_averaging(tensor1, tensor2, tensor3, tensor4, tensor5, False, None, weight)
    if tensor1 >= 96 or tensor2 >= 96 or tensor3 >= 96 or tensor4 >= 96 or tensor5 >= 96:
      V_ave, R_ave = vector_functions.tv_averaging(tensor1, tensor2, tensor3, tensor4, tensor5, True, averagingDatabaseFrame, weight)   


  v_ave_string = "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(V_ave[0],V_ave[1],V_ave[2],V_ave[3],V_ave[4],V_ave[5],V_ave[6],V_ave[7],V_ave[8],V_ave[9],V_ave[10],V_ave[11],V_ave[12],V_ave[13],V_ave[14],V_ave[15],V_ave[16],V_ave[17],V_ave[18],V_ave[19],V_ave[20],V_ave[21],V_ave[22],V_ave[23],V_ave[24],V_ave[25],V_ave[26],V_ave[27],V_ave[28],V_ave[29],V_ave[30],V_ave[31],V_ave[32],V_ave[33],V_ave[34],V_ave[35])
  return v_ave_string

# Maps each velocity tab to its plot type and whether it is a 3D plot
TAB_PLOT_TYPES = {
    'tab-1': ('Quiver', False),
    'tab-2': ('VP', False),
    'tab-3': ('VS1', False),
    'tab-4': ('VS2', False),
    'tab-13': ('VPVS1', False),
    'tab-5': ('3DQuiver', True),
    'tab-6': ('3DVP', True),
    'tab-7': ('3DVS1', True),
    'tab-8': ('3DVS2', True),
    'tab-9': ('3DVPVS1', True),
    'tab-10': ('BackAzimuthal', False),
    'tab-11': ('RadialPlots', False),
}

FOLD_TAB_PLOT_TYPES = {
    'fold-model-tab-1': ('Quiver', False),
    'fold-model-tab-2': ('VP', False),
    'fold-model-tab-3': ('VS1', False),
    'fold-model-tab-4': ('VS2', False),
    'fold-model-tab-5': ('VPVS1', False),
    'fold-model-tab-6': ('3DQuiver', True),
    'fold-model-tab-7': ('3DVP', True),
    'fold-model-tab-8': ('3DVS1', True),
    'fold-model-tab-9': ('3DVS2', True),
    'fold-model-tab-10': ('3DVPVS1', True),
}

FOLD_GRAPH_STYLE = {'text-align': 'center', 'top': '15px', 'max-height': '1000px', 'min-width': '1600px'}

# Database tensors are numbered 1-96; user-input tensors start at 97
FIRST_USER_TENSOR = 97


@app.callback(
[Output(component_id='graph-display', component_property='figure'),
Output(component_id='fold-model-display', component_property='children'),
Output(component_id='loading-output-1', component_property='children'),
Output(component_id='select_decomp_tensor1', component_property='value')],
[Input(component_id='fold-model-tabs-div', component_property='value'),
Input(component_id='tabs-div', component_property='value'),
Input(component_id='submit_button', component_property='n_clicks')],[
State(component_id='select_tensor', component_property='value'),
State(component_id='database_dataframe_user_input', component_property='children')
]
)
def render_content(folds_tab, tab, button_nclicks, tensor, userInputDataFrame):
  triggered = [p['prop_id'] for p in dash.callback_context.triggered]
  changed_id = triggered[0] if triggered else '.'
  print('this is changed id in the render_content tool', changed_id, 'this is the tab', tab)

  tensor = int(tensor)
  user_input = tensor >= FIRST_USER_TENSOR

  if changed_id == 'fold-model-tabs-div.value' and folds_tab in FOLD_TAB_PLOT_TYPES:
    plotType, is3D = FOLD_TAB_PLOT_TYPES[folds_tab]
    fig = vector_functions.tv_fold_model(
        tensor, user_input,
        None if is3D else plotType,
        plotType if is3D else None,
        userInputDataFrame)
    foldGraph = html.Div(children=[dcc.Graph(figure=fig)], style=FOLD_GRAPH_STYLE)
    return {}, foldGraph, '', tensor

  # Default branch: velocity plot for the current tab (also covers initial page load)
  plotType, is3D = TAB_PLOT_TYPES.get(tab, TAB_PLOT_TYPES['tab-2'])
  fig = vector_functions.calculate_tensor_symmetries(
      tensor, user_input,
      None if is3D else plotType,
      plotType if is3D else None,
      userInputDataFrame)
  return fig, '', '', tensor


app.layout = html.Div([navbar,body])

@app.callback(
  [Output(component_id='database_table', component_property='style'),
  Output(component_id='database-table-title', component_property='children')],
  [Input(component_id='databaseButton', component_property='n_clicks')])
def display_database_table(n_clicks):
  print('hello')
  sampleNumber = 1

  databaseTableTitle = ''
  databaseTableTitleStyle = {'display':'none'}
  

  changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
  
  if n_clicks % 2 == 0:
    database = {'display': 'none', 'text-align': 'center'}
  else:
    database = {'display': 'block', 'text-align': 'center'}
    breakdownDatabase = {'display': 'none', 'text-align': 'center'}

  return database,  databaseTableTitle



@app.callback(
  Output(component_id='json_hidden_table', component_property='title'),
  Output(component_id='loading-output-2', component_property='children'),
  Output(component_id='decomp_breakdown_table', component_property='style'),
  Output(component_id='decomp_breakdown_table', component_property='children'),
  Output(component_id='decomp_iso_table', component_property='style'),
  Output(component_id='decomp_iso_table', component_property='children'),
  Output(component_id='decomp_hexag_table', component_property='style'),
  Output(component_id='decomp_hexag_table', component_property='children'),
  Output(component_id='decomp_ortho_table', component_property='style'),
  Output(component_id='decomp_ortho_table', component_property='children'),
  Output(component_id='tensor_average_options', component_property='style'),
  Output(component_id='tensor_selection_options', component_property='style'),
  Output(component_id='averaged_tensors_div', component_property='style'),
  Output(component_id='tabs-div', component_property='style'),
  Output(component_id='fold-model-tabs-div', component_property='style'),
  Output(component_id='tab-content-display', component_property='style'),
  Output(component_id='fold-model-display', component_property='style'),
  Output(component_id='tensor_decomp_options', component_property='style'),
  Output(component_id='decomp-breakdown-displays', component_property='style'),
  Output(component_id='tool-title-display', component_property='children'),
  Output(component_id='database-table-title', component_property='style'),
  Output(component_id='graph-display', component_property='style'),
  Output(component_id='decomp_component_export', component_property='style'),
  Output(component_id='db_csv_table', component_property='data'),
  Output(component_id='tensor_input_table', component_property='style'),
  Output(component_id='select_tensor', component_property='options'),
  Output(component_id='select_decomp_tensor1', component_property='options'),
  Output(component_id='select_average_tensor1', component_property='options'),
  Output(component_id='select_average_tensor2', component_property='options'),
  Output(component_id='select_average_tensor3', component_property='options'),
  Output(component_id='select_average_tensor4', component_property='options'),
  Output(component_id='select_average_tensor5', component_property='options'),
  Output(component_id='tensors_not_filled', component_property='style'),
  Output(component_id='database_table', component_property='children'),
  Output(component_id='database_dataframe', component_property='children'),
  Output(component_id='database_dataframe_user_input', component_property='children'),
  Output(component_id='database_dataframe_user_averaging_input', component_property='children'),
  Output(component_id='select_tensor', component_property='value'),

  Output(component_id='select_average_tensor3', component_property='style'),
  Output(component_id='select_average_tensor4', component_property='style'),
  Output(component_id='select_average_tensor5', component_property='style'),

  Output(component_id='weight_selection_3', component_property='style'),
  Output(component_id='weight_selection_4', component_property='style'),
  Output(component_id='weight_selection_5', component_property='style'),

  Output(component_id='number_of_tensors_to_average_dropdown', component_property='value'),
  [Input(component_id='submit_decomp_tensors_button', component_property='n_clicks'),
  Input(component_id='decomp_breakdown', component_property='n_clicks'),
  Input(component_id='decomp_iso', component_property='n_clicks'),
  Input(component_id='decomp_hexag', component_property='n_clicks'),
  Input(component_id='decomp_ortho', component_property='n_clicks'),
  Input(component_id='average_tensors_button', component_property='n_clicks'),
  Input(component_id='calculate_velocities_button', component_property='n_clicks'),
  Input(component_id='calculate_fold_model_button', component_property='n_clicks'),
  Input(component_id='calculate_decomp_button', component_property='n_clicks'),
  Input(component_id='addRowButton', component_property='n_clicks'),
  Input(component_id='addTensorToDatabase', component_property='n_clicks'),
  Input(component_id='database_table', component_property='children'),
  Input(component_id='database_dataframe', component_property='children'),
  Input(component_id='number_of_tensors_to_average_dropdown', component_property='value')
  ],
  [State(component_id='select_decomp_tensor1', component_property='value'),
  State(component_id='select_tensor', component_property='options'),
  State(component_id='db_csv_table', component_property='data'),
  State(component_id='db_csv_table', component_property='columns'),
  State(component_id='custom_tensor_list', component_property='value')]
)
def generate_decomposition_switch_tools(n_clicks, decomp_breakdown, decomp_iso, decomp_hex, decomp_ortho,  average_button, calculate_velocities, calculate_fold_model, decomp, addDatabaseRowButtonClicks, addRowToDatabaseClicks, csvDatabase, databaseDataframe, numberOfAverageTensorsSelect, tensor1, selectDatabaseTensorOptions, databaseRows, databaseColumns, inputTensorList):
  print('hello')
  changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
  print("This is the changed ID - generated in decomp/styling tool", changed_id, type(changed_id))
  
  if(tensor1 <= 96):
    Line, isoComponent, hexagonalComponent, orthoComponent = vector_functions.decomposition(tensor1)
  if(tensor1 > 96):
    Line, isoComponent, hexagonalComponent, orthoComponent = vector_functions.decomposition(tensor1, True, 'None', databaseDataframe)
  
  breakdownFrame = vector_functions.print_breakdown(Line)
  isoFrame = vector_functions.print_iso_decomp(isoComponent[0])
  hexagFrame = vector_functions.print_hexag_decomp(hexagonalComponent[0])
  orthoFrame = vector_functions.print_ortho_decomp(orthoComponent[0])
  userInputDataFrame = databaseDataframe
  numberOfAverageTensorsSelect = numberOfAverageTensorsSelect

  
  breakdownDatabase = {'display': 'none', 'text-align': 'center'}
  isoDatabase = {'display': 'none', 'text-align': 'center'}
  hexagDatabase = {'display': 'none', 'text-align': 'center'}
  orthoDatabase = {'display': 'none', 'text-align': 'center'}

  tensor_selection_options_display = {'display':'flex'}
  display_tensor_average_options = {'display':'none'}
  display_tensor_average_cells = {'display':'none'}
  display_tensor_average_cells = {'display':'none'}
  display_decomp_options = {'display':'none'}
  display_decomp_database_options = {'display':'none'}

  third_average_option = {'display':'none'}
  fourth_average_option = {'display':'none'}
  fifth_average_option = {'display':'none'}


  third_weight_selection_option = {'display':'none'}
  fourth_weight_selection_option = {'display':'none'}
  fifth_weight_selection_option = {'display':'none'}

  graphDisplay = {'display':'flex', 'max-width':'1800px', 'margin':'auto'}
  decompComponentDisplay = {'display':'none'}
  decompOutput = ''
  
  toolTitleDisplay = 'Current Tool: Visualize Velocities'

  tabs_div_style = {'display':''}
  fold_model_tabs_div = {'display':'none'}
  tab_content_display = {'margin-top':'50px'}
  fold_model_display = {'display':'none'}

  if changed_id == 'decomp_breakdown.n_clicks':
    breakdownDatabase = {'display': 'block', 'text-align': 'center'}
    isoDatabase = {'display': 'none', 'text-align': 'center'}
    hexagDatabase = {'display': 'none', 'text-align': 'center'}
    orthoDatabase = {'display': 'none', 'text-align': 'center'}

    display_decomp_options = {'display':'flex'}
    display_decomp_database_options = {'display':'flex'}
    databaseTableTitle = 'Breakdown Table'
    databaseTableTitleStyle = {'display':'none'}
    tensor_selection_options_display = {'display':'none'}
    graphDisplay = {'display':'none'}
    tabs_div_style = {'display':'none'}
    fold_model_display = {'display':'none'}
    toolTitleDisplay = 'Current Tool: Decomposition'

    decompComponentDisplay = {'display':'flex'}



  if changed_id == 'decomp_iso.n_clicks':
    isoDatabase = {'display': 'block', 'text-align': 'center'}
    breakdownDatabase = {'display': 'none', 'text-align': 'center'}
    hexagDatabase = {'display': 'none', 'text-align': 'center'}
    orthoDatabase = {'display': 'none', 'text-align': 'center'}

    display_decomp_options = {'display':'flex'}
    display_decomp_database_options = {'display':'flex'}
    databaseTableTitle = 'Isotropic Decomposition Table'
    databaseTableTitleStyle = {'display':'none'}
    tensor_selection_options_display = {'display':'none'}
    graphDisplay = {'display':'none'}
    tabs_div_style = {'display':'none'}
    fold_model_display = {'display':'none'}
    toolTitleDisplay = 'Current Tool: Decomposition'
    decompComponentDisplay = {'display':'flex'}


  if changed_id == 'decomp_hexag.n_clicks':
    hexagDatabase = {'display': 'block', 'text-align': 'center'}
    breakdownDatabase = {'display': 'none', 'text-align': 'center'}
    isoDatabase = {'display': 'none', 'text-align': 'center'}
    orthoDatabase = {'display': 'none', 'text-align': 'center'}

    display_decomp_options = {'display':'flex'}
    display_decomp_database_options = {'display':'flex'}
    databaseTableTitle = 'Hexagonal Decomposition Table'
    databaseTableTitleStyle = {'display':'none'}
    tensor_selection_options_display = {'display':'none'}
    graphDisplay = {'display':'none'}
    tabs_div_style = {'display':'none'}
    fold_model_display = {'display':'none'}
    toolTitleDisplay = 'Current Tool: Decomposition'
    decompComponentDisplay = {'display':'flex'}


  if changed_id == 'decomp_ortho.n_clicks':
    orthoDatabase = {'display': 'block', 'text-align': 'center'}

    hexagDatabase = {'display': 'none', 'text-align': 'center'}
    breakdownDatabase = {'display': 'none', 'text-align': 'center'}
    isoDatabase = {'display': 'none', 'text-align': 'center'}

    display_decomp_options = {'display':'flex'}
    display_decomp_database_options = {'display':'flex'}

    databaseTableTitle = 'Orthrhombic Decomposition Table'
    databaseTableTitleStyle = {'display':'none'}
    tensor_selection_options_display = {'display':'none'}
    graphDisplay = {'display':'none'}
    fold_model_display = {'display':'none'}
    tabs_div_style = {'display':'none'}
    toolTitleDisplay = 'Current Tool: Decomposition'
    decompComponentDisplay = {'display':'flex'}



  databaseTitleDisplayStyle = {'display':'none'}


  if changed_id == 'average_tensors_button.n_clicks' or changed_id == 'number_of_tensors_to_average_dropdown.value':
    display_tensor_average_options = {'display':'flex'}

    tensor_selection_options_display = {'display':'none'}
    display_tensor_average_cells = {'display':'block'}
    display_decomp_options = {'display':'none'}
    display_decomp_database_options = {'display':'none'}
    fold_model_display = {'display':'none'}

    graphDisplay = {'display':'none'}
    tabs_div_style = {'display':'none'}

    tool = 'Averaging Tensors'
    toolTitleDisplay = 'Current Tool: ' + tool
    print(numberOfAverageTensorsSelect)
    if int(numberOfAverageTensorsSelect) == 3:
      third_average_option = {'height': '36px', 'width': '165px', 'display':'inline-block'}
      third_weight_selection_option = {'height': '36px', 'width': '165px', 'display':'inline-block'}

    if int(numberOfAverageTensorsSelect) == 4:
      third_average_option = {'height': '36px', 'width': '165px', 'display':'inline-block'}
      fourth_average_option = {'height': '36px', 'width': '165px', 'display':'inline-block'}
      third_weight_selection_option = {'height': '36px', 'width': '165px', 'display':'inline-block'}
      fourth_weight_selection_option = {'height': '36px', 'width': '165px', 'display':'inline-block'}
    if int(numberOfAverageTensorsSelect) == 5:
      third_average_option = {'height': '30x', 'width': '165px', 'display':'inline-block'}
      fourth_average_option = {'height': '35px', 'width': '165px', 'display':'inline-block'}
      fifth_average_option = {'height': '30x', 'width': '165px', 'display':'inline-block'}
      third_weight_selection_option = {'height': '36px', 'width': '165px', 'display':'inline-block'}
      fourth_weight_selection_option = {'height': '36px', 'width': '165px', 'display':'inline-block'}
      fifth_weight_selection_option = {'height': '36px', 'width': '165px', 'display':'inline-block'}

    
  if changed_id == 'calculate_velocities_button.n_clicks': 
    tensor_selection_options_display = {'display':'flex'}

    display_tensor_average_options = {'display':'none'}
    display_tensor_average_cells = {'display':'none'}
    display_decomp_database_options = {'display':'none'}
    fold_model_display = {'display':'none'}

    tool = 'Visualize Velocities'

    toolTitleDisplay = 'Current Tool: ' + tool

    display_decomp_options = {'display':'none'}

  if changed_id == 'calculate_fold_model_button.n_clicks':
    tensor_selection_options_display = {'display':'flex'}

    display_tensor_average_options = {'display':'none'}
    display_tensor_average_cells = {'display':'none'}
    display_decomp_options = {'display':'none'}
    display_decomp_database_options = {'display':'none'}


    tabs_div_style = {'display':'none'}
    fold_model_tabs_div = {'display':''}
    tab_content_display = {'margin-top':'50px', 'display':'none'}

    fold_model_display = {'display':'flex', 'max-width':'1800px', 'margin':'auto'}

    tool = 'Fold Model'
    toolTitleDisplay = 'Current Tool: ' + tool

    ### BEGINNING ADDING DATABASE ROW OPTIONS
  if addDatabaseRowButtonClicks % 2 != 0 :
    databaseTableDisplay = {'display':'block'}
  else:
    databaseTableDisplay = {'display':'none'}
  
  tensorListAdditions = selectDatabaseTensorOptions
  column_names = ['Paper Reference','tensor #' ,'Rock Type','density (g/cm3)','C11 (Mbar)','C12','C13','C14','C15','C16','C22','C23','C24','C25','C26','C33','C34','C35','C36','C44','C45','C46','C55','C56','C66']
  csvRow = ''
  tensorList= [x for x in range(96,122)]
  numberList = [x for x in range (1,27)]
  numberDict = [{'label': x, 'value':y} for x, y in zip(tensorList, numberList)]
  parseTensors = inputTensorList.replace('\t', ',').replace(' ', ',').split(",")
  tensorNotFilledDisplay = {'font-style': 'oblique', 'color':'black', 'text-decoration':'underline', 'display':'none'}
  
  if parseTensors[0] != '':
    g = float(parseTensors[0])
  if len(parseTensors) == 22:
    c11 = float(parseTensors[1])
    c12 = float(parseTensors[2])
    c13 = float(parseTensors[3])
    c14 = float(parseTensors[4])
    c15 = float(parseTensors[5])
    c16 = float(parseTensors[6])
    c22 = float(parseTensors[7])
    c23 = float(parseTensors[8])
    c24 = float(parseTensors[9])
    c25 = float(parseTensors[10])
    c26 = float(parseTensors[11])
    c33 = float(parseTensors[12])
    c34 = float(parseTensors[13])
    c35 = float(parseTensors[14])
    c36 = float(parseTensors[15])
    c44 = float(parseTensors[16])
    c45 = float(parseTensors[17])
    c46 = float(parseTensors[18])
    c55 = float(parseTensors[19])
    c56 = float(parseTensors[20])
    c66 = float(parseTensors[21])
    listOfTensorsToConfirm = [g, c11, c12, c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66]

  if changed_id == 'addTensorToDatabase.n_clicks':
    
    if len(parseTensors) != 22:
      tensorNotFilledDisplay = {'font-style': 'oblique', 'color':'black', 'text-decoration':'underline', 'display':'block'}
      return 1, '', breakdownDatabase, breakdownFrame, isoDatabase, isoFrame, hexagDatabase, hexagFrame, orthoDatabase, orthoFrame, display_tensor_average_options, tensor_selection_options_display, display_tensor_average_cells, tabs_div_style, fold_model_tabs_div, tab_content_display, fold_model_display, display_decomp_options, display_decomp_database_options, toolTitleDisplay, databaseTitleDisplayStyle, graphDisplay, decompComponentDisplay, databaseRows, databaseTableDisplay, tensorListAdditions, tensorListAdditions, tensorListAdditions, tensorListAdditions, tensorListAdditions, tensorListAdditions,  tensorListAdditions, tensorNotFilledDisplay,  csvDatabase, databaseDataframe, userInputDataFrame, userInputDataFrame, tensor1, third_average_option, fourth_average_option, fifth_average_option, third_weight_selection_option, fourth_weight_selection_option, fifth_weight_selection_option, numberOfAverageTensorsSelect
    
    else:
      tensorNotFilledDisplay = {'font-style': 'oblique', 'color':'black', 'text-decoration':'underline', 'display':'none'}
    userTensorNumber = tensorListAdditions[-1]['label']


    if userTensorNumber in tensorList:
      userTensorNumber = tensorList[(tensorList.index(userTensorNumber)+1)]
      dictAdd = {'label': userTensorNumber, 'value': userTensorNumber}

    tensorListAdditions.append(dictAdd)
    rowAddition = {'Paper Reference': 'User Input', 'tensor #': userTensorNumber, 'Rock Type': '1', 'density (g/cm3)': g, 'C11 (Mbar)': c11, 'C12': c12, 'C13': c13, 'C14':c14, 'C15': c15, 'C16':c16, 'C22': c22, 'C23': c23, 'C24': c24, 'C25': c25, 'C26': c26, 'C33': c33, 'C34': c34, 'C35': c35, 'C36': c36, 'C44': c44, 'C45': c45, 'C46': c46, 'C55': c55, 'C56': c56, 'C66': c66}
    csvDatabase, databaseDataframe = vector_functions.print_csv_table(databaseDataframe, rowAddition)
    userInputDataFrame = databaseDataframe
    tensor1 = userTensorNumber

  if decomp == 0:
      userInputDataFrame = databaseDataframe
      return 1, '', breakdownDatabase, breakdownFrame, isoDatabase, isoFrame, hexagDatabase, hexagFrame, orthoDatabase, orthoFrame, display_tensor_average_options, tensor_selection_options_display, display_tensor_average_cells, tabs_div_style, fold_model_tabs_div, tab_content_display, fold_model_display, display_decomp_options, display_decomp_database_options, toolTitleDisplay, databaseTitleDisplayStyle, graphDisplay, decompComponentDisplay, databaseRows, databaseTableDisplay, tensorListAdditions, tensorListAdditions, tensorListAdditions, tensorListAdditions, tensorListAdditions, tensorListAdditions, tensorListAdditions,  tensorNotFilledDisplay, csvDatabase, databaseDataframe, userInputDataFrame, userInputDataFrame, tensor1, third_average_option, fourth_average_option, fifth_average_option, third_weight_selection_option, fourth_weight_selection_option, fifth_weight_selection_option, numberOfAverageTensorsSelect


  if changed_id == 'calculate_decomp_button.n_clicks':
    
    tensor_selection_options_display = {'display':'none'}
    display_tensor_average_options = {'display':'none'}
    display_tensor_average_cells = {'display':'none'}
    display_decomp_options = {'display':'flex'}
    display_decomp_database_options = {'display':'flex'}
    

    tool = 'Decomposition'
    toolTitleDisplay = 'Current Tool: ' + tool
    databaseTitleDisplayStyle = {'display':'block'}
    graphDisplay = {'display':'none'}
    tabs_div_style = {'display':'none'}

  if changed_id == 'submit_decomp_tensors_button.n_clicks':
    tensor_selection_options_display = {'display':'none'}
    display_tensor_average_options = {'display':'none'}
    display_tensor_average_cells = {'display':'none'}
    display_decomp_options = {'display':'flex'}
    display_decomp_database_options = {'display':'flex'}
    ### TEST 1172024
    breakdownDatabase = {'display': 'block', 'text-align': 'center'}
    

    tool = 'Decomposition'
    toolTitleDisplay = 'Current Tool: ' + tool
    databaseTitleDisplayStyle = {'display':'block'}
    graphDisplay = {'display':'none'}
    tabs_div_style = {'display':'none'}

    dictColumns = ['C11 (Mbar)','C12','C13','C14','C15','C16','C22','C23','C24','C25','C26','C33','C34','C35','C36','C44','C45','C46','C55','C56','C66']
    userTensorNumber = tensorListAdditions[-1]['label']
    if userTensorNumber in tensorList:
      userTensorNumber = tensorList[(tensorList.index(userTensorNumber)+1)]
      dictAdd = {'label': userTensorNumber, 'value': userTensorNumber}
    tensorListAdditions.append(dictAdd)
    rowAddition = {'Paper Reference': 'Decomposed Orthorhombic Tensor', 'tensor #': userTensorNumber, 'Rock Type': '1', 'density (g/cm3)':'3'}
    counter = 0
    key = orthoFrame.data
    key = key[0].keys()
    
    for key in key:
      columnLabel = int(key)-1 
      rowAddition[dictColumns[columnLabel]] = str(orthoFrame.data[0][key])
    csvDatabase, databaseDataframe = vector_functions.print_csv_table(databaseDataframe, rowAddition)
    userInputDataFrame = databaseDataframe
  
    

  return 1, '', breakdownDatabase, breakdownFrame, isoDatabase, isoFrame, hexagDatabase, hexagFrame, orthoDatabase, orthoFrame, display_tensor_average_options, tensor_selection_options_display, display_tensor_average_cells, tabs_div_style, fold_model_tabs_div, tab_content_display, fold_model_display, display_decomp_options, display_decomp_database_options, toolTitleDisplay, databaseTitleDisplayStyle, graphDisplay, decompComponentDisplay, databaseRows, databaseTableDisplay, tensorListAdditions, tensorListAdditions,tensorListAdditions, tensorListAdditions, tensorListAdditions, tensorListAdditions, tensorListAdditions,   tensorNotFilledDisplay, csvDatabase, databaseDataframe, userInputDataFrame, userInputDataFrame, tensor1, third_average_option, fourth_average_option, fifth_average_option, third_weight_selection_option, fourth_weight_selection_option, fifth_weight_selection_option, numberOfAverageTensorsSelect






if __name__=='__main__':
    app.run_server(debug=True)



### matlab command for running the script - matlab -nodisplay -nosplash -nodesktop -r "numVRH=1;numMin=1;run('ChristoffelPlotDEV.m');exit;" 