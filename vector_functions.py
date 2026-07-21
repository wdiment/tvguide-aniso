import os 
import numpy as np
import math
import statistics as stats
import pandas as pd
import plotly as plt
import plotly.figure_factory as ff
import scipy as sp
import math
import dash
from dash import dcc
from dash import html
from dash import dash_table
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.figure_factory import create_quiver
import plotly.express as px
from scipy import interpolate
from scipy import linalg
from scipy import spatial
from iteration_utilities import flatten
from iteration_utilities import deepflatten
import ChristoffelDevParams
np.set_printoptions(threshold=np.inf)
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import plotly.tools

np.seterr(divide='ignore', invalid='ignore')


def calculate_tensor_symmetries(tensor_index, user_input=False, plotType='None', plotType3D='None', userInputDataFrame=None):
    christoffelParameters = ChristoffelDevParams.ChristoffelParameters()
    Res = 3
    pi = math.pi
    aggTensorColumns = christoffelParameters.columnNamesTensorCSV
    sphereColumnNames = [x for x in range(0, 121)]
    csv_tensor_number = [x for x in range(0, 50)]
    user_tensor_number = [x for x in range(97, 147)]

    aggTensorCSVLocation = 'databases/DB_Tensors_Rhos_OrigOrientation_CSV.csv'

    C = []
    
    returnedUserInputDataFrame = returnUserInputDataFrame(userInputDataFrame)

    minNum = christoffelParameters.minNum
    Samps = 3

    tensorDF = pd.read_csv(aggTensorCSVLocation, sep=',', engine='python')

    tensorIndex = int(tensor_index)

    if user_input == True:
        tensorDF = returnedUserInputDataFrame
        tensorList = (tensorDF.iloc[tensorIndex, 4:25])
        Rhos = (tensorDF.iloc[tensorIndex, 3])
    else:
        tensorList = (tensorDF.iloc[tensorIndex-1, 4:25])
        Rhos = (tensorDF.iloc[tensorIndex-1, 3])

    Rhos = float(Rhos)


    C = [[tensorList[0],tensorList[1],tensorList[2],tensorList[3],tensorList[4],tensorList[5]],
             [tensorList[1],tensorList[6],tensorList[7],tensorList[8],tensorList[9],tensorList[10]],
             [tensorList[2],tensorList[7],tensorList[11],tensorList[12],tensorList[13],tensorList[14]],
             [tensorList[3],tensorList[8],tensorList[12],tensorList[15],tensorList[16],tensorList[17]],
             [tensorList[4],tensorList[9],tensorList[13],tensorList[16],tensorList[18],tensorList[19]],
             [tensorList[5],tensorList[10],tensorList[14],tensorList[17],tensorList[19],tensorList[20]]
             ]

    
    
    

    C = np.asarray(C, dtype=float)
    print(C)

    voigtIndices = [[1,6,5], [6,2,4], [5,4,3]]

    Voigt_ijkl = pd.DataFrame(data=voigtIndices)
    Voigt_ijkl = Voigt_ijkl.to_numpy()

    rho = Rhos

    Xx = pd.read_csv('databases/xx.csv', names=sphereColumnNames, sep=',', engine='python')
    Yy = pd.read_csv('databases/yy.csv', names=sphereColumnNames, sep=',', engine='python')
    Zz = pd.read_csv('databases/zz.csv', names=sphereColumnNames, sep=',', engine='python')

    Xx = Xx.fillna(0)
    Yy = Yy.fillna(0)
    Zz = Zz.fillna(0)


    Xx = Xx.to_numpy()
    Yy = Yy.to_numpy()
    Zz = Zz.to_numpy()


    test1 = 0
    test2 = 0

    m = 0
    n = 0

    T = np.zeros((3,3), dtype=float)
    TT = np.zeros((3,3), dtype=float)
    velp = np.zeros((121,121), dtype=float)
    vels1 = np.zeros((121,121), dtype=float)
    vels2 = np.zeros((121,121), dtype=float)

    VpPolxyz = np.zeros((121,121,3), dtype=float)
    Vs1Polxyz = np.zeros((121,121,3), dtype=float)
    Vs2Polxyz = np.zeros((121,121,3), dtype=float)

    eigenVecs = 0
    eigenVals = 0
    for phi in range(0, len(Xx)):

        phi_index = phi
        for theta in range(0, len(Xx)):
            theta_index = theta

            X = np.array([(Xx[theta][phi]), (Yy[theta][phi]), (Zz[theta][phi])], dtype=float)
            for i in range(0,3):
                for k in range(0,3):
                    T[i,k] = 0.0
                    for j in range(0,3):
                        for l in range(0,3):

                            m = Voigt_ijkl[i][j]
                            n = Voigt_ijkl[k][l]
                            T[i,k] = T[i,k] + C[m-1,n-1] * X[j] * X[l]

            for i in range(0,3):
                for j in range(0,3):
                    TT[i,j] = 0.0
                    for k in range(0,3):
                        TT[i,j] = TT[i,j] + T[i,k] * T[j,k]

            

            eigenVals, eigenVecs = (np.linalg.eig(TT))
            idx = eigenVals.argsort()
            idy = eigenVals.argsort()
            eigenVals = eigenVals[idx]
            eigenVecs = eigenVecs[:,idy]

            vels2[theta,phi] = (math.sqrt((math.sqrt(eigenVals[0]))/rho))*10
            velp[theta,phi] = (math.sqrt((math.sqrt(eigenVals[2])) / rho))*10
            vels1[theta,phi] = (math.sqrt((math.sqrt(eigenVals[1])) / rho))*10

            VpPolxyz[theta,phi,:] = eigenVecs[:,2].transpose()
            Vs1Polxyz[theta,phi,:] = eigenVecs[:,1].transpose()
            Vs2Polxyz[theta,phi,:] = eigenVecs[:,0].transpose()



    Cij = C
    Sij = np.linalg.inv(Cij)
    Kvoigt = (1/9) * ((Cij[0,0] + Cij[1,1] + Cij[2,2]) + 2 * (Cij[0,1] + Cij[0,2] + Cij[1,2]))
    Kreuss =  1 / ((Sij[0,0] + Sij[1,1] + Sij[2,2]) + 2 *(Sij[0,1] + Sij[0,2] + Sij[1,2]))
    Kvrh = stats.mean([Kvoigt, Kreuss])

    Gvoigt = (1/15) * (Cij[0,0] + Cij[1,1] + Cij[2,2] + 3*(Cij[3,3] + Cij[4,4] + Cij[5,5]) - Cij[0,1] - Cij[0,2] - Cij[1,2])
    Greuss = 15 / (4 * (Sij[0,0] + Sij[1,1] + Sij[2,2] - Sij[0,1] - Sij[0,2] - Sij[1,2]) + 3 *(Sij[3,3]+Sij[4,4] + Sij[5,5]))
    Gvrh = stats.mean([Gvoigt, Greuss])

    VpIso = (((Kvrh+(4/3) * Gvrh) / (rho))**0.5) * 10
    #print(VpIso)
    VsIso = ((Gvrh/rho)**(0.5)) * 10
    VpVsIso = VpIso/VsIso

    ### time to begin plotting
    cmap1 = [[0.0, 1.0, 0.25, 0.0], [0.49,1.0,1.0, 0.0], [0.5, 0.5, 1.0,1.0],[0.51, 0.0, 1.0, 1.0],[1.0,0.5,0.0,1.0]]
    cmapA = [[0, 1, 1, 1],[0.0025, 1, 0.25, 0],[0.49, 1, 1,0],[0.5, 0.5, 1,1],[0.51, 0, 1,1],[1,0.5,0,1]]
    cmap = np.zeros((64,3), dtype=float)
    cmap2 = np.zeros((64,3), dtype=float)
    CMPTS = np.transpose(np.linspace(0.0,1.0, num=64, dtype=float))
    CMPTS2 = np.transpose(np.linspace(0.0,1.0, num=64, dtype=float))


    cmap1 = np.asarray(cmap1, dtype=float)
    cmapA = np.asarray(cmapA, dtype=float)
    f = 0



    for i in range(0,3):
        interpval = interpolate.interp1d(cmap1[:,0], cmap1[:,i+1])
        interpval2 = interpolate.interp1d(cmapA[:,0], cmapA[:,i+1])
        cmap[:,i] = interpval(CMPTS)
        cmap2[:,i] = interpval2(CMPTS2)


    NumQs = Res
    VpPol1 = np.zeros((41,41,3), dtype=float)
    VsPol1 = np.zeros((41,41,3), dtype=float)
    XxPol = np.zeros((41,41), dtype=float)
    YyPol = np.zeros((41,41), dtype=float)
    ZzPol1 = np.zeros((41,41), dtype=float)


    VpPol1 = VpPolxyz[0:len(VpPolxyz):NumQs,0:len(VpPolxyz):NumQs,:]
    Vs1Pol1 = Vs1Polxyz[0:len(Vs1Polxyz):NumQs,0:len(Vs1Polxyz):NumQs,:]
    XxPol = Xx[0:len(Xx):NumQs, 0:len(Xx):NumQs]
    YyPol = Yy[0:len(Yy):NumQs, 0:len(Yy):NumQs]
    ZzPol = Zz[0:len(Zz):NumQs, 0:len(Zz):NumQs]

    Vs1PolX = Vs1Pol1[:,:,0]
    Vs1PolY = Vs1Pol1[:,:,1]
    Vs1PolZ = Vs1Pol1[:,:,2]

    pxx = np.arange(pi*-1, pi, pi/180)
    pxx = np.append(pxx,pi)
    pxx = [math.cos(x) for x in pxx]



    py = np.arange(pi*-1, pi, pi/180)

    py = np.append(py,pi)
    py = [math.sin(y) for y in py]



    pxx = np.transpose(pxx)
    py = np.transpose(py)

    pxx = np.reshape(pxx, (len(pxx),1))
    py = np.reshape(py, (len(py),1))

    oX = np.reshape(np.transpose(np.array([-1.05, -1.05, 1.05, 1.05,-1.05])), (5, 1))
    oY = np.reshape(np.transpose(np.array([0.0, 1.05, 1.05,-1.05,-1.05])), (5,1))

    Pa = np.concatenate((pxx,py),axis=1)
    Pb = np.concatenate((oX, oY), axis=1)

    P1 = np.concatenate((Pa, Pb))


    XxReshapeSize = Xx.shape[0] * Xx.shape[1]
    YyReshapeSize = Yy.shape[0] * Yy.shape[1]
    ZzReshapeSize = Zz.shape[0] * Zz.shape[1]

    vpReshapeSize = velp.shape[0] * velp.shape[1]
    vels1ReshapeSize = vels1.shape[0] * vels1.shape[1]
    vels2ReshapeSize = vels2.shape[0] * vels2.shape[1]

    ### I WANT TO DIE
    Ux = np.zeros((XxReshapeSize,1), dtype=float)
    Ux = np.append(Ux, [Xx[:,x] for x in range(0,121)])
    Ux = np.delete(Ux, slice(0, XxReshapeSize))
    Ux = np.reshape(Ux, (XxReshapeSize,1))

    Uy = np.zeros((YyReshapeSize,1), dtype=float)
    Uy = np.append(Uy, [Yy[:,x] for x in range(0,121)])
    Uy = np.delete(Uy, slice(0, YyReshapeSize))
    Uy = np.reshape(Uy, (YyReshapeSize, 1))

    Uz = np.zeros((ZzReshapeSize,1), dtype=float)
    Uz = np.append(Uz, [Zz[:,x] for x in range(0,121)])
    Uz = np.delete(Uz, slice(0, ZzReshapeSize))
    Uz = np.reshape(Uz, (ZzReshapeSize,1))

    Vp = np.zeros((vpReshapeSize,1), dtype=float)
    Vp = np.append(Vp, [velp[:,x] for x in range(0,121)])
    Vp = np.delete(Vp, slice(0, vpReshapeSize))
    Vp = np.reshape(Vp, (vpReshapeSize, 1))

    Vs1 = np.zeros((vels1ReshapeSize,1), dtype=float)
    Vs1 = np.append(Vs1, [vels1[:,x] for x in range(0,121)])
    Vs1 = np.delete(Vs1, slice(0, vels1ReshapeSize))
    Vs1 = np.reshape(Vs1, (vels1ReshapeSize,1))

    Vs2 = np.zeros((vels2ReshapeSize,1), dtype=float)
    Vs2 = np.append(Vs2, [vels2[:,x] for x in range(0,121)])
    Vs2 = np.delete(Vs2, slice(0, vels2ReshapeSize))
    Vs2 = np.reshape(Vs2, (vels2ReshapeSize,1))


    Ux = Ux.flatten()
    Uy = Uy.flatten()
    Uz = Uz.flatten()

    Ux = Ux.tolist()
    Uy = Uy.tolist()
    Uz = Uz.tolist()


    Uxyz = np.array([(Ux), (Uy), (Uz)])

    Uxyz = np.transpose(Uxyz)

    df = pd.DataFrame(data=Uxyz)
    z1indicies = (df[df[2]<1].index.values)

    Uxyz = Uxyz[z1indicies,:]

    Vp = Vp[z1indicies]
    Vs1 = Vs1[z1indicies]
    Vs2 = Vs2[z1indicies]

    ### UXYZ is correct - 5/3/2023



    x, y = np.meshgrid(np.arange(-1,1.0,0.01), np.arange(-1,1,0.01))

    xaxis = np.multiply(XxPol, np.real(np.sqrt((1/(1-ZzPol)))))
    yaxis = np.multiply(YyPol, np.real(np.sqrt((1/(1-ZzPol)))))
    
    vels1X= 1 / vels1 #np.divide(1, vels1)
    vels2Y = 1/ vels2 #np.divide(1, vels2)
    
    quiver3DColorValues = (vels2Y-vels1X)
    bright_blue = [[0, '#7DF9FF'], [1, '#7DF9FF']]
    zero_pt = pd.Series([0])

    z1 = np.array([
    [1,1,1],
    [-1,-1,-1]
    ])
    y1 = [0, 0]
    x1 = [1,-1]

    XxPol = XxPol.flatten().tolist()
    YyPol = YyPol.flatten().tolist()
    ZzPol = ZzPol.flatten().tolist()
    Vs1PolX = Vs1PolX.flatten().tolist()
    Vs1PolY = Vs1PolY.flatten().tolist()
    Vs1PolZ = Vs1PolZ.flatten().tolist()

    xMeshInitiate, yMeshInitiate = np.arange(-1,1.01, 0.01), np.arange(-1,1.01, 0.01)
    ###2/13/2024 - store uxyz in field, access later on - split into two functions
    
    if plotType == 'Quiver' and plotType3D==None:
        firstInterSet = np.multiply(Uxyz[:,0], np.real(np.sqrt(np.divide(1, (1-Uxyz[:,2])))))

        secondInterSet = np.multiply(Uxyz[:,1], np.real(np.sqrt(np.divide(1, (1-Uxyz[:,2])))))
        
        thirdInterSet = np.divide(1, Vs2) - np.divide(1, Vs1)

        xTest = yTest = np.linspace(-1,1,len(firstInterSet))
        xTest, yTest = np.meshgrid(xMeshInitiate, yMeshInitiate)
    
        points = [firstInterSet.tolist(), secondInterSet.tolist()]
        values = [thirdInterSet.tolist()]
    
        data = interpolate.griddata((points[0], points[1]), values[0], (xTest, yTest), method='cubic')
    
        z_data = []
        for x in range(0, len(data)):
            combine = list(deepflatten(data[x], depth=1))
            z_data.append(combine)
        x_data = []
        combine = list(deepflatten(xTest, depth=1))
        x_data.append(combine)
        
    

        y_data = np.linspace(-1,1,200)
        xCos = []
        for x in range(0, 40000):
            xCos.append(math.cos(math.radians(x)))
        ySin = []
        for y in range(0, 200):
            ySin.append(math.sin(math.radians(y)))

        XxPolSquare =  XxPol * (np.real(np.sqrt(1 / (np.subtract(1, ZzPol)))))
        YyPolSquare = YyPol * (np.real(np.sqrt(1 / (np.subtract(1, ZzPol)))))


        fig = ff.create_quiver(XxPolSquare, YyPolSquare, Vs1PolX, Vs1PolY, arrow_scale=0.0001)
        #print(fig.data[0])
        fig.update_traces(name='Polarization Vectors (on/off)', visible='legendonly')
        fig.update_layout(legend=dict(
            yanchor="top",
            xanchor="left",
            y=0.99,
            x=0.01,
        ))
        fig.add_trace(go.Contour(x=x_data[0], y=y_data, z=z_data, line_smoothing=1, showscale=True,  contours_coloring='heatmap',
            colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        ))

        fig.update_layout(
            yaxis=dict(range=[-2,2]),
            xaxis=dict(range=[-3,3]), 
            showlegend=True,
            height=1000,
            margin=dict(l=20, r=20, t=20, b=20),
            modebar_add='drawclosedpath',
            newshape=dict(fillcolor='turquoise'),
            plot_bgcolor='rgba(0,0,0,0)',
            title={
                'text': 'VS1 Polarization and Splitting Time (s/km) <br> Max : {} Min {} '.format(np.around(np.nanmax(z_data), 4), np.around(np.nanmin(z_data), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
            )

        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[0],
            y=[0],
            marker_symbol='circle-open',
            showlegend=False, 
            marker_size=580,
            marker_line_width=162,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )
        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[-1.75],
            y=[0],
            marker_symbol='square',
            showlegend=False, 
            marker_size=300,
            marker_line_width=30,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )
        ### THESE NEXT TWO TRACES ARE TO COVER UP MARKS ON THE LEFT SIDE OF THE VS1 QUIVER/POLARIZATION PLOT - 1/16/2024

        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[-1.1],
            y=[1],
            marker_symbol='square',
            showlegend=False, 
            marker_size=50,
            marker_line_width=50,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )

        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[-1.1],
            y=[-1],
            marker_symbol='square',
            showlegend=False, 
            marker_size=50,
            marker_line_width=50,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )

        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[1.75],
            y=[0],
            marker_symbol='square',
            showlegend=False, 
            marker_size=300,
            marker_line_width=30,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )
        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)

        return fig


    if plotType == 'VP' and plotType3D==None:

        thirdInterSet = Vp
        xTest, yTest = np.meshgrid(xMeshInitiate, yMeshInitiate)

        fig = 1

        xx =  Uxyz[:,0] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))
        yy =  Uxyz[:,1] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))

        data = sp.interpolate.griddata((xx,yy), Vp, (xTest, yTest))
        x_data = []
        combine = list(deepflatten(xTest, depth=1))
        x_data.append(combine)

        y_data = np.linspace(-1,1,200)


        z_data = []
        for x in range(0, len(data)):
            combine = list(deepflatten(data[x], depth=1))
            z_data.append(combine)

        fig = go.Figure()
        fig.add_trace(go.Contour(x=x_data[0], y=y_data, z=z_data, line_smoothing=1, showscale=True, contours_coloring='heatmap',
            colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        ))
        
        fig.update_layout(
            yaxis=dict(range=[-2,2]),
            xaxis=dict(range=[-3,3]), 
            showlegend=True,
            height=1000,
            margin=dict(l=20, r=20, t=20, b=20),
            autosize=False,
            minreducedheight=750,
            minreducedwidth=1500,
            title={
                'text': 'VP <br> Max : {} Min {} '.format(np.around(np.nanmax(z_data), 4), np.around(np.nanmin(z_data), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
            )

        fig.update_layout(dragmode='drawrect',
                    # style of new shapes
                    newshape=dict(line_color='yellow',
                                fillcolor='turquoise',
                                opacity=0.5),
                    plot_bgcolor='rgba(0,0,0,0)')


        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[0],
            y=[0],
            marker_symbol='circle-open',
            showlegend=False, 
            marker_size=580,
            marker_line_width=162,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'

         )
        )

        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)

        return fig


    if plotType == 'VS1' and plotType3D==None:

        thirdInterSet = Vs1

        xTest, yTest = np.meshgrid(xMeshInitiate, yMeshInitiate)

        fig = 1

        xx =  Uxyz[:,0] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))
        yy =  Uxyz[:,1] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))

        data = sp.interpolate.griddata((xx,yy), Vs1, (xTest, yTest))
        x_data = []
        combine = list(deepflatten(xTest, depth=1))
        x_data.append(combine)
        y_data = np.linspace(-1,1,200)

        z_data = []
        for x in range(0, len(data)):
            combine = list(deepflatten(data[x], depth=1))
            z_data.append(combine)

        fig = go.Figure()
        fig.add_trace(go.Contour(x=x_data[0], y=y_data, z=z_data, line_smoothing=1,  showscale=True, contours_coloring='heatmap',
            colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        ))

        fig.update_layout(
            yaxis=dict(range=[-2,2]),
            xaxis=dict(range=[-3,3]), 
            showlegend=True,
            margin=dict(l=20, r=20, t=20, b=20),
            height=1000,
            plot_bgcolor='rgba(0,0,0,0)',
            title={
                'text': 'VS1 <br> Max : {} Min {} '.format(np.around(np.nanmax(z_data), 4), np.around(np.nanmin(z_data), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
            )


        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[0],
            y=[0],
            marker_symbol='circle-open',
            showlegend=False, 
            marker_size=580,
            marker_line_width=162,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )
        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)

        return fig
 

    if plotType == 'VS2' and plotType3D==None:

        thirdInterSet = Vs2

        xTest, yTest = np.meshgrid(xMeshInitiate, yMeshInitiate)

        fig = 1

        xx =  Uxyz[:,0] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))
        yy =  Uxyz[:,1] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))

        data = sp.interpolate.griddata((xx,yy), Vs2, (xTest, yTest))
        x_data = []
        combine = list(deepflatten(xTest, depth=1))
        x_data.append(combine)

        y_data = np.linspace(-1,1,200)


        z_data = []
        for x in range(0, len(data)):
            combine = list(deepflatten(data[x], depth=1))
            z_data.append(combine)

        
        #fig = ff.create_quiver(XxPol, YyPol, Vs1PolX, Vs1PolY)
        fig = go.Figure()
        fig.add_trace(go.Contour(x=x_data[0], y=y_data, z=z_data, line_smoothing=1,  showscale=True, contours_coloring='heatmap',
            colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        ))
        
        fig.update_layout(
            yaxis=dict(range=[-2,2]),
            xaxis=dict(range=[-3,3]),
            margin=dict(l=20, r=20, t=20, b=20), 
            showlegend=True,
            height=1000,
            plot_bgcolor='rgba(0,0,0,0)',
            title={
                'text': 'VS2 <br> Max : {} Min {} '.format(np.around(np.nanmax(z_data), 4), np.around(np.nanmin(z_data), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
            )

        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[0],
            y=[0],
            marker_symbol='circle-open',
            showlegend=False, 
            marker_size=580,
            marker_line_width=162,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )

        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)
   
        return fig

    if plotType == 'VPVS1' and plotType3D==None:

   
        thirdInterSet = Vp/Vs1

        xTest, yTest = np.meshgrid(xMeshInitiate, yMeshInitiate)

        fig = 1


        xx =  Uxyz[:,0] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))
        yy =  Uxyz[:,1] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))

        data = sp.interpolate.griddata((xx,yy), Vp/Vs1, (xTest, yTest))
        x_data = []
        combine = list(deepflatten(xTest, depth=1))
        x_data.append(combine)

        y_data = np.linspace(-1,1,200)


        z_data = []
        for x in range(0, len(data)):
            combine = list(deepflatten(data[x], depth=1))
            z_data.append(combine)

        
        #fig = ff.create_quiver(XxPol, YyPol, Vs1PolX, Vs1PolY)
        fig = go.Figure()
        fig.add_trace(go.Contour(x=x_data[0], y=y_data, z=z_data, line_smoothing=1, showscale=True, contours_coloring='heatmap',
            colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        ))
        

        fig.update_layout(
            yaxis=dict(range=[-2,2]),
            xaxis=dict(range=[-3,3]), 
            showlegend=True,
            height=1000,
            margin=dict(l=20, r=20, t=20, b=20),
            plot_bgcolor='rgba(0,0,0,0)',
            title={
                'text': 'VPVS1 <br> Max : {} Min {} '.format(np.around(np.nanmax(z_data), 4), np.around(np.nanmin(z_data), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
            )
        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[0],
            y=[0],
            marker_symbol='circle-open',
            showlegend=False, 
            marker_size=580,
            marker_line_width=162,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )

        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)
     
        return fig


    df = pd.DataFrame(list(zip(XxPol, YyPol, ZzPol, Vs1PolX, Vs1PolY, Vs1PolZ)), columns = ['x','y','z','u','v','w'])
    df = df.drop_duplicates()


    index_list = df.index


    if plotType3D=='3DQuiver' and plotType==None:
        #5/23/23 - x y z as origin point, u w v as end point, connect two 

        fig1 = make_subplots(rows=1, cols=1,
                            specs=[[{'is_3d': True}]],
                            subplot_titles=['Vs1 Polarization & splitting time (s/km)'],
                        )

        fig1.add_trace(go.Surface(x=Xx, y=Yy, z=Zz, surfacecolor=quiver3DColorValues, 
                colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]    
        ))

        Vs1PolXReversed = [x * -1 for x in Vs1PolX]
        Vs1PolYReversed = [x * -1 for x in Vs1PolY]
        Vs1PolZReversed = [x * -1 for x in Vs1PolZ]

        fig1.add_trace(go.Surface(z=z1, y=y1, x=x1, showscale=False))

        fig1.add_trace(go.Cone(
            x = XxPol,
            y = YyPol,
            z = ZzPol,
            u = Vs1PolX,
            v = Vs1PolY,
            w = Vs1PolZ,
            sizeref = 2.5,
            anchor= 'tail',
            showscale=False,
            colorscale="gray",
            cmax=2,
            cmin=1

        ))

        fig1.add_trace(go.Cone(
            x = XxPol,
            y = YyPol,
            z = ZzPol,
            u = Vs1PolXReversed,
            v = Vs1PolYReversed,
            w = Vs1PolZReversed,
            anchor = 'tail',
            sizeref = 2.5,
            showscale=False,
            colorscale="gray",
            cmax=2,
            cmin=1

        ))


        fig1.update_layout(
            scene=dict(
            xaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='X Axis')),
            yaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Y Axis')),
            zaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Z Axis')),
            camera=dict(eye=dict(x=-2, y=-2, z=1.25))
            ),
        )

        return fig1


        
    if plotType3D=='3DVP' and plotType==None:
        vpFig = make_subplots(rows=1, cols=1,
                        specs=[[{'is_3d': True}]],
                     )
        vpFig.add_trace(go.Surface(x=Xx, y=Yy, z=Zz, surfacecolor=velp,
                colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
         ), 1, 1)
        vpFig.add_trace(go.Surface(z=z1, y=y1, x=x1, showscale=False))

        vpFig.update_layout(
            scene=dict(
            xaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='X Axis')),
            yaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Y Axis')),
            zaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Z Axis')),
            camera=dict(eye=dict(x=-2, y=-2, z=1.25))
            ),
            title={
                'text': 'VP <br> Max : {} Min {} '.format(np.around(np.nanmax(velp), 4), np.around(np.nanmin(velp), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
        )

        return vpFig

    if plotType3D=='3DVS1' and plotType==None:
        vs1Fig = make_subplots(rows=1, cols=1,
                          specs=[[{'is_3d': True}]],
                     )
        vs1Fig.add_trace(go.Surface(x=Xx, y=Yy, z=Zz, surfacecolor=vels1,
                colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        
        ), 1, 1)
        vs1Fig.add_trace(go.Surface(z=z1, y=y1, x=x1, showscale=False))
        vs1Fig.update_layout(
            scene=dict(
            xaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='X Axis')),
            yaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Y Axis')),
            zaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Z Axis')),
            camera=dict(eye=dict(x=-2, y=-2, z=1.25))
            ),
            title={
                'text': 'VS1 <br> Max : {} Min {} '.format(np.around(np.nanmax(vels1), 4), np.around(np.nanmin(vels1), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
        )
        return vs1Fig

    if plotType3D=='3DVS2' and plotType==None:
        vs2Fig = make_subplots(rows=1, cols=1,
                            specs=[[{'is_3d': True}]],
                       )
        vs2Fig.add_trace(go.Surface(x=Xx, y=Yy, z=Zz, surfacecolor=vels2,
                colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        
        ), 1, 1)
        vs2Fig.add_trace(go.Surface(z=z1, y=y1, x=x1, showscale=False))
        vs2Fig.update_layout(
            scene=dict(
            xaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='X Axis')),
            yaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Y Axis')),
            zaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Z Axis')),
            camera=dict(eye=dict(x=-2, y=-2, z=1.25))
            ),
            title={
                'text': 'VS2 <br> Max : {} Min {} '.format(np.around(np.nanmax(vels2), 4), np.around(np.nanmin(vels2), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
        )
        return vs2Fig
    if plotType3D=='3DVPVS1' and plotType==None:
    
        vpvs1 = np.divide(velp, vels1)

        vpvs1Fig = make_subplots(rows=1, cols=1,
                        specs=[[{'is_3d': True}]],
                    )
        vpvs1Fig.add_trace(go.Surface(x=Xx, y=Yy, z=Zz, surfacecolor=vpvs1, 
                colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]), 1, 1)
                
        vpvs1Fig.add_trace(go.Surface(z=z1, y=y1, x=x1, showscale=False))
        vpvs1Fig.update_layout(
            scene=dict(
            xaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='X Axis')),
            yaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Y Axis')),
            zaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Z Axis')),
            camera=dict(eye=dict(x=-2, y=-2, z=1.25))
            ),
            title={
                'text': 'VPVS1 <br> Max : {} Min {} '.format(np.around(np.nanmax(vpvs1), 4), np.around(np.nanmin(vpvs1), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
        )

        return vpvs1Fig

    if plotType == 'BackAzimuthal' and plotType3D == None:
            
        VpPolX = VpPolxyz[:,:,0]
        VpPolY = VpPolxyz[:,:,1]
        VpPolZ = VpPolxyz[:,:,2]

        Vs1PolX = Vs1Polxyz[:,:,0]
        Vs1PolY = Vs1Polxyz[:,:,1]
        Vs1PolZ = Vs1Polxyz[:,:,2]

        Vs2PolX = Vs2Polxyz[:,:,0]
        Vs2PolY = Vs2Polxyz[:,:,1]
        Vs2PolZ = Vs2Polxyz[:,:,2]

        Ux = Xx.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
        Uy = Yy.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
        Uz = Zz.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
    
        Ux = Ux.flatten()
        Uy = Uy.flatten()
        Uz = Uz.flatten()

        Ux = Ux.tolist()
        Uy = Uy.tolist()
        Uz = Uz.tolist()  
        Uxyz = np.array([Ux, Uy, Uz]).transpose()

        VS1baz = vels1.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
        VS2baz = vels2.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
        VPbaz = velp.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()

        VpPolZbaz = VpPolZ.reshape(len(VpPolZ[0])*len(VpPolZ[1]), 1, order="F").copy()
        Vs1PolZbaz = Vs1PolZ.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
        Vs2PolZbaz = Vs2PolZ.reshape(len(Vs2PolZ[0])*len(Vs2PolZ[1]), 1, order="F").copy()

        VpPolXbaz = VpPolX.reshape(len(VpPolX[0])*len(VpPolX[1]), 1, order="F").copy()
        Vs1PolXbaz = Vs1PolX.reshape(len(Vs1PolX[0])*len(Vs1PolX[1]), 1, order="F").copy()
        Vs2PolXbaz = Vs2PolX.reshape(len(Vs2PolX[0])*len(Vs2PolX[1]), 1, order="F").copy()

        VpPolYbaz = VpPolY.reshape(len(VpPolY[0])*len(VpPolY[1]), 1, order="F").copy()
        Vs1PolYbaz = Vs1PolY.reshape(len(Vs1PolY[0])*len(Vs1PolY[1]), 1, order="F").copy()
        Vs2PolYbaz = Vs2PolY.reshape(len(Vs2PolY[0])*len(Vs2PolY[1]), 1, order="F").copy()

        HorzInt = 0.001
        Horzxyz = Uxyz[np.where(np.abs(Uxyz[:,2])<HorzInt),:]
        Horzxyz = Horzxyz[0,:]

        Lxy = np.sqrt(Horzxyz[:,0] **2+Horzxyz[:,1] **2 )


        VS1baz = VS1baz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        VS2baz = VS2baz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        VPbaz = VPbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]

        Vs1PolXbaz = Vs1PolXbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        Vs1PolYbaz = Vs1PolYbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        Vs1PolZbaz = Vs1PolZbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]

        Vs2PolXbaz = Vs2PolXbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        Vs2PolYbaz = Vs2PolYbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        Vs2PolZbaz = Vs2PolZbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]

        VpPolXbaz = VpPolXbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        VpPolYbaz = VpPolYbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        VpPolZbaz = VpPolZbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]

        BackAz = np.zeros((len(VPbaz), 1))


        for i in range(0, len(BackAz)):
            if Horzxyz[i, 1] >=0:
                BackAz[i] = np.rad2deg(np.arcsin(np.divide([Horzxyz[i,0]], Lxy[i])))
            elif Horzxyz[i,1] < 0:
                BackAz[i] = 180 - np.rad2deg(np.arcsin(np.divide([Horzxyz[i,0]], Lxy[i])))
        
        for x in range(0, len(BackAz)):
            if BackAz[x] < 0:
                BackAz[x] = BackAz[x] + 360
    

        Vs1Inc = np.absolute(np.rad2deg(np.arctan(Vs1PolZbaz / (Vs1PolXbaz ** 2 + Vs1PolYbaz ** 2)**0.5)))
        Vs2Inc = np.absolute(np.rad2deg(np.arctan(Vs2PolZbaz / (Vs2PolXbaz ** 2 + Vs2PolYbaz ** 2)**0.5)))
        Vs1Abs = []
        Vs2Abs = []

        randomColor = np.random.randn(500)
        ### 10/24/2023 - construct colorscale according to intensity - ranges are 0.0625
        rbColorScale = [
            [0, "rgb(255,0,0)"],
            [0.0625, "rgb(255,0,0)"],

            [0.0625, "rgb(255, 34, 34)"],
            [0.125, "rgb(255, 34, 34)"],

            [0.125, "rgb(255, 68, 68)"],
            [0.1875, "rgb(255, 68, 68)"],

            [0.1875, "rgb(255, 102, 102)"],
            [0.25, "rgb(255, 102, 102)"],

            [0.25, "rgb(255,136, 136)"],
            [0.3125, "rgb(255,136, 136)"],

            [0.3125, "rgb(255, 170, 170)"],
            [0.375, "rgb(255, 170, 170)"],

            [0.375, "rgb(255, 204, 204)"],
            [0.4375, "rgb(255, 204, 204)"],

            [0.4375, "rgb(255, 238, 238)"],
            [0.5, "rgb(255, 238, 238)"],

            [0.5, "rgb(238, 238, 255)"],
            [0.5625, "rgb(238, 238, 255)"],

            [0.5625, "rgb(204, 204, 255)"],
            [0.625, "rgb(204, 204, 255)"],
            
            [0.625, "rgb(170, 170, 255)"],
            [0.6875, "rgb(170, 170, 255)"],

            [0.6875, "rgb(136, 136, 255)"],
            [0.75, "rgb(136, 136, 255)"],

            [0.75, "rgb(102, 102, 255)"],
            [0.8125, "rgb(102, 102, 255)"],

            [0.8125, "rgb(68, 68, 255)"],
            [0.875, "rgb(68,68, 255)"],
            
            [0.875, "rgb(34,34,255)"],
            [0.9375, "rgb(34,34, 255)"],

            [0.9375, "rgb(0,0,255)"],
            [1, "rgb(0,0,255)"]
        ]

        for x in Vs2Inc:
            Vs2Abs.append(x[0])
        for x in Vs1Inc:
            Vs1Abs.append(x[0])


        fig = go.Figure()
        #10/4/2023 - iterate through, do every single point at once after trying to import colorscale manually
        fig.add_trace(go.Scatter(x=BackAz[:,0], y=VS2baz[:,0], mode='markers', marker=dict(
            colorscale=rbColorScale, color=Vs2Abs,  cmin=0, cmax=90, showscale=True
        ), name='Vs2 test'))

        fig.add_trace(go.Scatter(x=BackAz[:,0], y=VS1baz[:,0], mode='markers', marker=dict(
            colorscale=rbColorScale, color=Vs1Abs, cmin=0, cmax=90, showscale=True
            ), name='Vs1'))
        
        fig.update_layout(
            showlegend=False,
            height=800,
            width = 1500,
            plot_bgcolor='rgba(0,0,0,0)',
            title={
                'text': 'Vs with Backazimuth and Inclination from Horizontal',
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            },
            xaxis_title = 'Backazimuth',
            yaxis_title = 'Vs (km/s)'
            )




        return fig


    if plotType == 'RadialPlots':
        VpPolX = VpPolxyz[:,:,0]
        VpPolY = VpPolxyz[:,:,1]
        VpPolZ = VpPolxyz[:,:,2]

        Vs1PolX = Vs1Polxyz[:,:,0]
        Vs1PolY = Vs1Polxyz[:,:,1]
        Vs1PolZ = Vs1Polxyz[:,:,2]

        Vs2PolX = Vs2Polxyz[:,:,0]
        Vs2PolY = Vs2Polxyz[:,:,1]
        Vs2PolZ = Vs2Polxyz[:,:,2]

        Ux = Xx.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
        Uy = Yy.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
        Uz = Zz.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
    
        Ux = Ux.flatten()
        Uy = Uy.flatten()
        Uz = Uz.flatten()

        Ux = Ux.tolist()
        Uy = Uy.tolist()
        Uz = Uz.tolist()  
        Uxyz = np.array([Ux, Uy, Uz]).transpose()


        VS1baz = vels1.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
        VS2baz = vels2.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
        VPbaz = velp.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()

        VpPolZbaz = VpPolZ.reshape(len(VpPolZ[0])*len(VpPolZ[1]), 1, order="F").copy()
        Vs1PolZbaz = Vs1PolZ.reshape(len(Vs1PolZ[0])*len(Vs1PolZ[1]), 1, order="F").copy()
        Vs2PolZbaz = Vs2PolZ.reshape(len(Vs2PolZ[0])*len(Vs2PolZ[1]), 1, order="F").copy()

        VpPolXbaz = VpPolX.reshape(len(VpPolX[0])*len(VpPolX[1]), 1, order="F").copy()
        Vs1PolXbaz = Vs1PolX.reshape(len(Vs1PolX[0])*len(Vs1PolX[1]), 1, order="F").copy()
        Vs2PolXbaz = Vs2PolX.reshape(len(Vs2PolX[0])*len(Vs2PolX[1]), 1, order="F").copy()

        VpPolYbaz = VpPolY.reshape(len(VpPolY[0])*len(VpPolY[1]), 1, order="F").copy()
        Vs1PolYbaz = Vs1PolY.reshape(len(Vs1PolY[0])*len(Vs1PolY[1]), 1, order="F").copy()
        Vs2PolYbaz = Vs2PolY.reshape(len(Vs2PolY[0])*len(Vs2PolY[1]), 1, order="F").copy()

        HorzInt = 0.001
        Horzxyz = Uxyz[np.where(np.abs(Uxyz[:,2])<HorzInt),:]
        Horzxyz = Horzxyz[0,:]

        Lxy = np.sqrt(Horzxyz[:,0] **2+Horzxyz[:,1] **2 )


        VS1baz = VS1baz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        VS2baz = VS2baz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        VPbaz = VPbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
   

        Vs1PolXbaz = Vs1PolXbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        Vs1PolYbaz = Vs1PolYbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        Vs1PolZbaz = Vs1PolZbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]


        Vs2PolXbaz = Vs2PolXbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        Vs2PolYbaz = Vs2PolYbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]
        Vs2PolZbaz = Vs2PolZbaz[np.where(np.abs(Uxyz[:,2])<HorzInt)]

        VpPolXbaz = VpPolXbaz[np.where(np.abs(Uxyz[:,2]) < HorzInt)]
        VpPolYbaz = VpPolYbaz[np.where(np.abs(Uxyz[:,2]) < HorzInt)]
        VpPolZbaz = VpPolZbaz[np.where(np.abs(Uxyz[:,2]) < HorzInt)]

        BackAz = np.zeros((len(VPbaz), 1))


        for i in range(0, len(BackAz)):
            if Horzxyz[i, 1] >=0:
                BackAz[i] = np.rad2deg(np.arcsin(np.divide([Horzxyz[i,0]], Lxy[i])))
            elif Horzxyz[i,1] < 0:
                BackAz[i] = 180 - np.rad2deg(np.arcsin(np.divide([Horzxyz[i,0]], Lxy[i])))
        
        for x in range(0, len(BackAz)):
            if BackAz[x] < 0:
                BackAz[x] = BackAz[x] + 360
    

        Vs1Inc = np.absolute(np.rad2deg(np.arctan(Vs1PolZbaz / (Vs1PolXbaz ** 2 + Vs1PolYbaz ** 2)**0.5)))
        Vs2Inc = np.absolute(np.rad2deg(np.arctan(Vs2PolZbaz / (Vs2PolXbaz ** 2 + Vs2PolYbaz ** 2)**0.5)))

        R_PV1  = np.zeros((len(BackAz), 3))
        R_PV2  = np.zeros((len(BackAz), 3))
        R_PVP  = np.zeros((len(BackAz), 3))



        for i in range (0, len(BackAz)):
            cosdOne = np.cos((np.deg2rad(BackAz[i])))
            sindOne = np.sin((np.deg2rad(BackAz[i]))*-1)
            sindTwo = np.sin((np.deg2rad(BackAz[i])))
            cosdTwo = np.cos((np.deg2rad(BackAz[i])))
            R = [[cosdOne[0], sindOne[0], 0], [sindTwo[0], cosdTwo[0], 0], [0,0,1]]
            R = np.asarray(R)
            

            PV1 = np.array((Vs1PolXbaz[i],Vs1PolYbaz[i],Vs1PolZbaz[i]))

            PV2 = np.array((Vs2PolXbaz[i],Vs2PolYbaz[i],Vs2PolZbaz[i]))

            PVP = np.array((VpPolXbaz[i],VpPolYbaz[i],VpPolZbaz[i]))
            

            R_PV1[i,:] = np.transpose(np.matmul(R, PV1))
            R_PV2[i,:] = np.transpose(np.matmul(R, PV2))
            R_PVP[i,:] = np.transpose(np.matmul(R, PVP))
        VS1T = R_PV1[:,0]
        VS1R = R_PV1[:,1]
        VS1V = R_PV1[:,2]

        VS2T = R_PV2[:,0]
        VS2R = R_PV2[:,1]
        VS2V = R_PV2[:,2]

        VPT = R_PVP[:,0]
        VPR = R_PVP[:,1]
        VPV = R_PVP[:,2]

        BAZsort = np.sort(BackAz, axis=0)
        BAZargsort = np.argsort(BackAz, axis=0)


        VS1RPlot = np.zeros((121, 1))
        VS1VPlot = np.zeros((121, 1))
        VS1TPlot = np.zeros((121, 1))
        VS2RPlot = np.zeros((121, 1))
        VS2TPlot = np.zeros((121, 1))
        VS2VPlot = np.zeros((121, 1))
        #for x in BAZargsort:
            #y = VS1T[x[0]]
           #VS1TPlot[:,0] = y
        for x in range (0, len(BAZargsort)):
            y = BAZargsort[x]
            VS1Rholder = VS1R[y[0]]
            VS1Tholder = VS1T[y[0]]
            VS1Vholder = VS1V[y[0]]
            VS2Rholder = VS2R[y[0]]
            VS2Tholder = VS2T[y[0]]
            VS2Vholder = VS2V[y[0]]

            VS1RPlot[x,0] = np.abs(VS1Rholder)
            VS1TPlot[x,0] = np.abs(VS1Tholder)
            VS1VPlot[x,0] = np.abs(VS1Vholder)
            VS2RPlot[x,0] = np.abs(VS2Rholder)
            VS2TPlot[x,0] = np.abs(VS2Tholder)
            VS2VPlot[x,0] = np.abs(VS2Vholder)

        fig = go.Figure()
        #fig.add_trace(go.Scatter(x=BackAz[:,0], y=VS2baz[:,0], mode='markers', name='filled'))
        #print(BAZsort[:,0])
        fig.add_trace(go.Scatter(x=BAZsort[:,0], y=VS1RPlot[:,0], marker=dict(color="#222A2A"), mode='lines', name='VS1R'))
        fig.add_trace(go.Scatter(x=BAZsort[:,0], y=VS1TPlot[:,0], marker=dict(color="#1CBE4F"), mode='lines', name='VS1T'))
        fig.add_trace(go.Scatter(x=BAZsort[:,0], y=VS1VPlot[:,0], marker=dict(color="#3283FE"), mode='lines', name='VS1V'))
        fig.add_trace(go.Scatter(x=BAZsort[:,0], y=VS2RPlot[:,0], marker=dict(color="#222A2A"), mode='markers', name='VS2R'))
        fig.add_trace(go.Scatter(x=BAZsort[:,0], y=VS2TPlot[:,0], marker=dict(color="#1CBE4F"), mode='markers', name='VS2T'))
        fig.add_trace(go.Scatter(x=BAZsort[:,0], y=VS2VPlot[:,0], marker=dict(color="#3283FE"), mode='markers', name='VS2V'))

        #fig.add_trace(go.Scatter(x=mode='markers', name='filled'))




        fig.update_layout(
            showlegend=True,
            height=800,
            width = 1500,
            plot_bgcolor='rgba(0,0,0,0)',
            title={
                'text': 'Vs with Backazimuth and Inclination from Horizontal',
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            },
            xaxis_title = 'Backazimuth',
            yaxis_title = 'Vs (km/s)'
            )

        return fig


def tv_fold_model(tensor_index, user_input=False, plotType=None, plotType3D=None, userInputDataFrame=None):

    christoffelParameters = ChristoffelDevParams.ChristoffelParameters()

    pi = math.pi
    user_input = user_input

    aggTensorColumns = christoffelParameters.columnNamesTensorCSV
    sphereColumnNames = [x for x in range(0, 121)]
    csv_tensor_number = [x for x in range(0, 50)]
    user_tensor_number = [x for x in range(97, 147)]
    Res = 100

    user_to_csv_mapping = zip(user_tensor_number, csv_tensor_number)
    user_to_csv_mapping = dict(user_to_csv_mapping)
    returnedUserInputDataFrame = returnUserInputDataFrame(userInputDataFrame)
    


    aggTensorCSVLocation = 'databases/DB_Tensors_Rhos_OrigOrientation_CSV.csv'

    C = []

    tensorDF = pd.read_csv(aggTensorCSVLocation, sep=',', engine='python')


    tensorIndex = int(tensor_index)
    print(user_input)
    if user_input == True:
        tensorDF = returnedUserInputDataFrame
        tensorList = (tensorDF.iloc[tensorIndex, 4:25])
        rho = (tensorDF.iloc[tensorIndex, 3])
    else:
        tensorList = (tensorDF.iloc[tensorIndex-1, 4:25])
        rho = (tensorDF.iloc[tensorIndex-1, 3])
    print(tensorList)
    rho = float(rho)
    
    Tensors = [[tensorList[0],tensorList[1],tensorList[2],tensorList[3],tensorList[4],tensorList[5]],
             [tensorList[1], tensorList[6], tensorList[7], tensorList[8],tensorList[9],tensorList[10]],
             [tensorList[2],tensorList[7],tensorList[11],tensorList[12],tensorList[13],tensorList[14]],
             [tensorList[3],tensorList[8],tensorList[12],tensorList[15],tensorList[16],tensorList[17]],
             [tensorList[4],tensorList[9],tensorList[13],tensorList[16],tensorList[18],tensorList[19]],
             [tensorList[5],tensorList[10],tensorList[14],tensorList[17],tensorList[19],tensorList[20]]
             ]

    C = [[tensorList[0], tensorList[1], tensorList[2], tensorList[3], tensorList[4], tensorList[5]],
        [tensorList[1], tensorList[6], tensorList[7], tensorList[8], tensorList[9], tensorList[10]],
        [tensorList[2], tensorList[7], tensorList[11], tensorList[12], tensorList[13], tensorList[14]],
        [tensorList[3], tensorList[8], tensorList[12], tensorList[15], tensorList[16], tensorList[17]],
        [tensorList[4], tensorList[9], tensorList[13], tensorList[16], tensorList[18], tensorList[19]],
        [tensorList[5], tensorList[10], tensorList[14], tensorList[17], tensorList[19], tensorList[20]]
    ]


    Ang = []
    for ang in range(-180,181):
        Ang.append(ang)
    Amp = 25
    Per = 4
    ResN= 100
    Shape = [] 
    for x in Ang:
        Shape.append(Amp*math.sin(math.radians(Per*x)))
    
    TLoc = np.linspace(1, len(Ang)-1 ,ResN)
    TLoc = np.round(TLoc)


    diff = np.diff(Shape)
    dip_staging = -1*np.degrees(np.arctan(diff))
    dip = []

    for x in TLoc:
        dip.append(dip_staging[int(x)-1])
    RotatedCij = np.zeros((100,22))

    rake = np.zeros((100,1))
    
    trend = np.zeros((100,1))
    SL = 0
    L1 = 0
    L2 = 0 

    
    strike = np.zeros((100,1))
    for x in range(0, len(trend)):
        trend[x] = 90
    
    for i in range(0, len(dip)):
        if (trend[i] - strike[i]) == 180:
            rake[i] = 180
        elif trend[i] < strike[i]:
            trend[i] = trend[i] + 360
        if dip[i] == 0:
            rake[i] = trend[i]
        else:
            SL = np.cos(np.deg2rad(trend[i]-strike[i]))
            L1 = np.sin(np.deg2rad(trend[i]-strike[i]))
            L2 = L1 / np.cos(np.deg2rad(dip[i]))
            rake[i] = np.degrees(np.arctan(L2/SL))
            if rake[i] < 0:
                rake[i] = 180 + rake[i]
        rake[i] = -1*rake[i]

    ya = []
    xa = []
    za = []
    yrake = []

    for i in range(0, len(dip)):
        ya.append(dip[i])
        xa.append(0)
        za.append(strike[i][0])
        yrake.append(rake[i][0]-90)

    for s in range(0, len(dip)):
        zr = float(-1*za[s]*((math.pi)/180))
        xr = float(-1*xa[s]*((math.pi)/180))
        yr = float(ya[s]*((math.pi)/180))
        yrk = float((1*yrake[s])*((math.pi)/180))

        Rx = [[1,0,0,0,0,0], 
        [0, np.cos(xr)**2, np.sin(xr)**2, 2*np.cos(xr)*np.sin(xr), 0, 0], 
        [0, np.sin(xr)**2, np.cos(xr)**2, -2*np.sin(xr)*np.cos(xr), 0, 0], 
        [0, -np.cos(xr)*np.sin(xr), np.cos(xr)*np.sin(xr), (np.cos(xr)**2)-(np.sin(xr)**2),0,0],
        [0, 0, 0, 0, np.cos(xr), -np.sin(xr)],
        [0, 0, 0, 0, np.sin(xr), np.cos(xr)]
        ]



        Ry = [[np.cos(yr)**2, 0, np.sin(yr)**2, 0, 2*np.cos(yr)*np.sin(yr),0], 
        [0,1,0,0,0,0], 
        [np.sin(yr)**2, 0, np.cos(yr)**2, 0, -2*np.cos(yr)*np.sin(yr),0], 
        [0,0,0, np.cos(yr), 0, -1*np.sin(yr)],
        [-1*np.cos(yr)*np.sin(yr), 0, np.cos(yr)*np.sin(yr),0, (np.cos(yr)**2)-(np.sin(yr)**2),0],
        [0,0,0, np.sin(yr), 0, np.cos(yr)]
        ]


        Rz = [[np.cos(zr)**2, np.sin(zr)**2, 0,0,0, -2*np.cos(zr)*np.sin(zr)], 
        [np.sin(zr)**2, np.cos(zr)**2, 0,0,0, 2*np.cos(zr)*np.sin(zr)], 
        [0,0,1,0,0,0], 
        [0,0,0, np.cos(zr),np.sin(zr),0],
        [0,0,0,-1*np.sin(zr), np.cos(zr),0],
        [np.cos(zr)*np.sin(zr), -1*np.cos(zr)*np.sin(zr),0,0,0, (np.cos(zr)**2)-(np.sin(zr)**2)]
        ]
        
        
        Ryrake = [[np.cos(yrk)**2, np.sin(yrk)**2, 0,0,0, -2*(np.cos(yrk))*(np.sin(yrk))], 
        [np.sin(yrk)**2, np.cos(yrk)**2, 0, 0, 0, 2*np.cos(yrk)*np.sin(yrk)], 
        [0,0,1,0,0,0], 
        [0,0,0,np.cos(yrk),np.sin(yrk),0],
        [0,0,0, -1*np.sin(yrk), np.cos(yrk),0],
        [np.cos(yrk)*np.sin(yrk), -1*np.cos(yrk)*np.sin(yrk),0,0,0, (np.cos(yrk)**2)-(np.sin(yrk)**2)]
        ]
        


        CR = C
        CR = np.matmul(np.matmul(Ryrake, CR),np.transpose(Ryrake))
        CR = np.matmul(np.matmul(Rx, CR),np.transpose(Rx))
        CR = np.matmul(np.matmul(Ry, CR),np.transpose(Ry))
        CR = np.matmul(np.matmul(Rz, CR),np.transpose(Rz))

        RotatedCij[s][0] = rho
        RotatedCij[s][1:7] = CR[0]
        RotatedCij[s][7:12] = CR[1][1:6]
        RotatedCij[s][12:16] = CR[2][2:6]
        RotatedCij[s][16:19] = CR[3][3:6]
        RotatedCij[s][19:21] = CR[4][4:6]
        RotatedCij[s][21:23] = CR[5][5:7]

 
    tensorsV = []
    tensorsR = []
    #print(RotatedCij[0][1:])

    V_ave = np.zeros((36,1), dtype=float)
    R_ave = np.zeros((36,1), dtype=float)

    for array in range(0,100):
        tensorsV.append(np.zeros((6,6),dtype=float))
        tensorsR.append(np.zeros((6,6),dtype=float))
    for i in range(0,100):
        C = RotatedCij[i][1:]
        #all indicies below are subtracted by 1 for indexing purposes
        
        # second row had been C[6], C[7],C[8],C[9], C[10],C[11]
        tensorsV[i] = [[C[0], C[1], C[2], C[3], C[4], C[5]], 
        [C[1], C[6],C[7],C[8], C[9],C[10]], 
        [C[2], C[7], C[11], C[12], C[13], C[14]], 
        [C[3], C[8], C[12], C[15], C[16],C[17]],
        [C[4], C[9], C[13], C[16], C[18],C[19]],
        [C[5], C[10], C[14], C[17], C[19], C[20]]
        ]
        tensorsV[i] = np.asarray(tensorsV[i])
        
    for i in range(0,100):
        tensorsR[i] = np.linalg.pinv(tensorsV[i])

    for i in range(0,36):
        PV = np.zeros((100,1), dtype=float)
        PR = np.zeros((100,1), dtype=float)
        for j in range(0,100):
            flattenedVTensor = tensorsV[j].flatten()
            flattenedRTensor = tensorsR[j].flatten()
            PV[j,] = flattenedVTensor[i]
            PR[j,] = flattenedRTensor[i]
        V_ave[i] = np.mean(PV)
        R_ave[i] = np.mean(PR)

    


    
    test_array = np.zeros((6,6), dtype=float)

    #re-do average arrays
    
    V_ave = np.reshape(V_ave, (6,6))
    R_ave = np.reshape(R_ave, (6,6))

    #V_ave = V_ave.flatten()
    #V_ave = V_ave.tolist()
    
    #V_ave_rounded = []
    #R_ave_rounded = []

    #for x in V_ave:
        #y = round(x,4)
        #V_ave_rounded.append(y)
    VRH = (V_ave + linalg.inv(R_ave))/2

    C = np.reshape(VRH, (6,6))
    

    Res = 3

    voigtIndices = [[1,6,5], [6,2,4], [5,4,3]]

    Voigt_ijkl = pd.DataFrame(data=voigtIndices)
    Voigt_ijkl = Voigt_ijkl.to_numpy()

    rho = rho

    Xx = pd.read_csv('databases/xx.csv', names=sphereColumnNames, sep=',', engine='python')
    Yy = pd.read_csv('databases/yy.csv', names=sphereColumnNames, sep=',', engine='python')
    Zz = pd.read_csv('databases/zz.csv', names=sphereColumnNames, sep=',', engine='python')

    Xx = Xx.fillna(0)
    Yy = Yy.fillna(0)
    Zz = Zz.fillna(0)


    Xx = Xx.to_numpy()
    Yy = Yy.to_numpy()
    Zz = Zz.to_numpy()


    test1 = 0
    test2 = 0

    m = 0
    n = 0

    T = np.zeros((3,3), dtype=float)
    TT = np.zeros((3,3), dtype=float)
    velp = np.zeros((121,121), dtype=float)
    vels1 = np.zeros((121,121), dtype=float)
    vels2 = np.zeros((121,121), dtype=float)

    VpPolxyz = np.zeros((121,121,3), dtype=float)
    Vs1Polxyz = np.zeros((121,121,3), dtype=float)
    Vs2Polxyz = np.zeros((121,121,3), dtype=float)

    eigenVecs = 0
    eigenVals = 0
    for phi in range(0, len(Xx)):

        phi_index = phi
        for theta in range(0, len(Xx)):
            theta_index = theta

            X = np.array([(Xx[theta][phi]), (Yy[theta][phi]), (Zz[theta][phi])], dtype=float)
            for i in range(0,3):
                for k in range(0,3):
                    T[i,k] = 0.0
                    for j in range(0,3):
                        for l in range(0,3):

                            m = Voigt_ijkl[i][j]
                            n = Voigt_ijkl[k][l]
                            T[i,k] = T[i,k] + C[m-1,n-1] * (X[j] * X[l])

            for i in range(0,3):
                for j in range(0,3):
                    TT[i,j] = 0.0
                    for k in range(0,3):
                        TT[i,j] = TT[i,j] + T[i,k] * T[j,k]

            eigenVals, eigenVecs = (np.linalg.eig(TT))

            eigenVals, eigenVecs = (np.linalg.eig(TT))
            idx = eigenVals.argsort()
            idy = eigenVals.argsort()
            eigenVals = eigenVals[idx]
            eigenVecs = eigenVecs[:,idy]

            vels2[theta,phi] = (math.sqrt((math.sqrt(eigenVals[0]))/rho))*10
            velp[theta,phi] = (math.sqrt((math.sqrt(eigenVals[2])) / rho))*10
            vels1[theta,phi] = (math.sqrt((math.sqrt(eigenVals[1])) / rho))*10

            VpPolxyz[theta,phi,:] = eigenVecs[:,2].transpose()
            Vs1Polxyz[theta,phi,:] = eigenVecs[:,1].transpose()
            Vs2Polxyz[theta,phi,:] = eigenVecs[:,0].transpose()


    Cij = C
    Sij = np.linalg.inv(Cij)
    Kvoigt = (1/9) * ((Cij[0,0] + Cij[1,1] + Cij[2,2]) + 2 * (Cij[0,1] + Cij[0,2] + Cij[1,2]))
    Kreuss =  1 / ((Sij[0,0] + Sij[1,1] + Sij[2,2]) + 2 *(Sij[0,1] + Sij[0,2] + Sij[1,2]))
    Kvrh = stats.mean([Kvoigt, Kreuss])

    Gvoigt = (1/15) * (Cij[0,0] + Cij[1,1] + Cij[2,2] + 3*(Cij[3,3] + Cij[4,4] + Cij[5,5]) - Cij[0,1] - Cij[0,2] - Cij[1,2])
    Greuss = 15 / (4 * (Sij[0,0] + Sij[1,1] + Sij[2,2] - Sij[0,1] - Sij[0,2] - Sij[1,2]) + 3 *(Sij[3,3]+Sij[4,4] + Sij[5,5]))
    Gvrh = stats.mean([Gvoigt, Greuss])

    VpIso = (((Kvrh+(4/3) * Gvrh) / (rho))**0.5) * 10
    #print(VpIso)
    VsIso = ((Gvrh/rho)**(0.5)) * 10
    VpVsIso = VpIso/VsIso

    ### time to begin plotting
    cmap1 = [[0.0, 1.0, 0.25, 0.0], [0.49,1.0,1.0, 0.0], [0.5, 0.5, 1.0,1.0],[0.51, 0.0, 1.0, 1.0],[1.0,0.5,0.0,1.0]]
    cmapA = [[0, 1, 1, 1],[0.0025, 1, 0.25, 0],[0.49, 1, 1,0],[0.5, 0.5, 1,1],[0.51, 0, 1,1],[1,0.5,0,1]]
    cmap = np.zeros((64,3), dtype=float)
    cmap2 = np.zeros((64,3), dtype=float)
    CMPTS = np.transpose(np.linspace(0.0,1.0, num=64, dtype=float))
    CMPTS2 = np.transpose(np.linspace(0.0,1.0, num=64, dtype=float))


    cmap1 = np.asarray(cmap1, dtype=float)
    cmapA = np.asarray(cmapA, dtype=float)
    f = 0



    for i in range(0,3):
        interpval = interpolate.interp1d(cmap1[:,0], cmap1[:,i+1])
        interpval2 = interpolate.interp1d(cmapA[:,0], cmapA[:,i+1])
        cmap[:,i] = interpval(CMPTS)
        cmap2[:,i] = interpval2(CMPTS2)


    NumQs = Res
    VpPol1 = np.zeros((41,41,3), dtype=float)
    VsPol1 = np.zeros((41,41,3), dtype=float)
    XxPol = np.zeros((41,41), dtype=float)
    YyPol = np.zeros((41,41), dtype=float)
    ZzPol1 = np.zeros((41,41), dtype=float)


    VpPol1 = VpPolxyz[0:len(VpPolxyz):NumQs,0:len(VpPolxyz):NumQs,:]
    Vs1Pol1 = Vs1Polxyz[0:len(Vs1Polxyz):NumQs,0:len(Vs1Polxyz):NumQs,:]
    XxPol = Xx[0:len(Xx):NumQs, 0:len(Xx):NumQs]
    YyPol = Yy[0:len(Yy):NumQs, 0:len(Yy):NumQs]
    ZzPol = Zz[0:len(Zz):NumQs, 0:len(Zz):NumQs]

    Vs1PolX = Vs1Pol1[:,:,0]
    Vs1PolY = Vs1Pol1[:,:,1]
    Vs1PolZ = Vs1Pol1[:,:,2]

    px = np.arange(pi*-1, pi, pi/180)
    px = np.append(px,pi)
    px = [math.cos(x) for x in px]



    py = np.arange(pi*-1, pi, pi/180)

    py = np.append(py,pi)
    py = [math.sin(y) for y in py]



    px = np.transpose(px)
    py = np.transpose(py)

    px = np.reshape(px, (len(px),1))
    py = np.reshape(py, (len(py),1))

    oX = np.reshape(np.transpose(np.array([-1.05, -1.05, 1.05, 1.05,-1.05])), (5, 1))
    oY = np.reshape(np.transpose(np.array([0.0, 1.05, 1.05,-1.05,-1.05])), (5,1))

    Pa = np.concatenate((px,py),axis=1)
    Pb = np.concatenate((oX, oY), axis=1)

    P1 = np.concatenate((Pa, Pb))


    XxReshapeSize = Xx.shape[0] * Xx.shape[1]
    YyReshapeSize = Yy.shape[0] * Yy.shape[1]
    ZzReshapeSize = Zz.shape[0] * Zz.shape[1]

    vpReshapeSize = velp.shape[0] * velp.shape[1]
    vels1ReshapeSize = vels1.shape[0] * vels1.shape[1]
    vels2ReshapeSize = vels2.shape[0] * vels2.shape[1]

    ### I WANT TO DIE
    Ux = np.zeros((XxReshapeSize,1), dtype=float)
    Ux = np.append(Ux, [Xx[:,x] for x in range(0,121)])
    Ux = np.delete(Ux, slice(0, XxReshapeSize))
    Ux = np.reshape(Ux, (XxReshapeSize,1))

    Uy = np.zeros((YyReshapeSize,1), dtype=float)
    Uy = np.append(Uy, [Yy[:,x] for x in range(0,121)])
    Uy = np.delete(Uy, slice(0, YyReshapeSize))
    Uy = np.reshape(Uy, (YyReshapeSize, 1))

    Uz = np.zeros((ZzReshapeSize,1), dtype=float)
    Uz = np.append(Uz, [Zz[:,x] for x in range(0,121)])
    Uz = np.delete(Uz, slice(0, ZzReshapeSize))
    Uz = np.reshape(Uz, (ZzReshapeSize,1))

    Vp = np.zeros((vpReshapeSize,1), dtype=float)
    Vp = np.append(Vp, [velp[:,x] for x in range(0,121)])
    Vp = np.delete(Vp, slice(0, vpReshapeSize))
    Vp = np.reshape(Vp, (vpReshapeSize, 1))

    Vs1 = np.zeros((vels1ReshapeSize,1), dtype=float)
    Vs1 = np.append(Vs1, [vels1[:,x] for x in range(0,121)])
    Vs1 = np.delete(Vs1, slice(0, vels1ReshapeSize))
    Vs1 = np.reshape(Vs1, (vels1ReshapeSize,1))

    Vs2 = np.zeros((vels2ReshapeSize,1), dtype=float)
    Vs2 = np.append(Vs2, [vels2[:,x] for x in range(0,121)])
    Vs2 = np.delete(Vs2, slice(0, vels2ReshapeSize))
    Vs2 = np.reshape(Vs2, (vels2ReshapeSize,1))


    Ux = Ux.flatten()
    Uy = Uy.flatten()
    Uz = Uz.flatten()

    Ux = Ux.tolist()
    Uy = Uy.tolist()
    Uz = Uz.tolist()


    Uxyz = np.array([(Ux), (Uy), (Uz)])

    Uxyz = np.transpose(Uxyz)

    df = pd.DataFrame(data=Uxyz)
    z1indicies = (df[df[2]<1].index.values)
    Uxyz = Uxyz[z1indicies,:]
    Vp = Vp[z1indicies]
    Vs1 = Vs1[z1indicies]
    Vs2 = Vs2[z1indicies]

    ### UXYZ is correct - 5/3/2023



    x, y = np.meshgrid(np.arange(-1,1.0,0.01), np.arange(-1,1,0.01))

    xaxis = np.multiply(XxPol, np.real(np.sqrt((1/(1-ZzPol)))))
    yaxis = np.multiply(YyPol, np.real(np.sqrt((1/(1-ZzPol)))))

    vels1X= 1 / vels1
    vels2Y = 1 / vels2
    
    quiver3DColorValues = (vels2Y-vels1X)
    bright_blue = [[0, '#7DF9FF'], [1, '#7DF9FF']]
    zero_pt = pd.Series([0])

    z1 = np.array([
    [1,1,1],
    [-1,-1,-1]
    ])
    y1 = [0, 0]
    x1 = [1,-1]

    XxPol = XxPol.flatten().tolist()
    YyPol = YyPol.flatten().tolist()
    ZzPol = ZzPol.flatten().tolist()
    Vs1PolX = Vs1PolX.flatten().tolist()
    Vs1PolY = Vs1PolY.flatten().tolist()
    Vs1PolZ = Vs1PolZ.flatten().tolist()

    xMeshInitiate, yMeshInitiate = np.arange(-1,1.01, 0.01), np.arange(-1,1.01, 0.01)

    
    if plotType == 'Quiver' and plotType3D==None:
        print('in quiver in fold')
        firstInterSet = np.multiply(Uxyz[:,0], np.real(np.sqrt(np.divide(1, (1-Uxyz[:,2])))))

        secondInterSet = np.multiply(Uxyz[:,1], np.real(np.sqrt(np.divide(1, (1-Uxyz[:,2])))))
        
        thirdInterSet = np.divide(1, Vs2) - np.divide(1, Vs1)

        xTest = yTest = np.linspace(-1,1,len(firstInterSet))
        xTest, yTest = np.meshgrid(xMeshInitiate, yMeshInitiate)


    
        points = [firstInterSet.tolist(), secondInterSet.tolist()]
        values = [thirdInterSet.tolist()]
    
        data = interpolate.griddata((points[0], points[1]), values[0], (xTest, yTest), method='cubic')
    
        z_data = []
        for x in range(0, len(data)):
            combine = list(deepflatten(data[x], depth=1))
            z_data.append(combine)
        x_data = []
        combine = list(deepflatten(xTest, depth=1))
        x_data.append(combine)
        
    

        y_data = np.linspace(-1,1,200)
        xCos = []
        for x in range(0, 40000):
            xCos.append(math.cos(math.radians(x)))
        ySin = []
        for y in range(0, 200):
            ySin.append(math.sin(math.radians(y)))
 
        XxPolSquare =  XxPol * (np.real(np.sqrt(1 / (np.subtract(1, ZzPol)))))
        YyPolSquare = YyPol * (np.real(np.sqrt(1 / (np.subtract(1, ZzPol)))))


        fig = ff.create_quiver(XxPolSquare, YyPolSquare, Vs1PolX, Vs1PolY, arrow_scale=0.0001)
        #print(fig.data[0])
        fig.update_traces(name='Polarization Vectors (on/off)', visible='legendonly')
        fig.update_layout(legend=dict(
            yanchor="top",
            xanchor="left",
            y=0.99,
            x=0.01,
        ))
        fig.add_trace(go.Contour(x=x_data[0], y=y_data, z=z_data, line_smoothing=1, showscale=True,  contours_coloring='heatmap',
            colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        ))
        fig.update_layout(
            yaxis=dict(range=[-2,2]),
            xaxis=dict(range=[-3,3]), 
            showlegend=True,
            height=1000,
            margin=dict(l=20, r=20, t=20, b=20),
            modebar_add='drawclosedpath',
            newshape=dict(fillcolor='turquoise'),
            plot_bgcolor='rgba(0,0,0,0)',
            title={
                'text': 'VS1 Polarization and Splitting Time (s/km) <br> Max : {} Min {} '.format(np.around(np.nanmax(z_data), 4), np.around(np.nanmin(z_data), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
            )
        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[0],
            y=[0],
            marker_symbol='circle-open',
            showlegend=False, 
            marker_size=580,
            marker_line_width=162,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )

        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[-1.75],
            y=[0],
            marker_symbol='square',
            showlegend=False, 
            marker_size=300,
            marker_line_width=30,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )

        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[-1.1],
            y=[1],
            marker_symbol='square',
            showlegend=False, 
            marker_size=50,
            marker_line_width=50,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )

        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[-1.1],
            y=[-1],
            marker_symbol='square',
            showlegend=False, 
            marker_size=50,
            marker_line_width=50,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )

        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[1.75],
            y=[0],
            marker_symbol='square',
            showlegend=False, 
            marker_size=300,
            marker_line_width=30,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )
        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)



        return fig


    if plotType == 'VP' and plotType3D==None:

        thirdInterSet = Vp
        xTest, yTest = np.meshgrid(xMeshInitiate, yMeshInitiate)

        fig = 1

        xx =  Uxyz[:,0] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))
        yy =  Uxyz[:,1] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))

        data = sp.interpolate.griddata((xx,yy), Vp, (xTest, yTest))
        x_data = []
        combine = list(deepflatten(xTest, depth=1))
        x_data.append(combine)

        y_data = np.linspace(-1,1,200)

        z_data = []
        for x in range(0, len(data)):
            combine = list(deepflatten(data[x], depth=1))
            z_data.append(combine)

        fig = go.Figure()
        fig.add_trace(go.Contour(x=x_data[0], y=y_data, z=z_data, line_smoothing=1, showscale=True, contours_coloring='heatmap',
            colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        ))
        
        fig.update_layout(
            yaxis=dict(range=[-2,2]),
            xaxis=dict(range=[-3,3]), 
            showlegend=True,
            margin=dict(l=20, r=20, t=20, b=20),
            height=1000,
            minreducedheight=750,
            minreducedwidth=1500,
            title={
                'text': 'VP <br> Max : {} Min {} '.format(np.around(np.nanmax(z_data), 4), np.around(np.nanmin(z_data), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
            )

        fig.update_layout(dragmode='drawrect',
                    # style of new shapes
                    newshape=dict(line_color='yellow',
                                fillcolor='turquoise',
                                opacity=0.5),
                    plot_bgcolor='rgba(0,0,0,0)')


        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[0],
            y=[0],
            marker_symbol='circle-open',
            showlegend=False, 
            marker_size=580,
            marker_line_width=162,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'

         )
        )

        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)


        return fig


    if plotType == 'VS1' and plotType3D==None:

        thirdInterSet = Vs1

        xTest, yTest = np.meshgrid(xMeshInitiate, yMeshInitiate)

        fig = 1

        xx =  Uxyz[:,0] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))
        yy =  Uxyz[:,1] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))

        data = sp.interpolate.griddata((xx,yy), Vs1, (xTest, yTest))
        x_data = []
        combine = list(deepflatten(xTest, depth=1))
        x_data.append(combine)
        y_data = np.linspace(-1,1,200)

        z_data = []
        for x in range(0, len(data)):
            combine = list(deepflatten(data[x], depth=1))
            z_data.append(combine)

        fig = go.Figure()
        fig.add_trace(go.Contour(x=x_data[0], y=y_data, z=z_data, line_smoothing=1,  showscale=True, contours_coloring='heatmap',
            colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        ))
        fig.update_layout(
            yaxis=dict(range=[-2,2]),
            xaxis=dict(range=[-3,3]), 
            showlegend=True,
            margin=dict(l=20, r=20, t=20, b=20),
            height=1000,
            plot_bgcolor='rgba(0,0,0,0)',
            title={
                'text': 'VS1 <br> Max : {} Min {} '.format(np.around(np.nanmax(z_data), 4), np.around(np.nanmin(z_data), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
            )


        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[0],
            y=[0],
            marker_symbol='circle-open',
            showlegend=False, 
            marker_size=580,
            marker_line_width=162,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )
        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)


        return fig
 

    if plotType == 'VS2' and plotType3D==None:

        thirdInterSet = Vs2

        xTest, yTest = np.meshgrid(xMeshInitiate, yMeshInitiate)

        fig = 1

        xx =  Uxyz[:,0] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))
        yy =  Uxyz[:,1] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))

        data = sp.interpolate.griddata((xx,yy), Vs2, (xTest, yTest))
        x_data = []
        combine = list(deepflatten(xTest, depth=1))
        x_data.append(combine)

        y_data = np.linspace(-1,1,200)


        z_data = []
        for x in range(0, len(data)):
            combine = list(deepflatten(data[x], depth=1))
            z_data.append(combine)

        
        #fig = ff.create_quiver(XxPol, YyPol, Vs1PolX, Vs1PolY)
        fig = go.Figure()
        fig.add_trace(go.Contour(x=x_data[0], y=y_data, z=z_data, line_smoothing=1,  showscale=True, contours_coloring='heatmap',
            colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        ))
        fig.update_layout(
            yaxis=dict(range=[-2,2]),
            xaxis=dict(range=[-3,3]),
            margin=dict(l=20, r=20, t=20, b=20), 
            showlegend=True,
            height=1000,
            plot_bgcolor='rgba(0,0,0,0)',
            title={
                'text': 'VS2 <br> Max : {} Min {} '.format(np.around(np.nanmax(z_data), 4), np.around(np.nanmin(z_data), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
            )

        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[0],
            y=[0],
            marker_symbol='circle-open',
            showlegend=False, 
            marker_size=580,
            marker_line_width=162,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )

        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)
   
        return fig

    if plotType == 'VPVS1' and plotType3D==None:

   
        thirdInterSet = Vp/Vs1

        xTest, yTest = np.meshgrid(xMeshInitiate, yMeshInitiate)

        fig = 1


        xx =  Uxyz[:,0] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))
        yy =  Uxyz[:,1] * np.real(np.sqrt(1/(1-(Uxyz[:,2]))))

        data = sp.interpolate.griddata((xx,yy), Vp/Vs1, (xTest, yTest))
        x_data = []
        combine = list(deepflatten(xTest, depth=1))
        x_data.append(combine)

        y_data = np.linspace(-1,1,200)


        z_data = []
        for x in range(0, len(data)):
            combine = list(deepflatten(data[x], depth=1))
            z_data.append(combine)

        
        #fig = ff.create_quiver(XxPol, YyPol, Vs1PolX, Vs1PolY)
        fig = go.Figure()
        fig.add_trace(go.Contour(x=x_data[0], y=y_data, z=z_data, line_smoothing=1, showscale=True, contours_coloring='heatmap',
            colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        ))
        

        fig.update_layout(
            yaxis=dict(range=[-2,2]),
            xaxis=dict(range=[-3,3]), 
            showlegend=True,
            height=1000,
            margin=dict(l=20, r=20, t=20, b=20),
            plot_bgcolor='rgba(0,0,0,0)',
            title={
                'text': 'VPVS1 <br> Max : {} Min {} '.format(np.around(np.nanmax(z_data), 4), np.around(np.nanmin(z_data), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
            )
        fig.add_trace(
            go.Scatter(
            mode='markers',
            x=[0],
            y=[0],
            marker_symbol='circle-open',
            showlegend=False, 
            marker_size=580,
            marker_line_width=162,
            opacity=1.0,
            marker_line_color='white',
            marker_color='white',
            hoverinfo='skip'
         )
        )

        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)
        return fig


    df = pd.DataFrame(list(zip(XxPol, YyPol, ZzPol, Vs1PolX, Vs1PolY, Vs1PolZ)), columns = ['x','y','z','u','v','w'])
    df = df.drop_duplicates()

    #remove_n = 1300
    #drop_indices = np.random.choice(df.index, remove_n, replace=False)

    #df = df.drop(drop_indices)
    index_list = df.index

      
    #print(colorValues)
    if plotType3D=='3DQuiver' and plotType==None:
        #5/23/23 - x y z as origin point, u w v as end point, connect two 

        fig1 = make_subplots(rows=1, cols=1,
                            specs=[[{'is_3d': True}]],
                            subplot_titles=['Vs1 Polarization & splitting time (s/km)'],
                        )

        fig1.add_trace(go.Surface(x=Xx, y=Yy, z=Zz, surfacecolor=quiver3DColorValues, 
                colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]    
        ))

        Vs1PolXReversed = [x * -1 for x in Vs1PolX]
        Vs1PolYReversed = [x * -1 for x in Vs1PolY]
        Vs1PolZReversed = [x * -1 for x in Vs1PolZ]

        fig1.add_trace(go.Surface(z=z1, y=y1, x=x1, showscale=False))

        fig1.add_trace(go.Cone(
            x = XxPol,
            y = YyPol,
            z = ZzPol,
            u = Vs1PolX,
            v = Vs1PolY,
            w = Vs1PolZ,
            sizeref = 2.5,
            anchor= 'tail',
            showscale=False,
            colorscale="gray",
            cmax=2,
            cmin=1

        ))

        fig1.add_trace(go.Cone(
            x = XxPol,
            y = YyPol,
            z = ZzPol,
            u = Vs1PolXReversed,
            v = Vs1PolYReversed,
            w = Vs1PolZReversed,
            anchor = 'tail',
            sizeref = 2.5,
            showscale=False,
            colorscale="gray",
            cmax=2,
            cmin=1

        ))


        fig1.update_layout(
            scene=dict(
            xaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='X Axis')),
            yaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Y Axis')),
            zaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Z Axis')),
            camera=dict(eye=dict(x=-2, y=-2, z=1.25))
            ),
        )

        return fig1


        
    if plotType3D=='3DVP' and plotType==None:
        vpFig = make_subplots(rows=1, cols=1,
                        specs=[[{'is_3d': True}]],
                     )
        vpFig.add_trace(go.Surface(x=Xx, y=Yy, z=Zz, surfacecolor=velp,
                colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
         ), 1, 1)
        vpFig.add_trace(go.Surface(z=z1, y=y1, x=x1, showscale=False))

        vpFig.update_layout(
            scene=dict(
            xaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='X Axis')),
            yaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Y Axis')),
            zaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Z Axis')),
            camera=dict(eye=dict(x=-2, y=-2, z=1.25))
            ),
            title={
                'text': 'VP <br> Max : {} Min {} '.format(np.around(np.nanmax(velp), 4), np.around(np.nanmin(velp), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
        )
        return vpFig

    if plotType3D=='3DVS1' and plotType==None:
        vs1Fig = make_subplots(rows=1, cols=1,
                          specs=[[{'is_3d': True}]],
                     )
        vs1Fig.add_trace(go.Surface(x=Xx, y=Yy, z=Zz, surfacecolor=vels1,
                colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        
        ), 1, 1)
        vs1Fig.add_trace(go.Surface(z=z1, y=y1, x=x1, showscale=False))
        vs1Fig.update_layout(
            scene=dict(
            xaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='X Axis')),
            yaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Y Axis')),
            zaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Z Axis')),
            camera=dict(eye=dict(x=-2, y=-2, z=1.25))
            ),
            title={
                'text': 'VS1 <br> Max : {} Min {} '.format(np.around(np.nanmax(vels1), 4), np.around(np.nanmin(vels1), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
        )
        return vs1Fig

    if plotType3D=='3DVS2' and plotType==None:
        vs2Fig = make_subplots(rows=1, cols=1,
                            specs=[[{'is_3d': True}]],
                       )
        vs2Fig.add_trace(go.Surface(x=Xx, y=Yy, z=Zz, surfacecolor=vels2,
                colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]
        
        ), 1, 1)
        vs2Fig.add_trace(go.Surface(z=z1, y=y1, x=x1, showscale=False))
        vs2Fig.update_layout(
            scene=dict(
            xaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='X Axis')),
            yaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Y Axis')),
            zaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Z Axis')),
            camera=dict(eye=dict(x=-2, y=-2, z=1.25))
            ),
            title={
                'text': 'VS2 <br> Max : {} Min {} '.format(np.around(np.nanmax(vels2), 4), np.around(np.nanmin(vels2), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
        )
        return vs2Fig
    if plotType3D=='3DVPVS1' and plotType==None:
    
        vpvs1 = np.divide(velp, vels1)

        vpvs1Fig = make_subplots(rows=1, cols=1,
                        specs=[[{'is_3d': True}]],
                    )
        vpvs1Fig.add_trace(go.Surface(x=Xx, y=Yy, z=Zz, surfacecolor=vpvs1, 
                colorscale=[[0.0, "rgb(165,0,38)"],
                [0.1111111111111111, "rgb(215,48,39)"],
                [0.2222222222222222, "rgb(244,109,67)"],
                [0.3333333333333333, "rgb(253,174,97)"],
                [0.4444444444444444, "rgb(254,224,144)"],
                [0.5555555555555556, "rgb(224,243,248)"],
                [0.6666666666666666, "rgb(171,217,233)"],
                [0.7777777777777778, "rgb(116,173,209)"],
                [0.8888888888888888, "rgb(69,117,180)"],
                [1.0, "rgb(49,54,149)"]]), 1, 1)
                
        vpvs1Fig.add_trace(go.Surface(z=z1, y=y1, x=x1, showscale=False))
        vpvs1Fig.update_layout(
            scene=dict(
            xaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='X Axis')),
            yaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Y Axis')),
            zaxis= dict(nticks=5, range=[-2,2], visible=True, showbackground=False, title=dict(font=dict(size=12, color='#000000'), text='Z Axis')),
            camera=dict(eye=dict(x=-2, y=-2, z=1.25))
            ),
            title={
                'text': 'VPVS1 <br> Max : {} Min {} '.format(np.around(np.nanmax(vpvs1), 4), np.around(np.nanmin(vpvs1), 4)),
                'xanchor':'center',
                'yanchor':'top',
                'x': 0.5,
                'y': 0.9
            }
        )

        return vpvs1Fig


def tv_averaging(t1=None, t2=None, t3=None, t4=None, t5=None, user_input=False, userInputDataFrame=None, weight=None):
    christoffelParameters = ChristoffelDevParams.ChristoffelParameters()
    Res = 3
    pi = math.pi
    aggTensorColumns = christoffelParameters.columnNamesTensorCSV
    sphereColumnNames = [x for x in range(0, 121)]
    csv_tensor_number = [x for x in range(0, 50)]
    user_tensor_number = [x for x in range(97, 147)]
    user_to_csv_mapping = zip(user_tensor_number, csv_tensor_number)
    user_to_csv_mapping = dict(user_to_csv_mapping)


    returnedUserInputDataFrame = returnUserInputDataFrame(userInputDataFrame)

    #aggTensorCSVLocation = 'databases/AllAggTensorsCSV.csv'
    aggTensorCSVLocation = 'databases/DB_Tensors_Rhos_OrigOrientation_CSV.csv'
    tensors_to_average=[t1-1 ,t2-1]

    if t3 != None:
        tensors_to_average = [t1-1 , t2-1 , t3-1]

    if t4 != None:
        tensors_to_average = [t1-1, t2-1, t3-1,t4-1 ]

    if t5 != None:
        tensors_to_average = [t1-1, t2-1, t3-1, t4-1, t5-1]


    numberOfTensorsToAverage = len(tensors_to_average)
    tensorDF = pd.read_csv(aggTensorCSVLocation, sep=',', engine='python')
    #print(tensorDF)

    if user_input == True:
        print(user_input)
        tensorDF = returnedUserInputDataFrame
        Cijs = np.zeros((5,21), dtype=float)
        Rhos = np.zeros((5,1), dtype=float)
        print(tensors_to_average)
        for i in range(0,numberOfTensorsToAverage):
            for j in range(4,25):
                Cijs[i,j-4] = tensorDF.iloc[tensors_to_average[i]+1,j]
        for i in range(0,numberOfTensorsToAverage):
            for j in range(3,4):
                Rhos[i,j-3] = tensorDF.iloc[tensors_to_average[i]+1,j]

    else:
        Cijs = np.zeros((5,21), dtype=float)
        Rhos = np.zeros((5,1), dtype=float)

        for i in range(0,numberOfTensorsToAverage):
            for j in range(4,25):
                Cijs[i,j-4] = tensorDF.iloc[tensors_to_average[i],j]

        for i in range(0,numberOfTensorsToAverage):
            for j in range(3,4):
                Rhos[i,j-3] = tensorDF.iloc[tensors_to_average[i],j]



    Equal = 0
    #make sure to add customization for weight , options able to be selected from dropdown? input?
    wt = []
    for counter in range(0, numberOfTensorsToAverage):
        wt.append(weight[counter])
    
    #print(Cijs[0,:])

    wt = wt / np.sum(wt)


    if Equal == 1:
        wt = np.ones((3,1)*(1/3))

    tensorsV = []
    tensorsR = []
    V_ave = np.zeros((36,1), dtype=float)
    R_ave = np.zeros((36,1), dtype=float)


    
    for array in range(0,5):
        tensorsV.append(np.zeros((6,6),dtype=float))
        tensorsR.append(np.zeros((6,6),dtype=float))
    for i in range(0,5):
        C = Cijs[i,:]
        #all indicies below are subtracted by 1 for indexing purposes
        tensorsV[i] = [[C[0], C[1], C[2], C[3], C[4], C[5]],
        [C[1], C[6], C[7], C[8], C[9], C[10]],
        [C[2], C[7], C[11], C[12], C[13], C[14]],
        [C[3], C[8], C[12], C[15], C[16], C[17]],
        [C[4], C[9], C[13], C[16], C[18], C[19]],
        [C[5], C[10], C[14], C[17], C[19], C[20]]
        ]
        tensorsV[i] = np.asarray(tensorsV[i])

    # CHECK THE STUPID INVERSE FUNCTION
    for i in range(0,5):
        tensorsR[i] = np.linalg.pinv(tensorsV[i])

    
    for i in range(0,36):
        PV = np.zeros((5,1), dtype=float)
        PR = np.zeros((5,1), dtype=float)
        for j in range(0,numberOfTensorsToAverage):
            flattenedVTensor = tensorsV[j].flatten()
            flattenedRTensor = tensorsR[j].flatten()
            PV[j,] = (wt[j])*flattenedVTensor[i]
            PR[j,] = (wt[j])*flattenedRTensor[i]
        V_ave[i] = np.sum(PV)
        R_ave[i] = np.sum(PR)
    test_array = np.zeros((6,6), dtype=float)

    #re-do average arrays
    
    V_ave = np.reshape(V_ave, (6,6))
    R_ave = np.reshape(R_ave, (6,6))
    V_ave = V_ave.flatten()
    V_ave = V_ave.tolist()

    V_ave_rounded = []
    R_ave_rounded = []

    for x in V_ave:
        y = round(x,4)
        V_ave_rounded.append(y)

    return V_ave_rounded, R_ave





def decomposition(tensor_index, user_input=False, plotType='None', userInputDataFrame=None):



    christoffelParameters = ChristoffelDevParams.ChristoffelParameters()
    Res = 3
    pi = math.pi
    aggTensorColumns = christoffelParameters.columnNamesTensorCSV
    sphereColumnNames = [x for x in range(0, 121)]
    csv_tensor_number = [x for x in range(0, 50)]
    user_tensor_number = [x for x in range(97, 147)]
    user_to_csv_mapping = zip(user_tensor_number, csv_tensor_number)
    user_to_csv_mapping = dict(user_to_csv_mapping)

    #aggTensorCSVLocation = 'databases/AllAggTensorsCSV.csv'
    aggTensorCSVLocation = 'databases/DB_Tensors_Rhos_OrigOrientation_CSV.csv'

    sampleN = 1
    VRHNumber = 1
    C = []

    minNum = christoffelParameters.minNum
    Samps = 3

    tensorDF = pd.read_csv(aggTensorCSVLocation, sep=',', engine='python')
    FigName = 'Sample {} Phase={} VRH = {}'.format(str(sampleN), str(minNum), str(VRHNumber))


    tensorIndex = int(tensor_index)

    if user_input == True:
        returnedUserInputDataFrame = returnUserInputDataFrame(userInputDataFrame)
        tensorDF = returnedUserInputDataFrame
        tensorList = (tensorDF.iloc[tensorIndex, 4:25])
        Rhos = (tensorDF.iloc[tensorIndex, 3])
    else:
        tensorList = (tensorDF.iloc[tensorIndex-1, 4:25])
        Rhos = (tensorDF.iloc[tensorIndex-1, 3])

    HexA = []
    HexAA = []




    breakDown = np.zeros((1,32))
    symmetryAxis = np.zeros((1, 21))
    isoComponent = np.zeros((1, 21))
    hexagonalComponent = np.zeros((1, 21))
    orthorhombicComponent = np.zeros((1, 21))



    for s in range (0, 1):
        cInitial = [[tensorList[0],tensorList[1],tensorList[2],tensorList[3],tensorList[4],tensorList[5]],
             [tensorList[1],tensorList[6],tensorList[7],tensorList[8],tensorList[9],tensorList[10]],
             [tensorList[2],tensorList[7],tensorList[11],tensorList[12],tensorList[13],tensorList[14]],
             [tensorList[3],tensorList[8],tensorList[12],tensorList[15],tensorList[16],tensorList[17]],
             [tensorList[4],tensorList[9],tensorList[13],tensorList[16],tensorList[18],tensorList[19]],
             [tensorList[5],tensorList[10],tensorList[14],tensorList[17],tensorList[19],tensorList[20]]
             ]
            
        dij = [[cInitial[0][0]+cInitial[0][1]+cInitial[0][2], cInitial[0][5]+cInitial[1][5]+cInitial[2][5], cInitial[0][4]+cInitial[1][4]+cInitial[2][4] ],
            [cInitial[0][5]+cInitial[1][5]+cInitial[2][5], cInitial[0][1]+cInitial[1][1]+cInitial[2][1], cInitial[0][3]+cInitial[1][3]+cInitial[2][3] ],
            [cInitial[0][4]+cInitial[1][4]+cInitial[2][4], cInitial[0][3]+cInitial[1][3]+cInitial[2][3], cInitial[0][2]+cInitial[1][2]+cInitial[2][2] ] 
        ]
     

        vij = [[ cInitial[0][0]+cInitial[5][5]+cInitial[4][4], cInitial[0][5]+cInitial[1][5]+cInitial[3][4], cInitial[0][4]+cInitial[2][4]+cInitial[3][5]  ],
        [ cInitial[0][5]+cInitial[1][5]+cInitial[3][4], cInitial[5][5]+cInitial[1][1]+cInitial[3][3], cInitial[1][3]+cInitial[2][3]+cInitial[4][5] ],
        [ cInitial[0][4]+cInitial[2][4]+cInitial[3][5], cInitial[1][3]+cInitial[2][3]+cInitial[4][5], cInitial[4][4]+cInitial[3][3]+cInitial[2][2] ]
        ]



        eigenVals_dij, eigenVecs_dij = np.linalg.eig(dij)
        eigenVals_vij, eigenVecs_vij = np.linalg.eig(vij)


        idx_vij = eigenVals_vij.argsort()
        idx_dij = eigenVals_dij.argsort()

        idy_vij = eigenVals_vij.argsort()
        idy_dij = eigenVals_dij.argsort()




        eigenVals_vij = eigenVals_vij[idx_vij]
        eigenVals_dij = eigenVals_dij[idx_dij]


        eigenVecs_vij = eigenVecs_vij[:, idy_vij]
        eigenVecs_dij = eigenVecs_dij[:, idy_dij]




        
    
        eigs = eigenVals_dij
        eigsTwo = eigenVals_vij


        aOne = eigs[1]-eigs[0]
        aTwo = eigs[2]-eigs[1]
       

        aOneTwo = eigsTwo[1]- eigsTwo[0]
        aTwoTwo = eigsTwo[2]- eigsTwo[1]
        

        if aOne > aTwo:
            eigenVecs_dij = eigenVecs_dij[:, [2,1,0]]
            eigs = eigs[[2,1,0]]

        if aOneTwo > aTwoTwo:
            eigenVecs_vij  = eigenVecs_vij[:, [2,1,0]]
            eigsTwo = eigsTwo[[2,1,0]]
        

        #correct to here - 5/23/23 - 1241pm
   


        Ang = np.empty([3,3])
        AngTwo = np.empty([3,3])
        
        for i in range(0, 3):
            for j in range(0, 3):
                Ang[i,j] = np.rad2deg(np.arccos(np.dot(eigenVecs_vij[:,i], eigenVecs_dij[:,j])))
                AngTwo[i,j] = np.rad2deg(np.arccos(np.dot(eigenVecs_vij[:,i], -1*eigenVecs_dij[:,j])))
                #Ang[i,j] = math.degrees(math.acos(np.dot(eigenVecs_vij[:,i], eigenVecs_dij[:,j])))
                #AngTwo[i,j] = math.degrees(math.acos(np.dot(eigenVecs_vij[:,i], -1*eigenVecs_dij[:,j])))


        
        

        eVd = eigenVecs_dij
        AngF = Ang
        for i in range(0,3):
            if np.amin(AngTwo[:,i]) < np.amin(Ang[:,i]):
                eVd[:,i] = -1*eVd[:,i]
                AngF[:,i] = AngTwo[:,i]
        

        Bisectrix = np.empty([3,3])
        BA = np.empty([3,3])

        for i in range(0,3):
            I = np.where(AngF[:,i]==np.min(AngF[:,i]))
            Bisectrix[:,i] = np.divide((np.add(eigenVecs_vij[:,I[0][0]], eVd[:,i])), 2)
            Bisectrix[:,i] = Bisectrix[:,i] / np.linalg.norm(Bisectrix[:,i])

        Bisectrix = np.column_stack((Bisectrix[:,1], Bisectrix[:,0], Bisectrix[:,2]))
      

        
        

        

        BisectrixOne = Bisectrix
        BisectrixOne[:,1] = np.cross(BisectrixOne[:,2], Bisectrix[:,1])
        BisectrixOne[:,0] = np.cross(BisectrixOne[:,0], Bisectrix[:,2])
        
        BisectrixAll = np.array([BisectrixOne, BisectrixOne[:,[0,1,2]], BisectrixOne[:,[1,0,2]], BisectrixOne[:, [1,2,0]], BisectrixOne[:,[2,0,1]], BisectrixOne[:,[2,1,0]] ])


        ### line 106-111 in matlab code
        IsoA = np.zeros([6,1])
        #HexA = np.zeros([6,1])
        TetA = np.zeros([6,1])
        OrthoA = np.zeros([6,1])
        MonA = np.zeros([6,1])
        TricA = np.zeros([6,1])
       



        for btrix in range (0, len(BisectrixAll)):

            K2 = np.empty([3,3])
            K3 = np.empty([3,3])
            K4 = np.empty([3,3])
            Omega = BisectrixAll[btrix].conj().transpose()

            K1 = Omega**2

            for i in range (0,3):
                for j in range(0, 3):
                    K2[i,j] = Omega[i, round((j+1)%3)] * Omega[i, round((j+2)%3)]
                    K3[i,j] = Omega[round((i+1)%3), j ]*Omega[round((i+2)%3), j]
                    K4[i,j] = Omega[round((i+1)%3) ,round((j+1)%3) ]* Omega[round((i+2)%3), round((j+2)%3)]  + (Omega[round((i+1)%3),round((j+2)%3)]*Omega[round((i+2)%3), round((j+1)%3)])

                    
            K_a = np.concatenate((K1, 2*K2), axis=1)
            K_b = np.concatenate((K3, K4), axis=1)
            K = np.concatenate((K_a, K_b))
            #C = K * cInitial * K.conj().transpose()
            C = np.matmul(K, cInitial)
            C = np.matmul(C, K.conj().transpose())


            #symmetryAxis = np.empty((s,1))
            
            symmetry = [C[0][0:6], C[1][1:6], C[2][2:6],C[3][3:6], C[4][4:6]]
            symmetryAxis = []
            for x in symmetry:
                for y in x:
                    symmetryAxis.append(y)
            symmetryAxis.append(C[5][5])
            symmetryAxis = np.array(symmetryAxis)
            symmetryAxis = np.reshape(symmetryAxis, [1,21])



            
            

            ###Note 1/11/22 at office - check Vec_1 
            Vec_1 = [[C[0][0]], [C[1][1]], [C[2][2]], 
            [math.sqrt(2)*C[1][2]], [math.sqrt(2)*C[0][2]], [math.sqrt(2)*C[0][1]], 
            [2*C[3][3]], [2*C[4][4]], [2*C[5][5]], 
            [2*C[0][3]], [2*C[1][4]], [2*C[2][5]], [2*C[2][3]], [2*C[0][4]], 
            [2*C[1][5]], [2*C[1][3]], 
            [2*C[2][4]], [2*C[0][5]], [2*math.sqrt(2)*C[4][5]], [2*math.sqrt(2)*C[3][5]], 
            [2*math.sqrt(2)*C[3][4]]
            ]

            Vec_1 = np.array(Vec_1)

        
            Iso_proj = np.zeros([21,21])
            Iso_proj[0:3, 0:3] = 3/15
            

            Iso_proj[3:6, 3:6] = 4/15
    

            Iso_proj[6:9, 6:9] = 1/5


            Iso_proj[0:3, 6:9] = 2/15

    
            Iso_proj[6:9, 0:3] = 2/15


            Iso_proj[0:3, 3:6] = math.sqrt(2)/15


            Iso_proj[3:6, 0:3] = math.sqrt(2)/15

            Iso_proj[3:6, 6:9] = -1*math.sqrt(2)/15

            Iso_proj[6:9, 3:6] = -1*math.sqrt(2)/15


            Vec_iso = np.matmul(Vec_1.conj().transpose(), Iso_proj)

            



            MultsList = [1,1,1,1/math.sqrt(2), 1/math.sqrt(2) ,1/math.sqrt(2),1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2, 1/(2*math.sqrt(2)), 1/(2*math.sqrt(2)), 1/(2*math.sqrt(2))]
            Mults = np.zeros((1,21))
            for x in range (0,21):
                Mults[0,x] = MultsList[x]


            E_iso = Vec_iso * Mults
            


            Iso_tensor= [[E_iso[0][0], E_iso[0][5], E_iso[0][4], E_iso[0][9], E_iso[0][13], E_iso[0][17]],
            [E_iso[0][5], E_iso[0][1], E_iso[0][3], E_iso[0][15], E_iso[0][10], E_iso[0][14]],
            [E_iso[0][4], E_iso[0][3], E_iso[0][2], E_iso[0][12], E_iso[0][16], E_iso[0][11]],
            [E_iso[0][9], E_iso[0][15], E_iso[0][12], E_iso[0][6], E_iso[0][20], E_iso[0][19]],
            [E_iso[0][13], E_iso[0][10], E_iso[0][16], E_iso[0][20], E_iso[0][7], E_iso[0][18]],
            [E_iso[0][17], E_iso[0][14], E_iso[0][11], E_iso[0][19], E_iso[0][18], E_iso[0][8] ]]
            #print(Iso_tensor)


            Iso_tensor = np.array(Iso_tensor)
            isoComponent = np.zeros([1,21])
            #print(Iso_tensor[0][0:8])

        

            isoComponent[0,0:6] = Iso_tensor[0][0:6]
            isoComponent[0,6:11] = Iso_tensor[1][1:6]
            isoComponent[0, 11:15] = Iso_tensor[2][2:6]
            isoComponent[0, 15:18] = Iso_tensor[3][3:6] 
            isoComponent[0, 18:20] = Iso_tensor[4][4:6] 
            isoComponent[0, 20:21] = Iso_tensor[5][5:6]





            Vec_2 = np.subtract(Vec_1.conj().T, Vec_iso)

            Iso = (1-(np.divide((np.linalg.norm(Vec_2)), np.linalg.norm(Vec_1))))*100

            Hex_proj = np.zeros([21,21])

            Hex_proj[0:2, 0:2] = 3/8
            Hex_proj[0:2, 5] = 1 / (4 * math.sqrt(2))
            Hex_proj[0:2, 8] = 1/4
            Hex_proj[2, 2] = 1

            Hex_proj[3:5, 3:5] = 1/2

            Hex_proj[5, 0:2] = 1/ (4*math.sqrt(2))

            Hex_proj[5, 5] = 3/4

            Hex_proj[5, 8] = -1 * 1 / (2*math.sqrt(2))
            Hex_proj[6:8, 6:8] = 1/2
            Hex_proj[8, 0:2] = 1/4
            Hex_proj[8, 5] = -1 * 1/(2*math.sqrt(2))
            Hex_proj[8,8] = 1/2

            Vec_Hex = np.matmul(Vec_2,Hex_proj)

            E_hex = Vec_Hex * Mults

 


            Hex_tensor = [[E_hex[0][0], E_hex[0][5], E_hex[0][4], E_hex[0][9], E_hex[0][13], E_hex[0][17]],

            [E_hex[0][5], E_hex[0][1], E_hex[0][3], E_hex[0][15], E_hex[0][10], E_hex[0][14]],
            [E_hex[0][4], E_hex[0][3], E_hex[0][2], E_hex[0][12], E_hex[0][16], E_hex[0][11]],
            [E_hex[0][9], E_hex[0][15], E_hex[0][12], E_hex[0][6], E_hex[0][20], E_hex[0][19]],
            [E_hex[0][13], E_hex[0][10], E_hex[0][16], E_hex[0][20], E_hex[0][7], E_hex[0][18]],
            [E_hex[0][17], E_hex[0][14], E_hex[0][11], E_hex[0][19], E_hex[0][18], E_hex[0][8]]]

            hexagonalComponent[0,0:6] = Hex_tensor[0][0:6]
            hexagonalComponent[0,6:11] = Hex_tensor[1][1:6]
            hexagonalComponent[0, 11:15] = Hex_tensor[2][2:6]
            hexagonalComponent[0, 15:18] = Hex_tensor[3][3:6] 
            hexagonalComponent[0, 18:20] = Hex_tensor[4][4:6] 
            hexagonalComponent[0, 20:21] = Hex_tensor[5][5:6]


            Diff = np.empty([1, len(BisectrixAll)])

            Diff[:, btrix] = np.sqrt(np.sum( (hexagonalComponent[s,:] + isoComponent[s,:]) - symmetryAxis[s,:]) **2)
            
            Vec_3 = Vec_2 - Vec_Hex

            Hex = (1-(np.divide((np.linalg.norm(Vec_3)), np.linalg.norm(Vec_2))))*(100-Iso)
            HexA.append(Hex)
            #print(Hex)
            #HexA = np.empty([1,6])
            #HexA[0,:] = Hex

            Tet_proj = np.zeros([21,21])

            Tet_proj[0:2, 0:2] = 1/2
            Tet_proj[3:5, 3:5] = 1/2
            Tet_proj[6:8, 6:8] = 1/2
            Tet_proj[2,2] = 1
            Tet_proj[5,5] = 1
            Tet_proj[8,8] = 1

            Vec_Tet = np.matmul(Vec_3,Tet_proj)
            E_Tet = Vec_Tet * Mults
            Tet_tensor = [[E_Tet[0][0], E_Tet[0][5], E_Tet[0][4], E_Tet[0][9], E_Tet[0][13], E_Tet[0][17]],

            [E_Tet[0][5], E_Tet[0][1], E_Tet[0][3], E_Tet[0][15], E_Tet[0][10], E_Tet[0][14]],
            [E_Tet[0][4], E_Tet[0][3], E_Tet[0][2], E_Tet[0][12], E_Tet[0][16], E_Tet[0][11]],
            [E_Tet[0][9], E_Tet[0][15], E_Tet[0][12], E_Tet[0][6], E_Tet[0][20], E_Tet[0][19]],
            [E_Tet[0][13], E_Tet[0][10], E_Tet[0][16], E_Tet[0][20], E_Tet[0][7], E_Tet[0][18]],
            [E_Tet[0][17], E_Tet[0][14], E_Tet[0][11], E_Tet[0][19], E_Tet[0][18], E_Tet[0][8]]]
        
            Vec_4 = Vec_3 - Vec_Tet

            Tet = (1-(np.divide((np.linalg.norm(Vec_4)), np.linalg.norm(Vec_3))))*(100-(Iso+Hex))

            Ortho_proj = np.zeros([21,21])
            ### Need to preserve shape of original array - can't really use np.diag here (7/25/2023)
            for x in range(0, 9):
                Ortho_proj[x,x] = 1

            Vec_Ortho = np.matmul(Vec_4, Ortho_proj)
            E_Ortho = Vec_Ortho * Mults

            Ortho_tensor = [[E_Ortho[0][0], E_Ortho[0][5], E_Ortho[0][4], E_Ortho[0][9], E_Ortho[0][13], E_Ortho[0][17]],

            [E_Ortho[0][5], E_Ortho[0][1], E_Ortho[0][3], E_Ortho[0][15], E_Ortho[0][10], E_Ortho[0][14]],
            [E_Ortho[0][4], E_Ortho[0][3], E_Ortho[0][2], E_Ortho[0][12], E_Ortho[0][16], E_Ortho[0][11]],
            [E_Ortho[0][9], E_Ortho[0][15], E_Ortho[0][12], E_Ortho[0][6], E_Ortho[0][20], E_Ortho[0][19]],
            [E_Ortho[0][13], E_Ortho[0][10], E_Ortho[0][16], E_Ortho[0][20], E_Ortho[0][7], E_Ortho[0][18]],
            [E_Ortho[0][17], E_Ortho[0][14], E_Ortho[0][11], E_Ortho[0][19], E_Ortho[0][18], E_Ortho[0][8]]]
        
            Vec_5 = Vec_4 - Vec_Ortho

            Ortho = (1-(np.divide((np.linalg.norm(Vec_5)), np.linalg.norm(Vec_4))))*(100-(Iso+Hex+Tet))
            
            OrthoA = np.empty([1,6])
            OrthoA[0,:] = Ortho

            if np.linalg.norm(Vec_5) == 0:
                Mon = 0
                Tric = 0
            else:

                Mono_proj = np.zeros([21,21])
                for x in range(0, 9):
                    Mono_proj[x,x] = 1

                Mono_proj[11, 11] = 1
                Mono_proj[14,14] = 1
                Mono_proj[17, 17] = 1
                Mono_proj[20, 20] = 1

                Vec_Mono = np.matmul(Vec_5, Mono_proj)
                E_Mono = Vec_Mono * Mults
                Mono_tensor = [[E_Mono[0][0], E_Mono[0][5], E_Mono[0][4], E_Mono[0][9], E_Mono[0][13], E_Mono[0][17]],

                [E_Mono[0][5], E_Mono[0][1], E_Mono[0][3], E_Mono[0][15], E_Mono[0][10], E_Mono[0][14]],

                [E_Mono[0][4], E_Mono[0][3], E_Mono[0][2], E_Mono[0][12], E_Mono[0][16], E_Mono[0][11]],

                [E_Mono[0][9], E_Mono[0][15], E_Mono[0][12], E_Mono[0][6], E_Mono[0][20], E_Mono[0][19]],

                [E_Mono[0][13], E_Mono[0][10], E_Mono[0][16], E_Mono[0][20], E_Mono[0][7], E_Mono[0][18]],
                
                [E_Mono[0][17], E_Mono[0][14], E_Mono[0][11], E_Mono[0][19], E_Mono[0][18], E_Mono[0][8]]]
            
                Vec_6 = Vec_5 - Vec_Mono
                Mono = (1-(np.divide((np.linalg.norm(Vec_6)), np.linalg.norm(Vec_5))))*(100-(Iso+Hex+Tet+Ortho))
            Tric = 100 - (Iso+Hex+Tet+Ortho+Mono)

            

            I1 = np.where(HexA == max(HexA))
            I2 = np.where(HexA == min(HexA))

                
                

            Bisectrix[:,2] = BisectrixAll[I1[0][0]][:,2]
            Bisectrix[:,0] = BisectrixAll[I2[0][0]][:,2]
            Bisectrix[:,1] = np.cross(Bisectrix[:,2], Bisectrix[:,0])
            ### beginning decomp again - 7/26/2023
            K2 = np.empty([3,3])
            K3 = np.empty([3,3])
            K4 = np.empty([3,3])
            
            Omega = Bisectrix.conj().transpose()

            K1 = Omega**2

            for i in range (0,3):
                for j in range(0, 3):
                    K2[i,j] = Omega[i, round((j+1)%3)] * Omega[i, round((j+2)%3)]
                    K3[i,j] = Omega[round((i+1)%3), j ]*Omega[round((i+2)%3), j]
                    K4[i,j] = Omega[round((i+1)%3) ,round((j+1)%3) ]* Omega[round((i+2)%3), round((j+2)%3)]  + (Omega[round((i+1)%3),round((j+2)%3)]*Omega[round((i+2)%3), round((j+1)%3)])

                    
            K_a = np.concatenate((K1, 2*K2), axis=1)
            K_b = np.concatenate((K3, K4), axis=1)
            K = np.concatenate((K_a, K_b))
            #C = K * cInitial * K.conj().transpose()
            C = np.matmul(K, cInitial)
            C = np.matmul(C, K.conj().transpose())


            #symmetryAxis = np.empty((s,1))
            
            symmetry = [C[0][0:6], C[1][1:6], C[2][2:6],C[3][3:6], C[4][4:6]]
            symmetryAxis = []
            for x in symmetry:
                for y in x:
                    symmetryAxis.append(y)
            symmetryAxis.append(C[5][5])
            symmetryAxis = np.array(symmetryAxis)
            symmetryAxis = np.reshape(symmetryAxis, [1,21])


            ###Note 1/11/22 at office - check Vec_1 
            Vec_1 = [[C[0][0]], [C[1][1]], [C[2][2]], 
            [math.sqrt(2)*C[1][2]], [math.sqrt(2)*C[0][2]], [math.sqrt(2)*C[0][1]], 
            [2*C[3][3]], [2*C[4][4]], [2*C[5][5]], 
            [2*C[0][3]], [2*C[1][4]], [2*C[2][5]], [2*C[2][3]], [2*C[0][4]], 
            [2*C[1][5]], [2*C[1][3]], 
            [2*C[2][4]], [2*C[0][5]], [2*math.sqrt(2)*C[4][5]], [2*math.sqrt(2)*C[3][5]], 
            [2*math.sqrt(2)*C[3][4]]
            ]

            Vec_1 = np.array(Vec_1)

        
            Iso_proj = np.zeros([21,21])
            Iso_proj[0:3, 0:3] = 3/15
            

            Iso_proj[3:6, 3:6] = 4/15
    

            Iso_proj[6:9, 6:9] = 1/5


            Iso_proj[0:3, 6:9] = 2/15

    
            Iso_proj[6:9, 0:3] = 2/15


            Iso_proj[0:3, 3:6] = math.sqrt(2)/15


            Iso_proj[3:6, 0:3] = math.sqrt(2)/15

            Iso_proj[3:6, 6:9] = -1*math.sqrt(2)/15

            Iso_proj[6:9, 3:6] = -1*math.sqrt(2)/15


            Vec_iso = np.matmul(Vec_1.conj().transpose(), Iso_proj)
            
            MultsList = [1,1,1,1/math.sqrt(2), 1/math.sqrt(2) ,1/math.sqrt(2),1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2, 1/(2*math.sqrt(2)), 1/(2*math.sqrt(2)), 1/(2*math.sqrt(2))]
            Mults = np.zeros((1,21))
            for x in range (0,21):
                Mults[0,x] = MultsList[x]


            E_iso = Vec_iso * Mults
            Iso_tensor= [[E_iso[0][0], E_iso[0][5], E_iso[0][4], E_iso[0][9], E_iso[0][13], E_iso[0][17]],
            [E_iso[0][5], E_iso[0][1], E_iso[0][3], E_iso[0][15], E_iso[0][10], E_iso[0][14]],
            [E_iso[0][4], E_iso[0][3], E_iso[0][2], E_iso[0][12], E_iso[0][16], E_iso[0][11]],
            [E_iso[0][9], E_iso[0][15], E_iso[0][12], E_iso[0][6], E_iso[0][20], E_iso[0][19]],
            [E_iso[0][13], E_iso[0][10], E_iso[0][16], E_iso[0][20], E_iso[0][7], E_iso[0][18]],
            [E_iso[0][17], E_iso[0][14], E_iso[0][11], E_iso[0][19], E_iso[0][18], E_iso[0][8] ]]
            #print(Iso_tensor)


            Iso_tensor = np.array(Iso_tensor)
            isoComponent = np.zeros([1,21])
            #print(Iso_tensor[0][0:8])

        

            isoComponent[0,0:6] = Iso_tensor[0][0:6]
            isoComponent[0,6:11] = Iso_tensor[1][1:6]
            isoComponent[0, 11:15] = Iso_tensor[2][2:6]
            isoComponent[0, 15:18] = Iso_tensor[3][3:6] 
            isoComponent[0, 18:20] = Iso_tensor[4][4:6] 
            isoComponent[0, 20:21] = Iso_tensor[5][5:6]





            Vec_2 = np.subtract(Vec_1.conj().T, Vec_iso)

            Iso = (1-(np.divide((np.linalg.norm(Vec_2)), np.linalg.norm(Vec_1))))*100

            Hex_proj = np.zeros([21,21])

            Hex_proj[0:2, 0:2] = 3/8
            Hex_proj[0:2, 5] = 1 / (4 * math.sqrt(2))
            Hex_proj[0:2, 8] = 1/4
            Hex_proj[2, 2] = 1

            Hex_proj[3:5, 3:5] = 1/2

            Hex_proj[5, 0:2] = 1/ (4*math.sqrt(2))

            Hex_proj[5, 5] = 3/4

            Hex_proj[5, 8] = -1 * 1 / (2*math.sqrt(2))
            Hex_proj[6:8, 6:8] = 1/2
            Hex_proj[8, 0:2] = 1/4
            Hex_proj[8, 5] = -1 * 1/(2*math.sqrt(2))
            Hex_proj[8,8] = 1/2

            Vec_Hex = np.matmul(Vec_2,Hex_proj)

            E_hex = Vec_Hex * Mults

            Hex_tensor = [[E_hex[0][0], E_hex[0][5], E_hex[0][4], E_hex[0][9], E_hex[0][13], E_hex[0][17]],

            [E_hex[0][5], E_hex[0][1], E_hex[0][3], E_hex[0][15], E_hex[0][10], E_hex[0][14]],
            [E_hex[0][4], E_hex[0][3], E_hex[0][2], E_hex[0][12], E_hex[0][16], E_hex[0][11]],
            [E_hex[0][9], E_hex[0][15], E_hex[0][12], E_hex[0][6], E_hex[0][20], E_hex[0][19]],
            [E_hex[0][13], E_hex[0][10], E_hex[0][16], E_hex[0][20], E_hex[0][7], E_hex[0][18]],
            [E_hex[0][17], E_hex[0][14], E_hex[0][11], E_hex[0][19], E_hex[0][18], E_hex[0][8]]]

            hexagonalComponent[0,0:6] = Hex_tensor[0][0:6]
            hexagonalComponent[0,6:11] = Hex_tensor[1][1:6]
            hexagonalComponent[0, 11:15] = Hex_tensor[2][2:6]
            hexagonalComponent[0, 15:18] = Hex_tensor[3][3:6] 
            hexagonalComponent[0, 18:20] = Hex_tensor[4][4:6] 
            hexagonalComponent[0, 20:21] = Hex_tensor[5][5:6]
            hexagonalComponent = hexagonalComponent + isoComponent

            Vec_3 = Vec_2 - Vec_Hex

            Hex = (1-(np.divide((np.linalg.norm(Vec_3)), np.linalg.norm(Vec_2))))*(100-Iso)
            HexAA.append(Hex)

            Tet_proj = np.zeros([21,21])

            Tet_proj[0:2, 0:2] = 1/2
            Tet_proj[3:5, 3:5] = 1/2
            Tet_proj[6:8, 6:8] = 1/2
            Tet_proj[2,2] = 1
            Tet_proj[5,5] = 1
            Tet_proj[8,8] = 1

            Vec_Tet = np.matmul(Vec_3,Tet_proj)
            E_Tet = Vec_Tet * Mults
            Tet_tensor = [[E_Tet[0][0], E_Tet[0][5], E_Tet[0][4], E_Tet[0][9], E_Tet[0][13], E_Tet[0][17]],

            [E_Tet[0][5], E_Tet[0][1], E_Tet[0][3], E_Tet[0][15], E_Tet[0][10], E_Tet[0][14]],
            [E_Tet[0][4], E_Tet[0][3], E_Tet[0][2], E_Tet[0][12], E_Tet[0][16], E_Tet[0][11]],
            [E_Tet[0][9], E_Tet[0][15], E_Tet[0][12], E_Tet[0][6], E_Tet[0][20], E_Tet[0][19]],
            [E_Tet[0][13], E_Tet[0][10], E_Tet[0][16], E_Tet[0][20], E_Tet[0][7], E_Tet[0][18]],
            [E_Tet[0][17], E_Tet[0][14], E_Tet[0][11], E_Tet[0][19], E_Tet[0][18], E_Tet[0][8]]]
        
            Vec_4 = Vec_3 - Vec_Tet

            Tet = (1-(np.divide((np.linalg.norm(Vec_4)), np.linalg.norm(Vec_3))))*(100-(Iso+Hex))
            
            
            Ortho_proj = np.zeros([21,21])
            ### Need to preserve shape of original array - can't really use np.diag here (7/25/2023)
            for x in range(0, 9):
                Ortho_proj[x,x] = 1

            Vec_Ortho = np.matmul(Vec_4, Ortho_proj)
            E_Ortho = Vec_Ortho * Mults

            Ortho_tensor = [[E_Ortho[0][0], E_Ortho[0][5], E_Ortho[0][4], E_Ortho[0][9], E_Ortho[0][13], E_Ortho[0][17]],

            [E_Ortho[0][5], E_Ortho[0][1], E_Ortho[0][3], E_Ortho[0][15], E_Ortho[0][10], E_Ortho[0][14]],
            [E_Ortho[0][4], E_Ortho[0][3], E_Ortho[0][2], E_Ortho[0][12], E_Ortho[0][16], E_Ortho[0][11]],
            [E_Ortho[0][9], E_Ortho[0][15], E_Ortho[0][12], E_Ortho[0][6], E_Ortho[0][20], E_Ortho[0][19]],
            [E_Ortho[0][13], E_Ortho[0][10], E_Ortho[0][16], E_Ortho[0][20], E_Ortho[0][7], E_Ortho[0][18]],
            [E_Ortho[0][17], E_Ortho[0][14], E_Ortho[0][11], E_Ortho[0][19], E_Ortho[0][18], E_Ortho[0][8]]]

            orthorhombicComponent[0,0:6] = Ortho_tensor[0][0:6]
            orthorhombicComponent[0,6:11] = Ortho_tensor[1][1:6]
            orthorhombicComponent[0, 11:15] = Ortho_tensor[2][2:6]
            orthorhombicComponent[0, 15:18] = Ortho_tensor[3][3:6] 
            orthorhombicComponent[0, 18:20] = Ortho_tensor[4][4:6] 
            orthorhombicComponent[0, 20:21] = Ortho_tensor[5][5:6]
            orthorhombicComponent = hexagonalComponent + orthorhombicComponent
        
            Vec_5 = Vec_4 - Vec_Ortho

            Ortho = (1-(np.divide((np.linalg.norm(Vec_5)), np.linalg.norm(Vec_4))))*(100-(Iso+Hex+Tet))
            if np.linalg.norm(Vec_5) == 0:
                Mon = 0
                Tric = 0
            else:

                Mono_proj = np.zeros([21,21])
                for x in range(0, 9):
                    Mono_proj[x,x] = 1

                Mono_proj[11, 11] = 1
                Mono_proj[14,14] = 1
                Mono_proj[17, 17] = 1
                Mono_proj[20, 20] = 1

                Vec_Mono = np.matmul(Vec_5, Mono_proj)
                E_Mono = Vec_Mono * Mults
                Mono_tensor = [[E_Mono[0][0], E_Mono[0][5], E_Mono[0][4], E_Mono[0][9], E_Mono[0][13], E_Mono[0][17]],

                [E_Mono[0][5], E_Mono[0][1], E_Mono[0][3], E_Mono[0][15], E_Mono[0][10], E_Mono[0][14]],

                [E_Mono[0][4], E_Mono[0][3], E_Mono[0][2], E_Mono[0][12], E_Mono[0][16], E_Mono[0][11]],

                [E_Mono[0][9], E_Mono[0][15], E_Mono[0][12], E_Mono[0][6], E_Mono[0][20], E_Mono[0][19]],

                [E_Mono[0][13], E_Mono[0][10], E_Mono[0][16], E_Mono[0][20], E_Mono[0][7], E_Mono[0][18]],
                
                [E_Mono[0][17], E_Mono[0][14], E_Mono[0][11], E_Mono[0][19], E_Mono[0][18], E_Mono[0][8]]]
            
                Vec_6 = Vec_5 - Vec_Mono
                Mono = (1-(np.divide((np.linalg.norm(Vec_6)), np.linalg.norm(Vec_5))))*(100-(Iso+Hex+Tet+Ortho))

            Tric = 100 - (Iso+Hex+Tet+Ortho+Mono)
            Love_A = Iso_tensor[0, 0] + Hex_tensor[0][0]
            Love_C = Iso_tensor[2, 2] + Hex_tensor[2][2]
            Love_F = Iso_tensor[1, 2] + Hex_tensor[1][2]
            Love_L = Iso_tensor[3, 3] + Hex_tensor[3][3]
            Love_N = Iso_tensor[5, 5] + Hex_tensor[5][5]

            ACFLN = [Love_A, Love_C, Love_F, Love_L, Love_N]
            ACFLN = np.array(ACFLN)
            ACFLN = np.reshape(ACFLN, [1,5])

            Trix = np.array([[1, -1, 1, 0, 0], [1,1,1,0,0], [1, 0, -3,-2,-2], [0,0,0,1,1], [0,0,0,1,-1]])

            Backus = np.matmul(np.linalg.inv(Trix), ACFLN.conj().transpose())

            EtaK = (Love_F + Love_L)/(np.sqrt(Love_A-Love_L)*np.sqrt(Love_C-Love_L))
            Eta = Love_F/ (Love_A - (2*Love_L))
            BackusConj = Backus.conj().transpose()
            eigsTwoConj = eigsTwo.conj().transpose()
            bisectrixOneConj = Bisectrix[:,0].conj().transpose()
            bisectrixTwoConj = Bisectrix[:,1].conj().transpose()
            bisectrixThreeConj = Bisectrix[:,2].conj().transpose()

            Line = [Iso, Hex, Ortho, Tet, Mono, Tric, 0, 0,BackusConj[0][0], BackusConj[0][1], BackusConj[0][2], BackusConj[0][3], BackusConj[0][4], Love_A, Love_C, Love_F, Love_L, Love_N, Eta, EtaK, eigsTwoConj[0], eigsTwoConj[1], eigsTwoConj[2], bisectrixOneConj[0], bisectrixOneConj[1], bisectrixOneConj[2], bisectrixTwoConj[0], bisectrixTwoConj[1], bisectrixTwoConj[2], bisectrixThreeConj[0],bisectrixThreeConj[1], bisectrixThreeConj[2] ]
            Line = np.array(Line)
            Line = np.reshape(Line,[1, 32])
            Line = np.around(Line, 4)

            isoComponent = np.around(isoComponent, 4)
            hexagonalComponent = np.around(hexagonalComponent, 4)
            orthorhombicComponent = np.around(orthorhombicComponent, 4)
           


        

    #print(Line, '\n', isoComponent, '\n', hexagonalComponent, '\n', orthorhombicComponent)
    return Line, isoComponent, hexagonalComponent, orthorhombicComponent
    

def print_csv_table(databaseFrame = None, rowAddition=None):
    columnNames = ["Paper Reference", "tensor #", "Rock Type", "density (g/cm3)", "C11 (Mbar)", "C12", "C13","C14", "C15", "C16", "C22", "C23", "C24", "C25", "C26", "C33", "C34", "C35", "C36", "C44", "C45", "C46", "C55", "C56", "C66"]
    fileLoc = "databases/DB_Tensors_Rhos_OrigOrientation_CSV.csv"
    dbFile = pd.read_csv(fileLoc, names=columnNames, sep=',', engine='python')
    dbFrame = pd.DataFrame(data=dbFile, columns=columnNames)
    if rowAddition != None:
        databaseFrame = pd.DataFrame.from_dict(databaseFrame)
        databaseFrame = databaseFrame.append(rowAddition, ignore_index=True)
        #databaseFrame = databaseFrame.({"density (g/cm3)": float, "C11 (Mbar)": float})
        dbFrame = databaseFrame

    maxRows = len(dbFrame.index)
    
    dbDataTable = dash_table.DataTable(
        id='db_csv_table',
        columns=[{"name":i,"id": i, "deletable": True, "selectable": True} for i in dbFrame.columns],
        data=dbFrame.to_dict('records'),
        export_columns='visible',
        selected_rows=[],
        selected_columns=[],
        export_format='csv',
        export_headers='names',
        page_action='native',
        page_size=100,
        filter_action='native',
        row_selectable="multi",
        row_deletable=True,
        editable='true',
        style_cell_conditional=[{
            'if': {'column_id': 'Remark'},
            'textAlign': 'left'

        }],

        style_table={'height': '350px', 'overflowY': 'auto', 'text-align': 'left'},
    )
    
    #print(dbFrame)
    return dbDataTable, dbFrame.to_dict()



def print_breakdown(line=None):
    columnNames = ["Isotropic", "Hexagonal", "Orthorhombic", "Tetragonal", "Monoclinic", "Triclinic", "0", "0.", "Backus'.", "Backus'._", "Backus',", "Backus' ", "Backus'  ", "Love_A", "Love_C", "Love_F", "Love_L", "Love_N", "Eta", "EtaK", "Eigs2' ", "Eigs2'  ", "Eigs'2.", "Bisectrix[:,1]", "Bisectrix[:,1] ", "Bisectrix[:,1]  ", "Bisectrix[:,2]", "Bisectrix[:,2] ", "Bisectrix[:,2]  ", "Bisectrix{:,3]", "Bisectrix[:,3] ", "Bisectrix[,:3]  "]
    
    arrayCheck = np.any(line)
    if arrayCheck == False:
        return 1
    dbFrame = pd.DataFrame(data=line, columns=columnNames)

    maxRows = len(dbFrame.index)
    csvOutput = dbFrame.to_csv(index=False, encoding='utf-8')
    #print(dbFrame)

    dbDataTable = dash_table.DataTable(
        id='db_breakdown_table',
        columns=[{"name":i,"id": i, "deletable": True, "selectable": True} for i in dbFrame.columns],
        data=dbFrame.to_dict('records'),
        export_columns='visible',
        selected_rows=[],
        export_format='csv',
        export_headers='names',
        page_action='native',
        page_size=100,
        filter_action='native',
        row_selectable="multi",
        row_deletable=True,
        editable='true',
        style_cell_conditional=[{
            'if': {'column_id': 'Remark'},
            'textAlign': 'left'

        }],

        style_table={'height': '150px', 'overflowY': 'auto', 'text-align': 'left', 'overflowX':'scroll'},
    )

    return dbDataTable

def print_iso_decomp(iso=None):
    columnNames = ["test", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"]
    arrayCheck = np.any(iso)
    if arrayCheck == False:
        return 1
    dbFrame = pd.DataFrame(data=[iso], columns=columnNames)
    #print(dbFrame)

    maxRows = len(dbFrame.index)
    csvOutput = dbFrame.to_csv(index=False, encoding='utf-8')

    dbDataTable = dash_table.DataTable(
        id='isotropic_components_table',
        columns=[{"name":i,"id": i, "deletable": True, "selectable": True} for i in dbFrame.columns],
        data=dbFrame.to_dict('records'),
        export_columns='visible',
        selected_rows=[],
        export_format='csv',
        export_headers='names',
        page_action='native',
        page_size=100,
        filter_action='native',
        row_selectable="multi",
        row_deletable=True,
        editable='true',
        style_cell_conditional=[{
            'if': {'column_id': 'Remark'},
            'textAlign': 'left'

        }],

        style_table={'height': '150px', 'overflowY': 'auto', 'text-align': 'left', 'overflowX':'scroll'},
    )

    return dbDataTable


def print_hexag_decomp(hexag=None):
    columnNames = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"]
    arrayCheck = np.any(hexag)
    if arrayCheck == False:
        return 1
    dbFrame = pd.DataFrame(data=[hexag], columns=columnNames)

    maxRows = len(dbFrame.index)

    dbDataTable = dash_table.DataTable(
        id='hexagonal_components_table',
        columns=[{"name":i,"id": i, "deletable": True, "selectable": True} for i in dbFrame.columns],
        data=dbFrame.to_dict('records'),
        export_columns='visible',
        selected_rows=[],
        export_format='csv',
        export_headers='names',
        page_action='native',
        page_size=100,
        filter_action='native',
        row_selectable="multi",
        row_deletable=True,
        editable='true',
        style_cell_conditional=[{
            'if': {'column_id': 'Remark'},
            'textAlign': 'left'

        }],

        style_table={'height': '150px', 'overflowY': 'auto', 'text-align': 'left', 'overflowX':'scroll'},
    )
   

    return dbDataTable

def print_ortho_decomp(ortho=None):
    columnNames = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"]
    arrayCheck = np.any(ortho)
    if arrayCheck == False:
        return 1
    dbFrame = pd.DataFrame(data=[ortho], columns=columnNames)

    maxRows = len(dbFrame.index)

    dbDataTable = dash_table.DataTable(
        id='orthorhombic_components_table',
        columns=[{"name":i,"id": i, "deletable": True, "selectable": True} for i in dbFrame.columns],
        data=dbFrame.to_dict('records'),
        export_columns='visible',
        selected_rows=[],
        export_format='csv',
        export_headers='names',
        page_action='native',
        page_size=100,
        filter_action='native',
        row_selectable="multi",
        row_deletable=True,
        editable='true',
        style_cell_conditional=[{
            'if': {'column_id': 'Remark'},
            'textAlign': 'left'

        }],

        style_table={'height': '150px', 'overflowY': 'auto', 'text-align': 'left', 'overflowX':'scroll'},
    )
   

    return dbDataTable

def returnUserInputDataFrame(databaseFrame = None):
    
    databaseFrame = pd.DataFrame.from_dict(databaseFrame)
    dbFrame = databaseFrame
    return dbFrame

def averagingTensorDescription():
    description = "This is the tensor averaging tool. Enter the number of tensors you wish to average (between 2 and 5) in the left most input box, and then select the tensors in the following drop down boxes on the same line. Below, enter a weight for the given tensor above. The default is 1. Click 'Submit Average Tensor Options' to calculate the average tensor."
    return description

def calculateVelocitiesDescription():
    description = "Enter in a tensor whose velocity you want to calculate in the dropdown tab above. The tensor '1' is selected already by default. Then, hit the 'Submit Tensor' button to begin the calculation"
    return description



#tv_averaging(2,9,9)

#calculate_tensor_symmetries(1, plotType="RadialPlots", plotType3D=None)

#decomposition(1)
