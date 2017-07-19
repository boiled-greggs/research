#import matplotlib.pyplot as plt
#from matplotlib.figure import Figure
#from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
import pylab as plot
import plotly
import plotly.plotly as plt
import plotly.graph_objs as go
plotly.tools.set_credentials_file(username='boiled-greggs', api_key='vLwVsrAXGmlgJzki3Q4a')
bandf = input('enter name of bands file: ')
structure = input('enter name of output you\'d like: ')
f = open(bandf+'.txt')

kpt = []
energy = []
bandk = []
bande = []
fermi = float(input('enter fermi level: ')) #4.8252

for line in f:
    line = line.split()
    if len(line) == 2:
        line[0] = float(line[0])
        line[1] = float(line[1]) - fermi
        bandk.append(line[0])
        bande.append(line[1])
    else:
        kpt.append(bandk)
        energy.append(bande)
        bandk = []
        bande = []
data = []
for band in range(len(energy)):
    trace = go.Scatter(
        x = kpt[band],
        y = energy[band],
        line = dict(width=3,)
    )
    data.append(trace)
vlines = list()
node_x = [0.00000, 1.0000, 1.7071, 2.2071, 3.0731, 3.7802]
for i in range(len(node_x)):
    vlines.append(
        go.Scatter(
            x = [node_x[i], node_x[i]],
            y = [-100, 100],
            mode = 'lines',
            line = go.Line(color='#111111', width=1),
            showlegend=False
        )
    )
bandxaxis = go.XAxis(
    title = 'k path',
    range = [0, kpt[0][-1]],
    showgrid = True,
    showline = True,
    ticks = '',
    showticklabels = True,
    mirror = True,
    linewidth = 2,
    ticktext = [r'$\Gamma$', r'$H$', r'$N$', r'$P$', r'$\Gamma$', r'$N$'],
    tickvals = node_x
    )
bandyaxis = go.YAxis(
    title = '$E - E_f \quad (\\text{eV})$',
    range = [-6, 6],
    showgrid = True,
    showline = True,
    zeroline = False,
    mirror = 'ticks',
    ticks = 'inside',
    linewidth = 2,
    tickwidth = 2,
    )
layout = dict(title = 'P80 bands',
    width = 500,
    height = 700,
    xaxis = bandxaxis,
    yaxis = bandyaxis 
    )
fig = go.Figure(data = data + vlines, layout = layout) 
plt.iplot(fig, filename=structure+' bands')

fig, ax = plot.subplots(figsize=(3,4))

nodes = [0.0, 1.0000, 1.7071, 2.2071, 3.0731, 3.7802]
label = (r'$\Gamma$', r'$H$', r'$N$', r'$P$', r'$\Gamma$', r'$N$')

ax.set_xticks(nodes)
ax.set_xticklabels(label)

for n in range(len(nodes)):
    ax.axvline(x=nodes[n], linewidth=.3, color='black')

ax.set_title('P80 band structure')
ax.set_xlabel('Path in k-space')
ax.set_ylabel(r'$E$'+'–'+r'$E_f$')

ax.set_ylim([-6, 6])

for i in range(len(energy)):
    for j in range(len(energy[i])):
        energy[i][j] = energy[i][j] - fermi

for band in range(len(energy)):
    ax.plot(kpt[band], energy[band], linewidth=.5)

fig.savefig(structure+'.eps')
