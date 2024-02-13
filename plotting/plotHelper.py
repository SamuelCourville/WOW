#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator
import numpy as np
import sys
import os
import argparse

#script_dir = os.path.abspath(os.path.dirname(__file__)) + '/'
script_dir="/Users/samuelcourville/Documents/JPL/combinedModel/plotting/"

def gmt_colormap( cptfile=script_dir + 'JC_colormap.cpt' ):
    df = pd.read_csv(cptfile,
        delim_whitespace=True, #separator is whitespace
        header=None, #no header
        comment='#',
        # nrows=1e6,
        # usecols=[3, 4, 6], #parse only 3,4,6 columns
        names=['begin','b_red','b_green','b_blue', 'end', 'e_red','e_green','e_blue' ], #set columns names
        )

    first = df.begin.iloc[0]
    last = df.end.iloc[-1]
    span = last - first

    x = ( (df.begin - first) / span ).to_list()
    x.append( ( df.end.iloc[-1] - first ) / span )
    # print( x )

    cdict = dict( red=[], green=[], blue=[] )

    for c in cdict.keys():
        yleft = df[ 'b_' + c ].iloc[0] / 256
        yright = df[ 'b_' + c ].iloc[0] / 256
        cdict[c].append( [ x[0], yleft, yright ] )
        for i in range(1, len(x) - 1):
            yleft = df[ 'e_' + c ].iloc[i - 1] / 256
            yright = df[ 'b_' + c ].iloc[i] / 256
            cdict[c].append( [ x[i], yleft, yright ] )
        yleft = df[ 'e_' + c ].iloc[-1] / 256
        yright = df[ 'e_' + c ].iloc[-1] / 256
        cdict[c].append( [ x[-1], yleft, yright ] )

    newcmp = LinearSegmentedColormap('GMT_pal', segmentdata=cdict, N=11256)

    # for c in cdict.keys():
    #   print( c, cdict[c][0], cdict[c][1], '...', cdict[c][-2], cdict[c][-1] )
    # sys.exit()

    return first, last, newcmp


def load_data( infile='gmt.txt' ):

    df = pd.read_csv(infile,
        delim_whitespace=True, #separator is whitespace
        header=None, #no header
        # nrows=1e6,
        names=['time_my','radius_km','K'], #set columns names
        )

    # print( df )

    if True:
        # correct first samples ( given in meters while the rest is in km )
        print( max(df['radius_km']) )
        df['radius_km'] = df.apply( lambda r: r['radius_km'] / 1e3 if r['time_my'] == 0 else r['radius_km'], axis=1 )
        print( max(df['radius_km']) )
    else:  # the fix below is bad - don't use
        sys.exit()
        # df.drop(df[df.radius_km > 10000].index, inplace=True )
        df = df[df.radius_km.duplicated(keep=False)]  # remove rows with a radius_km value that show up only once (bad inout data)
        # df.to_hdf( 'gmt.hd5', key='df' )



    df['radius_km'] = round( df['radius_km'], 1 )  # rounding to nearest 100m, in order to keep the number of rows manageable

    # print( df )
    # sys.exit()


    if True:
        print( 'scatterplot' )

        try:
            vmin, vmax, new_colormap = gmt_colormap()
        except:
            vmin = 0
            vmax = 300
            new_colormap = 'jet'


        plt.clf()
        # ax = df.plot.scatter( x='time_my', y='radius_km', c='K', colormap='viridis' )
        ax = df.plot.scatter( x='time_my', y='radius_km', c='K', s=1, colormap=new_colormap, vmin=vmin, vmax=vmax )
        plt.savefig( 'gmt_datapoint.png' )

        # sys.exit()

        plt.clf()
        # ax = df.plot.scatter( x='time_my', y='radius_km', c='K', colormap='viridis' )
        ax = df.plot.scatter( x='time_my', y='radius_km', c='K', s=1, colormap=new_colormap, vmin=vmin, vmax=vmax )
        ax.set_xscale('log')
        plt.xlim( max( 1, min(df['time_my']) ), max(df['time_my']) )
        plt.savefig( 'gmt_datapoint_log.png' )

        # sys.exit()

    print( df )
    # sys.exit()

    pivoted = df.pivot( index='radius_km', columns='time_my', values='K')

    # pivoted.to_hdf( 'gmt_pivoted.hd5', key='df' )
    # print( pivoted )

    # change = pivoted.diff(axis=1)
    # print( change )

    # change = pivoted.round(0).diff(axis=1)
    # print( change )

    # print( np.sqrt(np.square(change).sum(axis=0)) )

    remove_those_colums = []
    log_bin_count = 1000
    log_bin_multiplier = np.power( 10, np.log10( pivoted.columns[-1] ) / log_bin_count )
    for i in range(len(pivoted.columns)):
        # print( i, pivoted.columns[i] )
        if i == 0:
            last_good_i = i
            continue
        if pivoted.columns[i] / pivoted.columns[last_good_i] > log_bin_multiplier:
            last_good_i = i
            continue
        remove_those_colums.append( pivoted.columns[i] )
    pivoted = pivoted.drop( remove_those_colums,  axis='columns' )

    return pivoted


def make_plot( pivoted, outfile='gmt.png', title='EUROPA', cptfile=script_dir + 'GMT_pal.cpt' ):

    zi = pivoted.to_numpy()
    xi = pivoted.columns.to_numpy()
    yi = pivoted.index.to_numpy()

    try:
        vmin, vmax, new_colormap = gmt_colormap( cptfile=cptfile )
    except:
        vmin = 0
        vmax = 300
        new_colormap = 'jet'

    print( 'plotting' )
    plt.clf()
    plt.figure( figsize=(8, 6), dpi=200 ) # figsize in inches
    plt.title( title, fontsize=14, fontweight='bold' )

    # major contour lines
    major_levels = range(100, 2300, 100)
    CS = plt.contour( xi,yi,zi, levels=major_levels, linewidths=0.5, colors='black', linestyles='dashed' )
    
    # labels on a subset of the major contour lines
    # labeled_levels = [ x for x in major_levels if x <= 600 ]
    labeled_levels = [ x for x in major_levels ]
    clabels = plt.clabel( CS, labeled_levels, fmt='%.0f', fontsize=6 )
    for label in clabels:
        label.set_rotation(-90)
    # [ txt.set_bbox(dict(facecolor='white', edgecolor='none', alpha=0.65, boxstyle='round')) for txt in clabels ]

    # minor contour lines
    minor_levels = [ x for x in range(100, 1800, 25) if x not in major_levels ]
    CS2 = plt.contour( xi,yi,zi, levels=minor_levels, linewidths=0.5, colors='white', linestyles='dashed' )

    # colored plot and colorbar
    CS3 = plt.contourf( xi,yi,zi, levels=500, cmap=new_colormap, vmin=vmin, vmax=vmax )  # 20220423: levels was 1000
    cbar = plt.colorbar(CS3, ticks=major_levels, orientation='horizontal') # draw colorbar
    cbar.ax.set_xlabel('Temperature (K)')

    # x-axis setiings
    plt.gca().set_xscale('log')
    plt.xlabel( 'Time (My)', fontsize=14, fontweight='bold' )
    plt.xlim( 10, max(xi) )


    # y-axis setiings
    plt.ylabel( 'Radius (km)', fontsize=14, fontweight='bold' )
    plt.gca().yaxis.set_minor_locator( MultipleLocator(25) )  # minor ticks every N km
    plt.gca().yaxis.set_major_locator( MultipleLocator(100) )  # major ticks every N km

    plt.savefig( outfile, bbox_inches='tight' )
    plt.close()



#def get_args():
#    parser = argparse.ArgumentParser(description='Demo', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#    parser.add_argument('--datafile', '-d', default='gmt.txt', help='file with 3 columns [ time_my, radius_km, K ]' )
#    parser.add_argument('--outfile', '-o', default='gmt.png', help='output image' )
#    parser.add_argument('--title', '-t', default='EUROPA', help='title of the plot' )
#
#    return parser.parse_args()


#def main():

#    args = get_args()

#    pivoted = load_data( infile=args.datafile )

#    # print( pivoted )

#    if False:
#        print( list(pivoted.index) )
#        sys.exit()

#    if False:
#        print( pivoted.iloc[-15:-1, 30:40] )  # indices, columns
#        sys.exit()

#    if True:
#        pivoted_interpolated = pivoted.interpolate( method='slinear' )


#        # print( pivoted.iloc[1000:1010, 60:70] )
#        # print( pivoted_interpolated.iloc[1000:1010, 60:70] )
#        # sys.exit()
#        pivoted = pivoted_interpolated



#    make_plot( pivoted, outfile=args.outfile, title=args.title )


#if __name__ == '__main__':
#    main()
