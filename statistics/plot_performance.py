#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio

formats = ['html','svg','png']
rgb = dict(
           blue   = '4480cf',
           orange = 'e69644',
           green  = '7bcf44',
          )

def round_log(x):
    O,o = np.divmod(np.log10(x),1.0)
    r = 1+np.arange(10)
    return np.power(10,O)*r[np.argmin(np.abs(np.log10(r)-o))]


def floor_log(x):
    O,o = np.divmod(np.log10(x),1.0)
    r = 1+np.arange(10)
    return np.power(10,O)*r[np.searchsorted(np.log10(r),o,'left')-1]


def ceil_log(x):
    O,o = np.divmod(np.log10(x),1.0)
    r = 1+np.arange(10)
    return np.power(10,O)*r[np.searchsorted(np.log10(r),o,'right')]


parser = argparse.ArgumentParser(Path(__file__).name,
                                 description='Plot historical runtime and memory footprint of DAMASK Grid solvers',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--commits', '-N', metavar='int', type=int,
                    help='number of recent commits to consider')
parser.add_argument('--since', metavar='date or period',
                    help='starting date for commits')
parser.add_argument('--linear', action='store_true',
                    help='use linear scale')
parser.add_argument('--template', choices=pio.templates.keys(),
                    default='ggplot2',
                    help='plotly template')
parser.add_argument('--outdir',
                    default='.',
                    help='output directory')
parser.add_argument('--format', choices=formats, nargs='*',
                    default=formats[:1],
                    help='output format(s)')
parser.add_argument('filename',nargs='?',
                    default='./performance.txt',
                    help='file containing performance data')

args = parser.parse_args()

if args.since is None:
    since = 0.0
else:
    duration = 0
    for what,property in dict(years={'key':'y','factor':365},
                              months={'key':'m','factor':30},
                              weeks={'key':'w','factor':7},
                              days={'key':'d','factor':1},
                              ).items():
        m = re.search('(\\d+)'+property['key'],args.since.lower())
        duration += property['factor']*int(m.group(1)) if m else 0

    since = (datetime.now(timezone.utc)-timedelta(days=duration) if duration > 0 else
             datetime.fromisoformat(args.since)).timestamp()

df = pd.read_csv(args.filename,sep='\t',index_col=None)

tags,idx = np.unique(df['tag'],return_index=True)
commits = tags[np.argsort(idx)] if args.commits is None else tags[np.argsort(idx)][-args.commits:]
df = df[(df['epoch'] > since) & (df['tag'].isin(commits))]

df['date'] = pd.to_datetime(df['epoch'],unit='s')
df['memory/GiB'] = df['memory/kB']/2**20
df['processor'] = df['processor'].replace(to_replace=r'\s+', value=' ', regex=True)
compiler = df['compiler'].str.extract(
    r'^(?P<vendor>\S+).*?\s+(?P<version>\d+\.[\d.-]+)'
)
df['compiler'] = compiler['vendor']+' '+compiler['version']
df['tag'] = df['tag'].str[:8]
df.rename(columns={'tag':'commit',
                   'OMP_NUM_THREADS':'OMP threads',
                   'processor':'CPU',
                   'N_proc':'CPU threads'},
          inplace=True)

logmargin = 1.5
marker = dict(size=4, opacity=0.5)
legend = dict(
    title='',
    orientation='h',
    yanchor='bottom',
    y=1.02,
    xanchor='right',
    x=1,
    itemsizing='constant',
)
params = dict(
    x='date',
    color='solver',
    hover_data={'date':'|%-b %d %Y',
                'commit':True,
                'OMP threads':True,
                'compiler':True,
                'PETSc':True,
                'OS':True,
                'CPU':True,
    },
    log_y=not args.linear,
    color_discrete_sequence=[f'#{rgb["blue"]}',f'#{rgb["orange"]}',f'#{rgb["green"]}',],
    template=args.template,
)

me = 'wall time/s'
params.update(dict(y=me,
                   range_y=[0,ceil_log(df[me].max())] if args.linear else
                           [floor_log(df[me].min()/logmargin),
                             ceil_log(df[me].max()*logmargin)],
                   labels={me:'Elapsed real time / s'},
                   title='DAMASK performance',
                   ))
perf = px.scatter(df,**params)
perf.update_traces(marker=marker)
perf.update_layout(legend=legend)

me = 'memory/GiB'
params.update(dict(y=me,
                   range_y=[0,ceil_log(df[me].max())] if args.linear else
                           [floor_log(df[me].min()/logmargin),
                             ceil_log(df[me].max()*logmargin)],
                   labels={me:'Max. resident set size / GiB'},
                   title='DAMASK memory footprint',
                   ))
mem = px.scatter(df,**params)
mem.update_traces(marker=marker)
mem.update_layout(legend=legend)

outdir = Path(args.outdir)
outdir.mkdir(parents=True,exist_ok=True)

specifier = ( ('_' if args.commits is not None or args.since is not None else '')
             +( str(args.commits) if df["commit"].iloc[0] == commits[0] and args.commits is not None
             else datetime.fromtimestamp(since).strftime("%Y-%m-%d") if args.since is not None
             else '')
             +('_linear' if args.linear else '')
)
perfname = outdir/('runtime'+specifier)
memname = outdir/('memory'+specifier)

for f in args.format:
    for fig,outname in zip([perf,mem],[perfname,memname]):
        if f == 'html':
            fig.write_html(outname.with_suffix('.html'),
                           full_html=False,
                           include_plotlyjs='cdn',
                           include_mathjax='cdn')
        else:
            fig.write_image(outname.with_suffix(f'.{f}'))
