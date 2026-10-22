#!/usr/bin/env python3

import argparse
import requests
import json
import re
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio

commits = defaultdict(lambda:0)
contributors = defaultdict(lambda:0)
pattern = re.compile('\\s+')

formats = ['html','svg','png']
rgb = dict(
           blue   = '4480cf',
           orange = 'e69644',
           green  = '7bcf44',
          )

def collect_commits(token,
                    reference,
                    N=1000,
                    params=dict()):
    """
    Query DAMASK GitLab instance for latest commit history.

    Parameters
    ----------
    token : string
        GitLab API access token.
    reference : string
        Commit reference to use as end point.
    N : int
        Number of latest commits to consider. Defaults to 1000.
    params : dict
        Additional parameters to pass to API call.

    Returns
    -------
    commits : json
        List of latest commits as JSON output from API call.

    """
    url_endpoint = 'https://git.damask-multiphysics.org/api/v4/projects/damask%2FDAMASK/repository/commits'
    params.update(
        all='false',
        ref_name=reference,
        per_page=100,
        access_token=token,
        )
    commits = []
    r = requests.Response()

    while len(commits) < N and r.headers.get('X-Next-Page') != '':
        params['page'] = r.headers.get('X-Next-Page',1)
        r = requests.get(url_endpoint,params=params)
        res = json.loads(r.text)
        if isinstance(res, dict) and 'error' in res.keys():
            raise RuntimeError(f'{res["error"]} error when loading commit: {res["error_description"]}')
        commits += res[:min(N-len(commits),params['per_page'])]
    return commits


def tag_dates(token,
              params=dict()):
    """
    Query DAMASK GitLab instance for creation date of all tags.

    Parameters
    ----------
    token : string
        GitLab API access token.
    params : dict
        Additional parameters to pass to API call.

    Returns
    -------
    dates : dict
        Creation date of each tag in %Y-%m-%d format.

    """
    url_endpoint = 'https://git.damask-multiphysics.org/api/v4/projects/damask%2FDAMASK/repository/tags'
    params.update(
        access_token=token,
        )
    response = requests.get(url_endpoint,params=params)
    return dict([(tag['name'],
                  datetime.fromisoformat(tag['commit']['created_at']).replace(tzinfo=timezone.utc).strftime('%Y-%m-%d'))
                  for tag in response.json()]) if response.status_code == 200 else None


parser = argparse.ArgumentParser(Path(__file__).name,
                                 description='Plot cumulative commit history of DAMASK GitLab repository',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('token', metavar='string', type=str,
                    help='GitLab API token')
parser.add_argument('--reference', metavar='string', default='HEAD', type=str,
                    help='commit reference to use as end point')
parser.add_argument('--commits', '-N', metavar='int', default=1000, type=int,
                    help='number of recent commits to consider')
parser.add_argument('--contributors', '-n', metavar='int', default=5, type=int,
                    help='number of top contributors to include')
parser.add_argument('--since', metavar='date or period',
                    help='starting date for commits')
parser.add_argument('--log', action='store_true',
                    help='use logarithmic scale')
parser.add_argument('--color', choices=rgb.keys(),
                    default=list(rgb)[0],
                    help='line color')
parser.add_argument('--template', choices=pio.templates.keys(),
                    default='ggplot2',
                    help='plotly template')
parser.add_argument('--outdir',
                    default='.',
                    help='output directory')
parser.add_argument('--format', choices=formats, nargs='*',
                    default=formats[:1],
                    help='output format(s)')

args = parser.parse_args()

params = {}

if args.since is not None:
    duration = 0
    for what,property in dict(years={'key':'y','factor':365},
                              months={'key':'m','factor':30},
                              weeks={'key':'w','factor':7},
                              days={'key':'d','factor':1},
                              ).items():
        m = re.search('(\\d+)'+property['key'],args.since.lower())
        duration += property['factor']*int(m.group(1)) if m else 0

    args.since = datetime.now(timezone.utc)-timedelta(days=duration) if duration > 0 else datetime.fromisoformat(args.since)
    params['since'] = args.since.strftime('%Y-%m-%d')

for c in collect_commits(args.token,args.reference,args.commits,params):
    date = datetime.fromisoformat(c['committed_date']).replace(tzinfo=timezone.utc).strftime('%Y-%m-%d')
    commits[date] += 1
    contributors[pattern.sub(' ',c['committer_name'].split('@')[0])] += 1

hst = pd.DataFrame({'date':commits.keys(),
                    'commits':commits.values()}).sort_values(by='date')
hst['cum'] = hst['commits'].cumsum()
commits_min = hst['cum'].iloc[0]
commits_max = hst['cum'].iloc[-1]
date_min = hst['date'].iloc[0]
date_max = hst['date'].iloc[-1]

dups = {}
contributors.pop('Testbot',None)

cnts = set(contributors)

for c in cnts:
    for co in cnts - set([c]):
        if      np.all([p in co.lower() for p in re.split('[ .]',c.lower())]) \
            and (c not in dups or len(dups[c]) < len(co)):
            dups[c] = co

for c,co in dups.items():
    contributors[co] += contributors.pop(c)

cnt = pd.DataFrame({'contributor':contributors.keys(),
                    'commits':contributors.values()}).sort_values(by='commits',ascending=False)

scale = (lambda x:np.log(x)) if args.log else (lambda x:x)

fig = px.line(hst,
              x='date',
              y='cum',
              hover_data=['commits'],
              range_y=[np.power(10,np.floor(np.log10(commits_min))) if args.log else 0,commits_max],
              log_y=args.log,
              line_shape='hv',
              labels={'cum':'cumulative commits'},
              title='DAMASK repository commits',
              color_discrete_sequence=[f'#{rgb[args.color]}'],
              template=args.template,
              )

for t,d in {k:v for k,v in tag_dates(args.token).items() if date_min <= v <= date_max}.items():
    fig.add_annotation(
        opacity=0.5,
        xref='x',axref='x',
        yref='y domain',ayref='y domain',
        xanchor='right',
        yanchor='middle',
        x =d,
        ax=d,
        y =1,
        ay=0.95,
        text=t,
        hovertext=d,
        textangle=-90,
        arrowhead=0,
        borderpad=10,
    )

for row in cnt.iloc[:args.contributors].itertuples():
    fig.add_annotation(
        opacity=0.5,
        xref='x domain',axref='x domain',
        yref='y',ayref='y',
        xanchor='right',
        yanchor='middle',
        x =1,
        ax=0.97,
        y =scale(row.commits),
        ay=scale(row.commits),
        text=row.contributor,
        hovertext=str(row.commits),
        arrowhead=0,
        borderpad=5,
    )

outdir = Path(args.outdir)
outdir.mkdir(parents=True,exist_ok=True)
outname = outdir/( 'commits'
                 +f'_{args.commits if args.commits == hst["cum"].iloc[-1] or args.since is None else hst["date"].iloc[0]}'
                 +f'_{args.contributors}'
                 +('_log' if args.log else ''))

for f in args.format:
    fig.write_html(outname.with_suffix('.html'),full_html=False,include_plotlyjs='cdn',include_mathjax='cdn') if f == 'html' else \
    fig.write_image(outname.with_suffix(f'.{f}'))
