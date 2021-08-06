import click
import time
import os
import re
import json

from regFilter.filter import *
from regFilter.spectra import *

def getSpectraIndexes(mgffile):
    spectrum = {}
    f = open(mgffile)
    lines = np.array(f.readlines())
    f.close()
    pos = [i for i, v in enumerate(lines) if 'SCANS=' in v]
    ind = [re.sub('[^0-9]', '', x) for x in lines[pos]]
    return ind

def matchFeats(query, ref, ppm=15, rtabs=20):
    mzdiff = ref['row m/z'].apply(lambda a: ((a-query['row m/z'])/query['row m/z'])*10**6).abs()
    rtdiff = (ref['row retention time']-query['row retention time']).abs()
    diff = (mzdiff < ppm) & (rtdiff < rtabs)
    ans = ref.loc[diff, ['row ID', 'row m/z', 'row retention time']]
    if len(ans):
        ans['q-id'] = query['row ID']
        ans['q-m/z'] = query['row m/z']
        ans['q-rt'] = query['row retention time']
        ans['ppm'] = mzdiff[diff]
        ans['rtabs'] = rtdiff[diff]
        return ans
    else:
        return pd.DataFrame()

@click.group()
def reg_filt():
    pass

@reg_filt.command()
@click.option("--feat",
              default='',
              help="MZmine feature table")
@click.option("--norm",
              default=1,
              help="Whether to normalize (TIC)")
@click.option("--score",
              default=0.9,
              help="Regression determination coefficient")
@click.option("--out",
              default='',
              help="Optional output file")
def filter(feat, norm, score, out):
    start = time.time()
    samp = Filter(feat)
    samp.formatFeatures(norm=norm)
    samp.filterFeatures(r2=score)
    end = time.time()
    print('Time:', end-start)

    if out!='':
        samp.feat_area[samp.sel[samp.sel==1].index.tolist()].to_csv(out, sep='\t')

@reg_filt.command()
@click.option("--quant",
              help="Average quant table")
@click.option("--idx",
              help="Feature index")
def plot(quant, idx):
    samp = Filter(quant)
    feat = pd.read_csv(quant, sep='\t', index_col='Groups')
    samp.feat_area = feat
    y = feat.index.str.replace('[^0-9]', '').str.replace('^$', '0').tolist()
    y = np.array([float(i) for i in y]).reshape(-1, 1)
    samp.y = y
    samp.filterReg(feat_idx=idx, plot=True)

@reg_filt.command()
@click.option("--queryq",
              help="Average quant query table")
@click.option("--queryf",
              help="Features query table")
@click.option("--refq",
              help="Average quant reference table")
@click.option("--reff",
              help="Features reference table")
@click.option("--mgf",
              help="Reference mgf file")
@click.option("--ppm",
              default=15,
              help="Part Per milion search difference")
@click.option("--rtabs",
              default=20,
              help="Absolute retention time search difference")
@click.option("--out",
              help="Output file name for match and stats")
def match_tables(queryq, queryf, refq, reff, mgf, ppm, rtabs, out):
    start = time.time()
    qfeat = pd.read_csv(queryf)
    qquant = pd.read_csv(queryq, sep='\t', index_col='Groups')
    rfeat = pd.read_csv(reff)
    rquant = pd.read_csv(refq, sep='\t', index_col='Groups')
    qlist = []

    for q in qquant.columns:
        qlist.append(matchFeats(qfeat.loc[int(q), ['row ID', 'row m/z', 'row retention time']],
                                rfeat, ppm=ppm, rtabs=rtabs))

    end = time.time()
    print('Time for comparison:', end-start)
    qtab = pd.concat(qlist).reset_index()
    ind = getSpectraIndexes(mgf)

    stats = {'N feat on query': qfeat.shape[0],
             'N feat on ref': rfeat.shape[0],
             'N feat on query quant': qquant.shape[1],
             'N feat on ref quant': rquant.shape[1],
             'N matching feats': len(qtab['q-id'].unique()),
             'N feat in spectra': len(ind),
             'Inter match spectra': len(set(ind).intersection(set(qtab['row ID'].astype(str))))
            }
    with open('%s_stats.json' % out) as f:
        json.dump(stats, f)

    qtab.to_csv('%s_match.tsv' % out, sep='\t', index=None)

@reg_filt.command()
@click.option("--feat",
              help="Feature table")
@click.option("--idx",
              help="Feature index")
@click.option("--path",
              help="Path to raw files")
@click.option("--ppm",
              default=15,
              help="Part Per milion search difference")
@click.option("--rtabs",
              default=20,
              help="Absolute retention time search difference")
@click.option("--tp",
              default='TI',
              help="Type, either TI or BP")
def plot_xic(feat, idx, path, ppm, rtabs, tp):
    samp = pd.read_csv(feat)

    fnames = samp.columns
    fnames = fnames[fnames.str.contains('Peak area')].str.replace(' Peak area', '').tolist()
    for fn in fnames:
        mzml = re.sub('.mzXML$', '.mzML', fn)
        if not os.path.exists(f'{path}/{mzml}'):
            mzXML2mzML(f'{path}/{fn}')

    cm = plt.get_cmap('gist_rainbow')

    cdict = {}
    for i in range(len(fnames)):
        cdict[re.sub('_MS1.*', '', fnames[i]).split('_')[2]] = i

    for k,v in cdict.items():
        cdict[k] = cm(v/3*3.0/len(cdict))

    for fn in fnames:
        leg = re.sub('_MS1.*', '', fn)
        cid = leg.split('_')[2]
        fn.replace('mzXML', 'mzML')
        exp = loadExp(f'{path}/{fn}')
        tic = getTIC(exp)
        mz = samp.loc[int(idx), 'row m/z']
        rt = samp.loc[int(idx), 'row retention time']
        xic = getXIC(mz, rt*60, tic, type=tp, rttol=rtabs,
                     ppm=ppm)
        plt.plot(xic['rt'], xic['xic'], color=cdict[cid], label=leg)
        exp.reset()

    plt.legend(fontsize=10)
    plt.show()

if __name__ == '__main__':
    reg_filt()
