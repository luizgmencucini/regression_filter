import click
import time

from regFilter.filter import *
from regFilter.spectra import *

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
    samp.filterReg(feat_idx=int(idx))

if __name__ == '__main__':
    reg_filt()
