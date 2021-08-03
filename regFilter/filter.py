import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

class Filter:
    def __init__(self, features):
        """ Filter by regression
        Returns
        -------
        feature  :  feature file or table
        """
        self.features = features


    def formatFeatures(self, col_pattern='Peak area', meta=None,
                       norm=True):
        if type(self.features)==str:
            feat_tab = pd.read_csv(self.features)
        else:
            feat_tab = self.features

        feat_area = feat_tab.copy()[feat_tab.columns[feat_tab.columns.str.contains(col_pattern)]]
        if meta is None:
            meta = feat_area.columns.str.replace('.mzXML Peak area', '').str.split('_').tolist()
            # assumes DILX dilution factor is on position 3
            meta = [x[2] for x in meta]

        # Normalize samples
        # scale?
        if norm:
            feat_area = feat_area.apply(lambda a: a/a.sum())

        # move features to columns and adds metadata
        feat_area = feat_area.T
        feat_area.reset_index(inplace=True, drop=True)
        feat_area['Groups'] = meta

        # Average technical replicates
        feat_area = feat_area.groupby('Groups').mean()

        # Assumes number describing increasing or decreacing ibdex for dilution
        # here ORI is used for original sample
        y = feat_area.index.str.replace('[^0-9]', '').str.replace('^$', '0').tolist()
        y = np.array([float(i) for i in y]).reshape(-1, 1)

        self.feat_tab = feat_tab
        self.feat_area = feat_area
        self.y = y


    def filterReg(self, feat_idx, r2=0.9, sig='-',
                  plot=False):
        X = self.feat_area[feat_idx].to_numpy()
        y = self.y

        X = X.reshape(-1, 1)
        reg = LinearRegression()
        reg.fit(X, y)
        sc = reg.score(X, y)
        b1 = reg.coef_[0][0]

        if plot:
            plt.scatter(y, X)

            plt.plot(reg.predict(X), X)

            plt.ylabel('Normalized Area')
            plt.xlabel('Dilution')

            plt.xticks([0, 1, 2, 3, 4], ['0', '1', '2', '3', '4'])
            plt.text(y[2,0], X[1,0], r"$R^2=%s$" % np.round(sc, 3))

        if sig=='-':
            if sc > r2 and b1 < 0:
                return 1
            else:
                return 0
        else:
            if sc > r2 and b1 > 0:
                return 1
            else:
                return 0


    def filterFeatures(self, r2=0.9, sig='-'):
        sel = self.feat_area.apply(lambda a: self.filterReg(a.name, r2=r2, sig=sig))

        self.sel = sel

        print('Numeber of fearures:', len(sel))
        print('Numeber of fearures after filtering:', sel.sum())

