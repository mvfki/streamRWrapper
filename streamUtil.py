'''In this script there are some customized functions for using stream to do 
single-cell analysis, though stream itself is already a function wrapper for
performing everything
Reference link for stream: https://github.com/pinellolab/STREAM
'''
import sys
import stream as st
import seaborn as sns
import numpy as np
from sklearn.cluster import KMeans
from collections import OrderedDict, defaultdict
import copy
import anndata

def rgb216(color):
    '''Mainly for converting seaborn RGB color code (from 0 to 1) to 
    hexadecimal, which stream plotting methods support.'''
    c = '#%02x%02x%02x' % (int(round(color[0] * 256)), 
                           int(round(color[1] * 256)), 
                           int(round(color[2] * 256)))
    return c

class streamSingleCellSamples(object):
    """
    A streamSamples object stores information read and processed by stream 
    package and wraps standard pipeline commands as method functions
    Values:
        adata - anndata.AnnData Object. All stream original function that takes
                an anndata.AnnData Object as an input argument can work by 
                running st.function(streamSingleCellSamples.adata, arg2, ...). 
    Initialization:
        When initializing, basic filtration will be done, including removing 
        duplicated obs and vars, removing mitochondria genes, and optional 
        normalization. 
        fileName       - str, required. 
                         The file name of the cellranger MTX output. Only 
                         support cellranger MTX format for now. The 
                         "barcodes.tsv" and "genes.tsv" file should be located 
                         at the exactly same path and named exactly as 
                         mentioned here. 
        cellLabel      - str, required. 
                         A TSV file with a label each line, and each line 
                         corresponds to a cell in "barcodes.tsv".
        cellLabelColor - str, required. 
                         A TSV file with each line written in "<label>\\t<hex-
                         color-code>" format.
        rawCount       - Boolean, optional (default: True). 
                         If True, st.normalize_per_cell() and 
                         st.log_transform() will be performed. Note that this 
                         normalization produces RPM values rather than z-score.
    """
    def __init__(self, fileName, cellLabel, cellLabelColor, rawCount = True):
        self.adata = st.read(file_name = fileName, file_format = 'mtx')
        st.add_cell_labels(self.adata, file_name = cellLabel)
        st.add_cell_colors(self.adata, file_name = cellLabelColor)
        self.adata.var_names_make_unique()
        self.adata.obs_names_make_unique()
        print('Raw input parsed...')
        print(self.adata)
        self.nCells = self.adata.n_obs
        self.nGenes = self.adata.n_vars
        self._keepCurrentRecords()
        st.remove_mt_genes(self.adata)
        if rawCount:
            st.normalize_per_cell(self.adata)
            st.log_transform(self.adata)
        self.backupDict = {}
        self.backupKey = 0
        self.backup(0)
        print('Initial backup saved with key: 0')
        print('Restore with self.restoreFromBackup()')

    def backup(self, key = None):
        '''
        Make a copy of the current anndata.AnnData Object. Can be restored by
        streamSingleCellSamples.restoreFromBackup(key).
        Arguments:
        key - Hashable value (e.g. str or int), optional.  
              self.backupDict will have keys to identify different backup 
              version. If not specified, it will be an auto incremented int.
        '''
        sys.stderr.write('WARNING: Maintaining backup takes memory. \n')
        if key == None:
            self.backupKey += 1
        tmp = copy.deepcopy(self.adata)
        self.backupDict[self.backupKey] = tmp
        if type(key) == str:
            key = '"' + key + '"'
        print('Current results backuped with key:', key)
        print('Restore with self.restoreFromBackup(%s)' % str(key))

    def restoreFromBackup(self, key = 0):
        self.adata = copy.deepcopy(self.backupDict[key])

    def preprocess(self, min_num_genes = None, min_num_cells = None, 
                   expr_cutoff = 1, loess_frac = 0.01, n_genes = None,  
                   nb_pct = 0.1, n_cluster = 15):
        '''
        Do the filtration, variable gene selection and dimension reduction to 
        the dataset. By default it also performs a KMeans clustering for 
        labeling small outlier groups. 
        Arguments:
        min_num_genes - int, optional (default: None). 
                        Minimum number of genes expressed. For filtering 
                        cells. 
        min_num_cells - int, optional (default: max(5, (nCells * 0.001)) ). 
                        Minimum number of cells expressing one gene. For 
                        filtering genes. 
        expr_cutoff   - float, optional (default: 1). 
                        Expression cutoff. If greater than expr_cutoff, the 
                        gene is considered 'expressed'. 
        loess_frac    - float, optional (default: 0.01). Between 0 and 1.
                        The fraction of the data used when estimating each 
                        y-value in LOWESS function. For selecting variable 
                        genes. 
        n_genes       - int, optional (default: None). 
                        Number of variable genes to use. Will be automatically 
                        detected by stream if not specified. 
        nb_pct        - float, optional (default: 0.1). 
                        The percentage neighbor cells. For doing MLLE 
                        dimension reduction. Smaller number such as 0.01 is 
                        recommended when dealing with large dataset with more 
                        than 10,000 cells.
        n_cluster     - int, optional (default: 15). 
                        The number of cluster when it by default run a KMeans 
                        clustering to help remove the small outlier groups.  
        '''
        st.filter_cells(self.adata, min_num_genes = min_num_genes)
        if min_num_cells == None:
            min_num_cells = int(round(self.adata.shape[0] * 0.001))
        st.filter_genes(self.adata, min_num_cells = min_num_cells, expr_cutoff = expr_cutoff)
        st.select_variable_genes(self.adata, loess_frac = loess_frac, n_genes = n_genes)
        st.dimension_reduction(self.adata, nb_pct = nb_pct)
        print('Initial UMAP visualization')
        st.plot_visualization_2D(self.adata, fig_legend_ncol = 4)
        self.kmeans = KMeans(n_cluster = n_cluster, init = 'k-means++')
        self.kmeans.fit(self.adata.obsm['X_vis_umap'])
        self._labelByCluster()
        print('Clustered UMAP visualization')
        st.plot_visualization_2D(self.adata, fig_legend_ncol = 4)
        self.backup(key = 'kmeans')
        print('If you want to remove some small groups of cells indicated by ')
        print('the clustering, Run streamSingleCellSamples.removeSmallCluster([\'cN\', ...])')
        print('And then run streamSingleCellSamples.confirmRemoval()')

    def removeSmallCluster(self, clusters):
        '''
        Remove the specified groups of cells from the dataset. Since the 
        coloring can be confusing sometimes, you can run this function 
        multiple times until you want to confirm.
        Argument:
        clusters - a str or a list of str, Required.
                   The str(s) is/are the cluster label(s) which is/are shown  
                   in the clustered UMAP visualization.
        '''
        self.restoreFromBackup('kmeans')
        cellToRemove = []
        if type(clusters) == str:
            cellToRemove = list(self.adata.obs['label'][self.adata.obs['label'] == clusters].index)
        elif type(clusters) == list:
            for c in clusters:
                cellToRemove.extend(list(self.adata.obs['label'][self.adata.obs['label'] == c].index))
        cellToRemove = set(cellToRemove)
        remain = set(self.adata.obs.index).difference(cellToRemove)
        Idx = []
        for i in remain:
            Idx.append(list(self.adata.obs.index).index(i))
        self.adata._inplace_subset_obs(Idx)
        print('Filtered and clustered UMAP visualization')
        st.plot_visualization_2D(self.adata, fig_legend_ncol = 4)
        print('If you are good with the result, run streamSingleCellSamples.confirmRemoval()')
        
    def confirmRemoval(self):
        '''
        Run this only when you are confirmed with the remove-by-clustering. 
        '''
        self.recoverRecords()
        print('Confirmed UMAP visualization')
        st.plot_visualization_2D(self.adata, fig_legend_ncol = 4)

    def _labelByCluster(self):
        nC = self.kmeans.n_clusters
        palettes = sns.color_palette("husl", nC)
        palettes = [rgb216(c) for c in palettes]
        self.adata.obs['label'] = 'unknown'
        self.adata.obs['label_color'] = 'gray'
        self.adata.uns['label_color'] = {}
        self.adata.uns['label_color']['unknown'] = 'gray'
        for i in range(nC):
            c0 = np.where(self.kmeans.labels_ == i)[0]
            label = 'c' + str(i)
            self.adata.obs['label'][c0] = label
            self.adata.obs['label_color'][c0] = palettes[i]
            self.adata.uns['label_color'][label] = palettes[i]

    def _keepCurrentRecords(self):
        self.allCells = self.adata.obs.index.to_list()
        self.label2cells = defaultdict(set)
        self.label2color = {}
        for cell, info in self.adata.obs.iterrows():
            self.label2cells[info['label']].add(cell)
            if info['label'] in self.label2color:
                if self.label2color[info['label']] != info['label_color']:
                    sys.stderr.write('In cell %s - label %s - color %s, does not match with %s. Probably it is not a unique cell. \n' % (cell, info['label'], info['label_color'], self.label2color[info['label']]))
                    sys.stderr.write('Ingnoring cell %s. \n' % (cell))
            else:
                self.label2color[info['label']] = info['label_color']

    def _recoverRecords(self):
        self.adata.uns['label_color'] = self.label2color
        existingCells = set(self.adata.obs.index.to_list())
        for label, cellSet in self.label2cells.items():
            cellSet = cellSet.intersection(existingCells)
            for cell in cellSet:
                self.adata.obs.at[cell, 'label'] = label
                self.adata.obs.at[cell, 'label_color'] = self.label2color[label]
    
    def __repr__(self):
        return self.adata

class testClass(object):
    """docstring for testClass"""
    def __init__(self, arg):
        super(testClass, self).__init__()
        self.arg = arg

    def func1(self):
        return self.arg

    def func2(self, arg2):
        self.arg += arg2
        return self.arg

def testFunc(integer):
    integer = int(integer)
    return (type(integer), integer)