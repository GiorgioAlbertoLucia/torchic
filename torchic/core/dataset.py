import pandas as pd
import uproot
from ROOT import TH1F, TH2F

from torchic.core.histogram import AxisSpec, build_TH1, build_TH2
from torchic.utils.overload import overload, signature
from torchic.utils.terminal_colors import TerminalColors as tc

class SubsetDict:
    '''
        A dictionary to access DataFrame subsets
    '''

    def __init__(self):

        self._subsets = {}
        self.add_subset('full', lambda: slice(None))

    def add_subset(self, name, condition):

        if name in self._subsets.keys():
            raise ValueError(tc.RED+'[ERROR]: '+tc.RESET+f'Subset {name} already exists')
        self._subsets[name] = condition

    def __getitem__(self, key):
        return self._subsets[key]()

class Dataset:

    def __init__(self, data, **kwargs):
        '''
            Constructor for the Dataset class.
            
            Args:
                data (str, list, or pd.DataFrame): The input data to be loaded. If a string, it should be the path to a single file. If a list, it should be a list of paths to multiple files. If a pd.DataFrame, it should be the data itself.
                **kwargs: Additional keyword arguments to be passed to the pandas read_csv or read_parquet functions.
                    - columns (list): The list of columns to read from the file.
                    - folder_name (str): The name of the folder in the root file.
                    - tree_name (str): The name of the tree in the root file.
        '''
        
        self._data = None
        self._open(data, **kwargs)
        self._subsets = SubsetDict()

    #def __getitem__(self, key):
    #    return Dataset(self._subsets[key])
        
    # ADJUST TO HANDLE SUBSETS AS WELL
        
    def __getitem__(self, key):
        if ':' in key:
            key1, key2 = key.split(':')
            subset = self._subsets[key1]
            return subset[key2] if key2 else subset
        return self._data[key]
    
    def __setitem__(self, key, value):
        self._data[key] = value

    def _open(self, data, **kwargs):
        
        self._data = pd.DataFrame()

        if isinstance(data, pd.DataFrame):
            self._data = data
        elif isinstance(data, str) or (isinstance(data, list) and all(isinstance(file, str) for file in data)):
            self._files = data if isinstance(data, list) else [data]
            for file in self._files:
                if file.endswith('.root'):
                    self._data = pd.concat([self._data, self._open_root(file, **kwargs)], ignore_index=True, copy=False)
                elif file.endswith('.csv'):
                    print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+file+tc.RESET)
                    self._data = pd.concat([self._data, pd.read_csv(file, **kwargs)], ignore_index=True, copy=False)
                elif file.endswith('.parquet'):
                    print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+file+tc.RESET)
                    self._data = pd.concat([self._data, pd.read_parquet(file, **kwargs)], ignore_index=True, copy=False)
                else:
                    raise ValueError(tc.RED+'[ERROR]: '+tc.RESET+'Input data must be a list of .root or .csv files.')
        else:
            raise ValueError(tc.RED+'[ERROR]: '+tc.RESET+'Input data must be a string, a list of strings, or a pandas DataFrame.')
        
    def _open_root(self, file, **kwargs) -> pd.DataFrame:
        
        tree_name = kwargs.get('tree_name', None)
        uproot_kwargs = {key: value for key, value in kwargs.items() if key not in ['tree_name', 'columns', 'folder_name']}
        if tree_name is None:  
            print(tc.RED+'[ERROR]: '+tc.RESET+'tree_name must be specified when using a .root file.')
            return

        columns = kwargs.get('columns', None)
        folder_name = kwargs.get('folder_name', None)
        if folder_name is None:
            return uproot.open(f'{file}:{tree_name}').arrays(filter_name=columns, library='pd', **uproot_kwargs)

        if folder_name[-1] != '*':
            print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+f'{file}:{folder_name}/{tree_name}'+tc.RESET)
            return uproot.open(f'{file}:{folder_name}/{tree_name}').arrays(filter_name=columns, library='pd', **uproot_kwargs)

        file_folders = uproot.open(file).keys()
        tree_path_list = []
        for folder in file_folders:
            if folder_name[:-1] in folder and tree_name in folder:
                tree_path_list.append(folder)

        tmp_data = pd.DataFrame()
        for tree_path in tree_path_list:
            print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+f'{file}:{tree_path}'+tc.RESET)
            tmp_data = pd.concat([tmp_data, uproot.open(f'{file}:{tree_path}').arrays(filter_name=columns, library='pd', **uproot_kwargs)], ignore_index=True, copy=False)
        return tmp_data
    
    @property
    def data(self):
        return self._data
    
    @property
    def subsets(self):
        return self._subsets
    
    def add_subset(self, name, condition):
        self._subsets.add_subset(name, lambda: self._data[condition])

    def query(self, expr: str, *, inplace: bool = True, **kwargs) -> pd.DataFrame | None:
        '''
            Query the dataset using a string expression.
            
            Args:
                expr (str): The query expression.
                inplace (bool): Whether to modify the dataset in place.
                **kwargs: Additional keyword arguments to be passed to the pandas query function.
        '''
        
        if inplace:
            self._data.query(expr, inplace=True, **kwargs)
        else:
            tmp_data = self._data.query(expr, inplace=False, **kwargs).copy()
            return Dataset(tmp_data)

    def concat(self, other: 'Dataset', **kwargs) -> 'Dataset':
        '''
            Concatenate two datasets.

            Args:
                other (Dataset): The other dataset to concatenate.
                **kwargs: Additional keyword arguments to be passed to the pandas concat function.
        '''
        
        return Dataset(pd.concat([self._data, other.data], **kwargs))

    def describe(self, **kwargs) -> pd.DataFrame:
        '''
            Generate descriptive statistics of the dataset.

            Args:
                **kwargs: Additional keyword arguments to be passed to the pandas describe function.
        '''
        
        return self._data.describe(**kwargs)
    
    def eval(self, expr: str, **kwargs) -> pd.DataFrame:
        '''
            Evaluate an expression in the dataset.

            Args:
                expr (str): The expression to evaluate.
                **kwargs: Additional keyword arguments to be passed to the pandas eval function.
        '''
        
        return self._data.eval(expr, **kwargs)
    
    def apply(self, func, **kwargs) -> pd.DataFrame:
        '''
            Apply a function to the dataset.

            Args:
                func (function): The function to apply.
                **kwargs: Additional keyword arguments to be passed to the pandas apply function.
        '''
        
        return self._data.apply(func, **kwargs)
    
    def head(self, n: int = 5) -> pd.DataFrame:
        '''
            Return the first n rows of the dataset.

            Args:
                n (int): The number of rows to return.
        '''
        
        return self._data.head(n)

    @overload
    @signature('str', 'AxisSpec')  
    def build_hist(self, column: str, axis_spec_x: AxisSpec, **kwargs) -> TH1F:
        '''
            Build a histogram with one axis

            Args:
                column (str): The column to be histogrammed
                axis_spec_x (AxisSpec): The specification for the x-axis

            Returns:
                TH1F: The histogram
        '''
        subset = kwargs.get('subset', None)
        if subset:
            return build_TH1(self._subsets[subset][column], axis_spec_x)
        else:
            return build_TH1(self._data[column], axis_spec_x)
    
    @build_hist.overload
    @signature('str', 'str', 'AxisSpec', 'AxisSpec')
    def build_hist(self, column_x: str, column_y: str, axis_spec_x: AxisSpec, axis_spec_y: AxisSpec, **kwargs) -> TH2F:
        '''
            Build a histogram with two axes

            Args:
                column_x (str): The column to be histogrammed on the x-axis
                column_y (str): The column to be histogrammed on the y-axis
                axis_spec_x (AxisSpec): The specification for the x-axis
                axis_spec_y (AxisSpec): The specification for the y-axis

            Returns:
                TH2F: The histogram
        '''
        subset = kwargs.get('subset', None)
        if subset:
            return build_TH2(self._subsets[subset][column_x], self._subsets[subset][column_y], axis_spec_x, axis_spec_y)
        else:
            return build_TH2(self._data[column_x], self._data[column_y], axis_spec_x, axis_spec_y)
    
    