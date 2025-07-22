import pandas as pd
import uproot
import boost_histogram as bh
from ROOT import TH1F, TH2F

from torchic.core.histogram import AxisSpec, build_TH1, build_TH2, build_boost1, build_boost2
from torchic.utils.terminal_colors import TerminalColors as tc

class SubsetDict:
    '''
        A dictionary to access DataFrame subsets
    '''

    def __init__(self):

        self._subsets = {}

    def add_subset(self, name, condition):

        if name in self._subsets.keys():
            raise ValueError(tc.RED+'[ERROR]: '+tc.RESET+f'Subset {name} already exists')
        self._subsets[name] = condition

    def __getitem__(self, key):
        return self._subsets[key]()

#############################################################

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
        
        self._data = pd.DataFrame()
        self._open(data, **kwargs)
        self._subsets = SubsetDict()

    def __getitem__(self, key):
        if ':' in key:
            key1, key2 = key.split(':')
            subset = self._subsets[key1]
            return subset[key2] if key2 else subset
        return self._data[key]
    
    def __setitem__(self, key, value):
        self._data[key] = value

    def _open(self, data, **kwargs):

        if isinstance(data, pd.DataFrame):
            self._data = data
        elif isinstance(data, str) or (isinstance(data, list) and all(isinstance(file, str) for file in data)):
            self._files = data if isinstance(data, list) else [data]
            for file in self._files:
                if file.endswith('.csv'):
                    print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+file+tc.RESET)
                    self._data = pd.concat([self._data, pd.read_csv(file, **kwargs)], ignore_index=True, copy=False)
                elif file.endswith('.parquet'):
                    print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+file+tc.RESET)
                    self._data = pd.concat([self._data, pd.read_parquet(file, **kwargs)], ignore_index=True, copy=False)
                else:
                    raise ValueError(tc.RED+'[ERROR]: '+tc.RESET+'Input data must be a list of .root or .csv files.')
        else:
            raise ValueError(tc.RED+'[ERROR]: '+tc.RESET+'Input data must be a string, a list of strings, or a pandas DataFrame.')
        
    @classmethod
    def from_root(cls, files, tree_name: str, folder_name: str = None, columns: list = None, **kwargs) -> 'Dataset':
        """
        Improved from_root: collects DataFrames in a list, uses uproot.concatenate when possible,
        handles missing trees, and concatenates once at the end.
        """

        files_list = []
        if isinstance(files, str) or (isinstance(files, list) and all(isinstance(file, str) for file in files)):
            files_list = files if isinstance(files, list) else [files]
        else: 
            raise ValueError(tc.RED+'[ERROR]: '+tc.RESET+'Input data must be a string or a list of strings.')
        
        uproot_kwargs = {k: v for k, v in kwargs.items() if k not in ['tree_name', 'columns', 'folder_name']}

        # If no folder_name wildcard and multiple files, try uproot.concatenate in one shot
        if folder_name is None:
            data = uproot.concatenate(
                {file: tree_name for file in files_list},
                filter_name=columns,
                library="pd",
                **uproot_kwargs
            )
            return cls(data)

        else:
            dfs = []
            wildcard = folder_name.endswith('*')
            base = folder_name[:-1] if wildcard else folder_name
            for file in files_list:
                f = uproot.open(file)
                
                if wildcard:
                    keys = list(f.keys())
                    _file_folders = [folder for folder in keys if (folder.startswith(base) and '/' not in folder)]
                    file_folders_duplicated = [folder.split(';')[0] for folder in _file_folders] # list with potentially duplicated folders
                    seen = {}
                    for idx, val in enumerate(file_folders_duplicated):
                        if val not in seen:
                            seen[val] = idx
                    file_folders = [_file_folders[idx] for idx in seen.values()]

                    if not file_folders:
                        print(tc.RED+'[WARNING]: '+tc.RESET+f'No folders matching "{base}*" in {file}, skipping.')
                    for folder in file_folders:
                        path = f'{file}:{folder}/{tree_name}'
                        print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+path+tc.RESET)
                        df_temp = uproot.open(path).arrays(filter_name=columns, library='pd', **uproot_kwargs)
                        dfs.append(df_temp)
                        
                else:
                    path = f'{file}:{folder_name}/{tree_name}'
                    print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+path+tc.RESET)
                    df_temp = uproot.open(path).arrays(filter_name=columns, library='pd', **uproot_kwargs)
                    dfs.append(df_temp)

            try:
                data = pd.concat(dfs, ignore_index=True, copy=False)
            except Exception as e:
                print(tc.RED+'[ERROR]: '+tc.RESET+f'Concatenation failed: {e}, falling back to default concat.')
                data = pd.concat(dfs, ignore_index=True)
            return cls(data)      
    
    @property
    def columns(self):
        return self._data.columns
    
    @property
    def shape(self):
        return self._data.shape
    
    @property
    def loc(self):
        return self._data.loc
    
    @property
    def data(self):
        return self._data
    
    @property
    def subsets(self):
        return self._subsets
    
    def add_subset(self, name, condition):
        '''
            Add a subset to the dataset by applying a condition.
        '''
        self._subsets.add_subset(name, lambda: self._data.loc[condition])

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
    
    def drop(self, labels=None, *, axis=0, index=None, columns=None, level=None, inplace=True, errors='raise') -> None | pd.DataFrame:
        '''
            Drop specified labels from the dataset.

            Args:
                labels: The labels to drop.
                **kwargs: Additional keyword arguments to be passed to the pandas drop function.
        '''
        
        if inplace:
            self._data.drop(labels=labels, axis=axis, index=index, columns=columns, level=level, inplace=True, errors=errors)
        else:
            return self._data.drop(labels=labels, axis=axis, index=index, columns=columns, level=level, inplace=False, errors=errors)

    def head(self, n: int = 5) -> pd.DataFrame:
        '''
            Return the first n rows of the dataset.

            Args:
                n (int): The number of rows to return.
        '''
        
        return self._data.head(n)

    def build_th1(self, column: str, axis_spec_x: AxisSpec, **kwargs) -> TH1F:
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
        
    def build_th2(self, column_x: str, column_y: str, axis_spec_x: AxisSpec, axis_spec_y: AxisSpec, **kwargs) -> TH2F:
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

    def build_boost1d(self, column: str, axis_spec_x: AxisSpec, **kwargs) -> bh.Histogram:
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
            return build_boost1(self._subsets[subset][column], axis_spec_x)
        else:
            return build_boost1(self._data[column], axis_spec_x)
        
    def build_boost2d(self, column_x: str, column_y: str, axis_spec_x: AxisSpec, axis_spec_y: AxisSpec, **kwargs) -> bh.Histogram:
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
            return build_boost2(self._subsets[subset][column_x], self._subsets[subset][column_y], axis_spec_x, axis_spec_y)
        else:
            return build_boost2(self._data[column_x], self._data[column_y], axis_spec_x, axis_spec_y)
