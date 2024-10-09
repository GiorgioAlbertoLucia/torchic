import pandas as pd
import uproot

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
        if tree_name is None:  
            print(tc.RED+'[ERROR]: '+tc.RESET+'tree_name must be specified when using a .root file.')
            return

        columns = kwargs.get('columns', None)
        folder_name = kwargs.get('folder_name', None)
        if folder_name is None:
            return uproot.open(f'{file}:{tree_name}').arrays(filter_name=columns, library='pd', **kwargs)

        if folder_name[-1] != '*':
            print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+f'{file}:{folder_name}/{tree_name}'+tc.RESET)
            return uproot.open(f'{file}:{folder_name}/{tree_name}').arrays(filter_name=columns, library='pd', **kwargs)

        file_folders = uproot.open(file).keys()
        tree_path_list = []
        for folder in file_folders:
            if folder_name[:-1] in folder and tree_name in folder:
                tree_path_list.append(folder)

        tmp_data = pd.DataFrame()
        for tree_path in tree_path_list:
            print(tc.GREEN+'[INFO]: '+tc.RESET+'Opening file: '+tc.UNDERLINE+tc.BLUE+f'{file}:{tree_path}'+tc.RESET)
            print(f"Reading {file}:{tree_path}")
            tmp_data = pd.concat([tmp_data, uproot.open(f'{file}:{tree_path}').arrays(filter_name=columns, library='pd', **kwargs)], ignore_index=True, copy=False)
        return tmp_data
    
    @property
    def data(self):
        return self._data
    
    @property
    def subsets(self):
        return self._subsets
    
    def add_subset(self, name, condition):
        self._subsets.add_subset(name, lambda: self._data[condition])