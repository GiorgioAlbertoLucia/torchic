'''
    Abstract class for an analysis flow. 
    Implements the base functionalities for an analysis.
'''

from abc import ABC, abstractmethod
from torchic.core.dataset import Dataset
from torchic.core.histogram import AxisSpec
from torchic.utils.terminal_colors import TerminalColors as tc
import numpy as np
import yaml
import uproot

class AnalysisFlow(ABC):

    def __init__(self, dataset: Dataset, outfile_path: str = None, visual_config_file: str = None) -> None:
        '''
            - dataset (Dataset):  data to be preprocessed
            - dataset['full'] (Dataset): all entries in the dataset

        '''
        
        self._dataset = dataset
        self._available_subsets = ['full']

        self._outfile = None
        if outfile_path is not None:
            self._open_outfile(outfile_path)

        self._visual_config = None
        if visual_config_file is not None:
            with open(visual_config_file, 'r') as file:     
                self._visual_config = yaml.safe_load(file)
            
        
    def _open_outfile(self, outfile_path: str) -> None:
        '''
            Open the output file for writing.
            Parameters:
            - outfile_path (str): path to the output file
        '''
        print(f'{tc.GREEN}[INFO]{tc.RESET}: Creating output file {tc.UNDERLINE}{tc.CYAN}{outfile_path}{tc.RESET}')
        self._outfile = uproot.recreate(outfile_path)

    @property
    def dataset(self) -> Dataset:
        return self._dataset
    
    @property
    def available_subsets(self) -> list:
        return self._available_subsets

    def add_subset(self, name: str, condition: np.ndarray) -> None:
        '''
            Add a new subset to the dataset.
        '''
        self._dataset.add_subset(name, condition)
        self._available_subsets.append(name)

    def apply_cut(self, col_exp:str) -> None:
        '''
            Apply a cut on the dataset
            Parameters:
            - col_exp (str): column: expression to be evaluated (e.g. 'fPtHe3:fPtHe3 > 1.6')
        '''
        column, expression = col_exp.split(':')
        if column not in self._dataset.columns:
            print(tc.MAGENTA+'[WARNING]:'+tc.RESET+' Column',column,'not present in dataset!')
            return
        self._dataset.query(expression, inplace=True)

    def _visualize_plot(self, config: dict) -> None:
        '''
            Visualize a plot based on the configuration provided.
            Parameters:
            - config (dict): configuration dictionary for the plot
        '''
        if self._outfile is None:
            raise ValueError(f'{tc.RED}[ERROR]{tc.RESET}: Output file not opened! Please provide a valid output file path.')
        
        subset = config.get('opt', 'full')
        if subset not in self._available_subsets:
            print(f'{tc.MAGENTA}[WARNING]{tc.RESET}: Subset {subset} not available!')
            return
                
        if 'TH1' in config['type']:
            self._visualize_h1(config, subset)

        if 'TH2' in config['type']:
            self._visualize_h2(config, subset)

    def _visualize_h1(self, config: dict, subset: str) -> None:
        
        if config['xVariable'] not in self._dataset.columns:
            print(f'{tc.MAGENTA}[WARNING]{tc.RESET}: {config["xVariable"]} not present in dataset!')
            return 
        axis_spec_x = AxisSpec(config['nXBins'], config['xMin'], config['xMax'], config['name'], config['title'])
        hist = self._dataset.build_th1(config['xVariable']+config, axis_spec_x, subset=subset)
        
        hist_name = config['name']
        dirname = config.get('dir', 'None')
        if dirname != 'None':    
            self._outfile[f'{dirname}/{hist_name}'] = hist
        else:
            self._outfile[hist_name] = hist

    def _visualize_h2(self, config: dict, subset: str) -> None:
        
        if config['xVariable'] not in self._dataset.columns:
            print(f'{tc.MAGENTA}[WARNING]{tc.RESET}: {config["xVariable"]} not present in dataset!')
            return
        elif config['yVariable'] not in self._dataset.columns:
            print(f'{tc.MAGENTA}[WARNING]{tc.RESET}: {config["yVariable"]} not present in dataset!')
            return
        axis_spec_x = AxisSpec(config['nXBins'], config['xMin'], config['xMax'], config['name'], config['title'])
        axis_spec_y = AxisSpec(config['nYBins'], config['yMin'], config['yMax'], config['name'], config['title'])
        hist = self._dataset.build_th2(config['xVariable'], config['yVariable'], axis_spec_x, axis_spec_y, subset=subset)

        hist_name = config['name']
        dirname = config.get('dir', 'None')
        if dirname != 'None':
            self._outfile[f'{dirname}/{hist_name}'] = hist
        else:
            self._outfile[hist_name] = hist

    def visualize(self, plots: list) -> None:
        ''' 
            Visualization of data.

            - plots (list): list of plots names to produce. The plots are defined in a configuration file.
        '''

        if self._outfile is None:
            raise ValueError(f'{tc.RED}[ERROR]{tc.RESET}: Output file not opened! Please provide a valid output file path.')
        if self._visual_config is None:
            raise ValueError(f'{tc.RED}[ERROR]{tc.RESET}: Visual configuration not loaded! Please provide a valid visual configuration file.')

        for plot in plots:
            if plot not in self._visual_config:
                print(tc.MAGENTA+'[WARNING]:'+tc.RESET+' Plot',plot,'not found in visual configuration!')
                continue
            
            config = self._visual_config[plot]
            self._visualize_plot(config)
