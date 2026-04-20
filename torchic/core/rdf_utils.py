from ROOT import TChain, TFile

from torchic.utils.terminal_colors import TerminalColors as tc

def prepare_input_tchain(input_data, tree_name, mode='DF', tchain_name='tchain'):

    file_data_list = input_data if isinstance(input_data, list) else [input_data]
    chain_data = TChain(tchain_name)

    for file_name in file_data_list:
      fileData = TFile(file_name)

      if mode == 'DF':
        for key in fileData.GetListOfKeys():
          key_name = key.GetName()
          if 'DF_' in key_name :
              print(f'Adding {tc.CYAN+tc.UNDERLINE}{file_name}/{key_name}/{tree_name}{tc.RESET} to the chain')

              chain_data.Add(f'{file_name}/{key_name}/{tree_name}')
      elif mode == 'tree':
        print(f'Adding {tc.CYAN+tc.UNDERLINE}{file_name}/{tree_name}{tc.RESET} to the chain')
        chain_data.Add(f'{file_name}/{tree_name}')
    
    
    return chain_data