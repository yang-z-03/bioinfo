
import os
from ansi import error, warning, info

def display(registry, n_row, show_status = False):

    config = {}

    current_dataset = '*'
    current_sample = '*'
    
    if os.path.exists('configs/samples'):
        with open('configs/samples', 'r') as fr:
            contents = fr.read().splitlines()
            for line in contents:
                if len(line.strip()) == 0: continue
                if line.strip().startswith('#'): continue

                tokens = line.strip().split()
                if len(tokens) < 2: 
                    error('syntax error:', line)
                    continue

                if tokens[0] == 'dataset':
                    
                    if tokens[1] in config.keys():
                        warning(f'duplicated dataset [{tokens[1]}]')
                        
                    current_dataset = tokens[1]
                    config[current_dataset] = {}
                    current_sample = '*'
                    continue

                elif tokens[0] == 'sample':
                    
                    if tokens[1] in config[current_dataset].keys():
                        warning(f'duplicated sample [{current_dataset}/{tokens[1]}]')
                    
                    current_sample = tokens[1]
                    config[current_dataset][current_sample] = {}
                    continue

                elif tokens[0] == 'prop':
                    if len(tokens) <= 2: continue
                    if current_dataset not in config.keys():
                        config[current_dataset] = {}
                    if current_sample not in config[current_dataset].keys():
                        config[current_dataset][current_sample] = {}
                    config[current_dataset][current_sample][tokens[1]] = tokens[2]
    
    else: error('configs/samples do not exist.')
    pass
