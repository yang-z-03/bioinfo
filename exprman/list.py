
import os
from printtbl import printtbl
from ansi import error

def styles_acctype(x: str):
    if x == 'gse':
        return '\033[32m' + x + '\033[0m'
    return '\033[31m' + x + '\033[0m'

def styles_number(x: str):
    if x.isdecimal():
        return '\033[33;1m' + x + '\033[0m'
    return '\033[31;0m' + x + '\033[0m'

def display(registry, nrow, show_status = False):

    # print the basic table.

    if not show_status:
        printtbl(
            registry, styles = { 'type': styles_acctype, 'accession': styles_number },
            appends = { 'type': 9, 'accession': 11 }
        )
    
    else:

        from download import status_download, styles_download
        from download import status_data, styles_data
        from metagen import styles_dtype, status_dtype
        from qc import status_qc, styles_qc

        registry['download'] = status_download(registry)
        registry['data'] = status_data(registry)
        registry['dtype'] = status_dtype(registry)
        registry['qc'] = status_qc(registry)

        printtbl(
            registry, styles = { 
                'type': styles_acctype, 'accession': styles_number,
                'download': styles_download, 'data': styles_data,
                'dtype': styles_dtype, 'qc': styles_qc
            },
            appends = { 
                'type': 9, 'accession': 11,
                'download': 11, 'data': 11, 'dtype': 11, 'qc': 11
            }
        )

    pass