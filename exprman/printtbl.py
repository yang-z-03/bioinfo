
def printtbl(table, styles = {}, appends = {}, nfrom = 0, nto = 999):
    
    vars = []
    widths = []
    cols = []
    objs = []
    fmtstring = '{_row:-4d} '
    fmtheader = '  no '
    rows = 10000 # some number that are excessively large.

    for k in table.keys():
        vars += [table[k]]
        width = max(len(k) + 2, len(max(table[k], key = len)))
        widths += [width]
        cols += [k]
        if rows >= len(table[k]): rows = len(table[k])
        fmtstring += '{' + k + ':' + str(width + (2 if k not in appends.keys() else (appends[k] + 2))) + '}'
        fmtheader += '{' + k + ':' + str(width + 2) + '}'

    # construct row data
    for i in range(rows):
        dicti = {}
        for col in cols:
            dicti[col] = table[col][i]
            if col in styles.keys(): dicti[col] = styles[col](dicti[col])
        objs += [dicti]

    # print header
    header_row = {}
    for col in cols: header_row[col] = col
    print(fmtheader.format(**header_row), end = '\n')
    total_width = 0
    print('---- ', end = '')
    for w in widths:
        print('-' * (w + 1) + ' ', end = '')
        total_width += (w + 2)
    
    print('', end = '\n')

    # print rows
    for i in range(max(nfrom, 0), 1 + min(nto, len(objs) - 1)):
        print(fmtstring.format(_row = i + 1, **objs[i]), end = '\n')
    
    return total_width + 5
