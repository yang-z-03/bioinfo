
import requests
import json
import os
import time
import openai
from ansi import fore_green, fore_cyan, ansi_reset

# model = 'deepseek-r1-distill-qwen-14b'
# model = 'deepseek-r1-distill-llama-8b'
model = 'deepseek-ai/DeepSeek-R1'
apikey = 'sk-tydjnkebvgdqwipaiowmuqtpneahwqziyjaqbqzdwrmobngq'

# client = openai.OpenAI(base_url = "http://127.0.0.1:1234/v1", api_key = "lm-studio")
client = openai.OpenAI(base_url = 'https://api.siliconflow.cn/v1/', api_key = apikey)

url = "https://api.siliconflow.cn/v1/chat/completions"

payload = {
    "model": model,
    "messages": [],
    "stream": False,
    "max_tokens": 16384,
    "temperature": 0.1,
    "top_p": 0.7,
    "top_k": 50,
    "frequency_penalty": 2,
    "n": 1,
    "response_format": {"type": "text"}
}

headers = {
    "authorization": f"Bearer {apikey}",
    "content-type": "application/json"
}

session_system = [{
    "role": "system",
    "content": 
        'You should give your result in JSON format according to the given JSON schema. '
        'You should give your result entirely in English. '
        'You should not try to fix typos, just leave it as it is. '
        'Do not contain anything other than ["accession", "dtype", "strain", "tissue", "genotype", "age", "sex", "sort"] in your output JSON.'
}]

fprefix = []
facc = []
with open('configs/lookup', 'r') as flookup:
    lines = flookup.read()
    for line in lines.splitlines():
        if len(line) == 0: continue
        prefx, acc = line.split('\t')
        fprefix += [prefx]
        facc += [acc]

preface = [
    # file format description
    'Here, I will supply a database record from NCBI GEO. It is a dataset of high-throughput '
    'sequencing of mice. The record contains the summary information of the data series, and '
    'followed by a list of samples. (or just one sample). The summary of the series includes '
    'its accesssion code starting with GSE, the title of the study, an abstract, and experimental '
    'designs. The record of the sample includes the sample accession starting with GSM, sample '
    'name, sample characteristics which may contains the details of this animal and experimental '
    'groupings, a more detailed treatment protocol, and the protocols used for constructing the '
    'cDNA library for sequencing. You should place the inferred sample GSM accession into the'
    'output (field "accession")',

    # file type criteria
    'The dataset may contain data from several different high-throughput sequencing methods. '
    '(1) single cell RNA sequencing using 3\' test kit by 10x genomics. (notated by "3") '
    '(2) single cell RNA sequencing using 5\' test kit by 10x genomics. (notated by "5") '
    '(3) others. You should distinguish between these cases, and determine the data type for '
    'this sample. (field "dtype")',

    # sorting description
    'The cells may be sorted using one or more surface markers. The sorting process is often '
    'conducted using flow cytometry. In the process, people pick out the cells of interest '
    'by surface expression of certain proteins. These criteria should be stored in a list of '
    'surface marker names and their expression status. The expression status is denoted by '
    '(1) positive, or high expression, with "+", (2) negative, or low expression, with "-", '
    '(3) intermediate, with "int". You should not try to alter the marker name, just stick to '
    'what is given exactly in the text. You need not to check whether this protein name is '
    'valid, because these naming criteria are constantly changing. '
    # 'In most cases, the sorting criteria is stated using "[{marker}{+,-,int}]. '
    # 'and may be written without spacing. The status goes after the marker name, for example, '
    # 'CD3+CD4-CD8- indicates positive for CD3, while negative for CD4 and CD8. '
    'The JSON schema requires you to give the pair of marker names and its status in a JSON '
    'list. Surface marker names should be lowercased word. If no sorting is applied, leave this '
    'field blank. (field "sort")',

    'The animal taken in the experiment may be wild type, or genetically modified.'
    'Transgenic strains should be notated with tg(gene_symbol) where "gene_symbol" represent '
    'lowercased gene name, mutants should be notated with gene_symbol^(mutation_status) '
    'and "mutation_status" may describe coding alternations (for example "k18s") or'
    'allele heterozygosity (for example -/-, +/-). Do not try to alter the gene name, just '
    'turn them to lowercase of what stated exactly in the text. (field "genotype")',

    'The age (field "age") of the animal should be written in units of weeks. We assume a month '
    'to be 4 weeks. Embryonic age should be written exactly as the original text. '
    'If not specified, place "0" as age.',

    'The tissue source of the sample should be abbreviated using my nomenclature. We will abbreviate '
    'tissue source names as follows: (1) "bm" for bone marrow (2) "thy" for thymus (3) "ln" for '
    'lymph node (4) "lv" for liver (5) "spl" for spleen (6) "spi" for spinal cord (7) "int" for '
    'small intestine (8) "col" for colon (9) "bat" for brown adipose (10) "wat" for white adipose '
    '(11) "tumor:*" for specific tumor cell lines, the * stands for the name of cell line. You should'
    'place the abbreviated tissue string in (field "tissue").',

    'The gender of mice is encoded with "f" for female, "m" for male, and "?" for unknown. (field "sex")',

    'The strain of the mice should be formatted to fully-lowercased letters. For example, all '
    'strains derived from C57BL/6 mice should be notated as "c57bl/6". The strain name should '
    'be the base name of the strain without genotype information. These should be recorded in '
    'field "genotype" but not here. If the strain is not given, place "na" here. (field "strain")',

    # 'Last, you should re-organize a shortened identifier of the experiment groupings from the '
    # 'characteristics and processing protocol. This identifier (field "expm") '
    # 'indicates what differs in this sample. You should use at most 4 words to construct this '
    # 'identifier, and finally transform it to fully lowercased.',

    'The data record:'
]

epilog_gse = [
    'The file {filename} is one of the data files from the given dataset. You should infer from its '
    'name which sample does it belongs to, and extract the corresponding information from that sample. '
    'You should extract these information from the record in the format specified by the JSON schema. '
]

epilog_gsm = [
    'You should extract these information for the file "{filename}" from the only sample in the format '
    'specified by the JSON schema. Additional grouping information may be inferred from the file name. '
    'These abbreviations in the file name is not complex. You should not think too much.'
]

json_schema = {
    'type': 'json_schema',
    'json_schema': {
        "name": "sample.schema",
        "strict": True,
        "schema": {
            "type": "object",
            "properties": {
                "accession": { 
                    "type": "string",
                    "description": "Accession codes starting with GSM."
                },
                "dtype": {
                    "enum": ["3", "5", 'other'],
                    "description": "Data type for this high-throughput sequencing dataset."
                },
                "strain": { 
                    "type": "string",
                    "description": "Lowercased mouse strain name."
                },
                "tissue": { 
                    "type": "string",
                    "description": "Abbreviated tissue source."
                },
                "genotype": { 
                    "type": "string",
                    "description": "Genotype for the mouse."
                },
                "age": { 
                    "type": "integer",
                    "description": "Age in weeks."
                },
                "sex": {
                    "enum": ["f", "m", "?"],
                    "description": "Gender code."
                },
                "sort": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "marker": { "type": "string" },
                            "status": { "enum": ["+", "-", "int" ] }
                        },
                        "required": ["marker", "status"]
                    }
                }
            },
            "required": [
                "accession", "dtype", "strain", "tissue", "genotype",
                "age", "sex", "sort"
            ]
        }
    }
}

# here, gse dataset level queries.
id = 0
for prefx in fprefix:
    id += 1
    acctype, acc = facc[fprefix.index(prefx)].split(':')

    if os.path.exists(f'query/{acctype}/{acc}.json'):
        continue
    print(f'querying for sample [{prefx}] - {acctype}:{acc} ({id}/{len(fprefix)})')
    summary_fname = f'summary/{acctype}/{acc}.sum'
    with open(summary_fname, 'r', encoding = 'utf-8') as f:
        db = f.read()

    json_export = {
        'reasoning': '',
        'content': '',
        'usage': {
            'prompt': 0,
            'completion': 0,
            'total': 0
        }
    }

    json_dataset = {
        'accession': None,
        'dtype': None,
        'strain': None,
        'tissue': None,
        'genotype': None,
        'age': None,
        'expm': None,
        'sort': None,
        'sex': None
    }

    # each dataset is answered in one query session.
    session = []
    session += session_system
    text = '\n\n'.join(preface) + '\n\n' + db + '\n' + (
        epilog_gse[0].format(filename = prefx) if acctype == 'gse' else
        epilog_gsm[0].format(filename = prefx))
    
    session += [{ 'role': 'user', 'content': text }]
    payload['messages'] = session
    payload['response_format'] = json_schema

    raw_content = ''

    try:
        response = client.chat.completions.create(
            messages = session,
            model = model, 
            frequency_penalty = 2,
            max_tokens = 16384,
            response_format = json_schema,
            temperature = 0.1,
            stream = True,
            timeout = 900 # 15 minutes
        )

        is_message = False
        fore_green()

        for chunk in response:
            chunk_message = chunk.choices[0].delta.content
            chunk_reasoning = chunk.choices[0].delta.reasoning_content
            
            if chunk_message != None:
                if not is_message:
                    is_message = True
                    fore_cyan()
                print(chunk_message, end = '', flush = True)
                raw_content += chunk_message
            
            if chunk_reasoning != None:
                print(chunk_reasoning, end = '', flush = True)

        ansi_reset()
        
        pass
    
    except openai.OpenAIError as e:
        print(f'  (request error) - {acctype}:{acc}')
        print(e)
        continue

    import json
    try:
        # this content may not be proper json. we just write them down.
        with open(f'query/{acctype}/{acc}.json', 'w') as fo:
            fo.write(raw_content)
            print(f'  (saved) - {acctype}:{acc}')

    except:
        print(f'  (json parse fail) - {acctype}:{acc}')

    time.sleep(1) # sleep for 10s.

