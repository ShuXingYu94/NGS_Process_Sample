import pandas as pd
import numpy as np
from functools import reduce
from tqdm import tqdm
import sys

def process_bar(num, total):
    rate = float(num)/total
    ratenum = int(100*rate)
    r = '\r[{}{}]{}%'.format('*'*ratenum,' '*(100-ratenum), ratenum)
    sys.stdout.write(r)
    sys.stdout.flush()

def uniq(ls):
    out = list(dict.fromkeys(ls))
    return out
def array_unique(lst):
    return dict(zip(*np.unique(lst, return_counts=True)))
def read_table(name,prefix='test_',subfix='.txt',folder='./stacks/',usecol=range(0,11)):
    address='{0}{1}{2}{3}'.format(folder,prefix,name,subfix)
    df=pd.read_table(address,header=None,on_bad_lines='warn',usecols=usecol,engine='python')
    if df.index[0]!=0:
        df=df.reset_index()
        df=df.iloc[:,:-1]
    return df

def process_table(df):
    # Get start-pos and end-pos
    tmp=df.copy()
    # tmp.columns=['Qname','FLAG','RNAME','POS','MAPQ','CIGAR','MRNM','MPOS','ISIZE','SEQ','QUAL',11,12,13,14,15,16,17]
    tmp.columns=['Qname','FLAG','RNAME','POS','MAPQ','CIGAR','MRNM','MPOS','ISIZE','SEQ','QUAL']
    len_ls=[]
    for record in list(tmp['SEQ']):
        len_ls.append(len(record))
    tmp.loc[:,'SEQ_LEN']=len_ls

    end_pos=[]
    start_pos=list(tmp['POS'])
    for ind in tmp.index:
        start=tmp.loc[ind,'POS']

        end=start+abs(tmp.loc[ind,'ISIZE'])
        end_pos.append(end)
    tmp.loc[:,'START_POS']=start_pos
    tmp.loc[:,'END_POS']=end_pos
    return tmp

def filter_df(tmp,ISIZE_upper=500,ISIZE_lower=-500,MAPQ=30,MRNM='='):
    processed_data=tmp.copy()
    processed_data=processed_data.loc[:,['RNAME','POS','MAPQ','MRNM','MPOS','ISIZE','SEQ_LEN','SEQ_LEN', 'START_POS', 'END_POS']]

    print('Before filter: {} records.'.format(processed_data.shape[0]))
    processed_data=processed_data[processed_data['MRNM']==MRNM]
    print('After filter with MRNM: {} records.'.format(processed_data.shape[0]))
    processed_data=processed_data[processed_data['ISIZE']<=ISIZE_upper]
    processed_data=processed_data[processed_data['ISIZE']>=ISIZE_lower]
    print('After filter with ISIZE: {} records.'.format(processed_data.shape[0]))
    processed_data=processed_data[processed_data['MAPQ']>=MAPQ]
    print('After filter with MAPQ: {} records.'.format(processed_data.shape[0]))
    return processed_data

def get_coverage_array(filtered_data,chromosome,sample_name):
    data=filtered_data[filtered_data['RNAME']==chromosome].copy()
    start_pos_unique = uniq(list(data['START_POS']))
    count_reads=0
    ary_result=np.arange(0,0,1) # Refresh the array for new chromosome

    for a in start_pos_unique:
        if data[data['START_POS']==a]['END_POS'].min() <a:
            print('End Position Error on {},{}bp.'.format(chromosome,a))
        count_reads+=1
        ary=np.arange(a,data[data['START_POS']==a]['END_POS'].max())
        ary_result=np.concatenate((ary_result,ary))
    ary_result=np.unique(ary_result)

    try:
        mean_len=len(ary_result)/count_reads
    except ZeroDivisionError:
        mean_len=0

    dic={'SAMPLE': sample_name, 'CHROMOSOME': chromosome,'UNIQ_START_COUNT':count_reads,'LENGTH_ARRAY':len(ary_result),'MEAN_LENGTH_PER_START':mean_len, 'ARRAY': ary_result}
    return dic


def operate(sample_names, ISIZE_upper=500, ISIZE_lower=-500, MAPQ=30, MRNM='=', overlap_rate=0.8):
    print('Operating with {} samples with {} overplapping rate:\n{}'.format(len(sample_names), overlap_rate,
                                                                            sample_names))
    condition = len(sample_names) * overlap_rate
    # Create dataframe
    ary_df = pd.DataFrame(
        columns=['SAMPLE', 'CHROMOSOME', 'UNIQ_START_COUNT', 'LENGTH_ARRAY', 'MEAN_LENGTH_PER_START', 'ARRAY'])
    # Get chromosome list
    print('Fteching Chromosome list...')
    #     for sample_name in sample_names:
    #         df=read_table(sample_name,prefix='',subfix='.txt',folder='./') # Read file
    #         for element in uniq(list(df.iloc[:,-1])):
    #             if element not in chr_ls:
    #                 chr_ls.append(element)
    chr_df = pd.read_table('./stacks/catalog.chrs.tsv')
    chr_ls = list(chr_df['# Chrom'])
    print('{} chromosomes are found.'.format(len(chr_ls)))

    # Filter reads with given condition.
    print('Start filtering reads with condition:\nInsert size:{1}~{0}\nMapping Quality:{2}\nIf paired:{3}'.format(
        ISIZE_upper, ISIZE_lower, MAPQ, MRNM))
    for sample_name in sample_names:
        print('Processing read data of {}'.format(sample_name))
        df = read_table(sample_name, prefix='', subfix='.txt', folder='./stacks/')  # Read file
        df = process_table(df)  # Acquire start-pos, end-pos, seq_length
        filtered_data = filter_df(df, ISIZE_upper, ISIZE_lower, MAPQ, MRNM)

        # Loop through
        n = 0
        for chromosome in chr_ls:
            n += 1
            process_bar(int(n / len(chr_ls) * 100), 100)
            dict_result = get_coverage_array(filtered_data, chromosome, sample_name)
            ary_df = ary_df.append([dict_result])
        print('Filtering Complete.')

    # Get consensus of certain overlapping rate
    print('Start to Calculate Overlapping-Rate...')
    consensus_len = {}
    consensus = {}
    for chromosome in chr_ls:
        tmp = np.concatenate(list(ary_df[ary_df['CHROMOSOME'] == chromosome]['ARRAY']))
        d1 = array_unique(tmp)
        d3 = {k: v for k, v in d1.items() if v >= condition}
        consensus[chromosome] = d3.keys()
        consensus_len[chromosome] = len(d3)
    ary_df['CONSENSUS_ARRAY'] = ary_df.apply(lambda x: consensus_len[x.CHROMOSOME], axis=1)
    ary_df['LENGTH_CONSENSUS'.format(overlap_rate)] = ary_df.apply(lambda x: consensus_len[x.CHROMOSOME], axis=1)
    ary_df['CONSENSUS_RATE'] = ary_df.apply(lambda x: x.LENGTH_CONSENSUS / x.LENGTH_ARRAY if x.LENGTH_ARRAY != 0 else 0,
                                            axis=1)
    print('Overlapping-Rate Calculation Complete.')
    print('Sorting DataFrame Based On Chromosome and Sample Name...')
    ary_df = ary_df.sort_values(by=['CHROMOSOME', 'SAMPLE'], ascending=True)  # Sort values (optional)
    print('Sorting Complete.')
    return ary_df


def calculate_mean(ary_df):
    print('Start to Calculate Overlapping-Rate Mean...')
    tmp1 = ary_df[['CHROMOSOME', 'LENGTH_ARRAY']].groupby('CHROMOSOME').mean().copy()
    tmp2 = ary_df[['CHROMOSOME', 'LENGTH_CONSENSUS', 'CONSENSUS_RATE']].groupby('CHROMOSOME').mean().copy()
    ary_mean = pd.concat([tmp1, tmp2], axis=1)
    tmp1 = ary_df[['CHROMOSOME', 'LENGTH_ARRAY']].groupby('CHROMOSOME').std().copy()
    tmp2 = ary_df[['CHROMOSOME', 'CONSENSUS_RATE']].groupby('CHROMOSOME').std().copy()
    ary_mean = pd.concat([ary_mean, tmp1, tmp2], axis=1)
    ary_mean.columns = ['LENGTH_ARRAY_mean', 'LENGTH_CONSENSUS_mean', 'CONSENSUS_RATE_mean', 'LENGTH_ARRAY_std',
                        'CONSENSUS_RATE_std']
    print('Overlapping-Rate Mean Calculation Complete.')
    return ary_mean

with open('./statistics/file_names.txt') as f:
    text=f.read()
    f.close()
filelist=[]
for a in text.split('\n'):
    if a != '':
        filelist.append(a)

filelist=[a.replace('.txt','') for a in filelist]
# General filtering conditions
overlap_rate=0.8
ISIZE_upper=1000 # Upper range of insetion size
ISIZE_lower=-1000    # Lower range of insetion size
MAPQ=30 # Map quality
MRNM='='    # If pair end reads are mapped on the same chromosome

result=operate(filelist,ISIZE_upper,ISIZE_lower,MAPQ,MRNM,overlap_rate)
result.to_csv('./statistics/compare_output.csv')
result_stat=calculate_mean(result)
result_stat.to_csv('./statistics/compare_statistic.csv')