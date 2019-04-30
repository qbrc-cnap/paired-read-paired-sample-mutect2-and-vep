import argparse
import itertools
import pandas as pd
import sys

def read_annotations(annotation_filepath):
    '''
    Tries to parse the annotation file.  If it cannot, issue return some sensible
    message
    '''
    generic_problem_message = '''A problem occurred when trying to parse your annotation 
        file, which was inferred to be in %s format.  
        Please ensure it follows our expected formatting.'''
    file_extension = annotation_filepath.split('.')[-1].lower()
    if file_extension == 'tsv':
        try:
            df = pd.read_csv(annotation_filepath, header=None, sep='\t')
        except Exception as ex:
            return (None, [generic_problem_message % 'tab-delimited'])
    elif file_extension == 'csv':
        try:
            df = pd.read_csv(annotation_filepath, header=None, sep=',')
        except Exception as ex:
            return (None, [generic_problem_message % 'comma-separated'])
    elif ((file_extension == 'xlsx') or (file_extension == 'xls')):    
        try:
            df = pd.read_excel(annotation_filepath, header=None)
        except Exception as ex:
            return (None, [generic_problem_message % 'MS Excel'])


    # now that we have successfully parsed something.  
    if df.shape[1] != 2:
        return (None, 
                ['The file extension of the annotation file was %s, but the' 
                ' file reader parsed %d column(s).  Please check your annotation file.  Could it have the wrong file extension?' % (file_extension, df.shape[1])])

    # check for NAs:
    if df.dropna().shape != df.shape:
        return (None, ['There were missing inputs in the annotation table.  Look for blank cells in particular.'])

    # we have two cols and had no NAs.  Now name the cols:
    return (df, [])

def write_matched_fastq(fastq_lookup, annot_df):
    '''
    Writes full matched fastq to stdout as tab delimited.
    '''
    for idx, row in annot_df.iterrows():
        fastqs = list(itertools.chain.from_iterable([fastq_lookup[col]
                                                     for col in row]))
        sys.stdout.write('\t'.join(fastqs) + '\n')


def get_commandline_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r1', required=True, nargs='+', dest='r1')
    parser.add_argument('-r2', required=True, nargs='+', dest='r2')
    parser.add_argument('-annot', required=True, dest='annot')
    args = parser.parse_args()
    return vars(args)


def main():
    '''
    Matches the ordered R1 and R2 files to the sample matching annotation file.
    '''
    arg_dict = get_commandline_args()
    fastq_lookup = dict(zip([f.split('/')[-1] for f in arg_dict['r1']],
                            zip(arg_dict['r1'], arg_dict['r2'])))
    annot_df, _ = read_annotations(arg_dict['annot'])
    write_matched_fastq(fastq_lookup, annot_df)


if __name__ == "__main__":
    main()
