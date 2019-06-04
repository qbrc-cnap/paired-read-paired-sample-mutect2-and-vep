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
    else:
        return (None, ['Your annotation file did not have the expected extension.  We found an extension of "%s", but expected one of: csv, tsv, or Excel.' % file_extension])


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

def check_fastq_existence(fastq_lookup, annot_df):
    '''
    Writes full matched fastq to stdout as tab delimited.
    '''
    errs = []
    err_message = 'The match annotation file indicates a file not uploaded:'
    for idx, row in annot_df.iterrows():
        for col in row:
            if col not in fastq_lookup:
                errs.append(' '.join([err_message, col]))
    return errs


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
    fastq_lookup = set([f.split('/')[-1] for f in arg_dict['r1']])
    annot_df, err_list = read_annotations(arg_dict['annot'])
    err_list.extend(check_fastq_existence(fastq_lookup, annot_df))
    if len(err_list) > 0:
        sys.stderr.write('#####'.join(err_list)) # the 5-hash delimiter since some stderr messages can be multiline
        sys.exit(1) # need this to trigger Cromwell to fail


if __name__ == "__main__":
    main()
