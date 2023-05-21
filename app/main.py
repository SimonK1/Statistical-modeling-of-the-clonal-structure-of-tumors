
from sklearn.base import TransformerMixin
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import io
import os

"""
Custom loading function for VCF file loading. 
This functions skips entire head section and categorizes data into columns.
Function creates readable dataframe.
The function is based on the latest VCF standard.
Created with an emphasis on low time complexity.
"""
def load_vcf_file(path):
    with open(path, 'r') as f:
        # Skipping lines starting with ## - Head
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(io.StringIO(u''.join(lines)), dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str}, sep='\t').rename(columns={'#CHROM': 'CHROM'})

class VcfDataExtractionTransformer(TransformerMixin):
    """
    A tranformer for extracting needed information from VCF file.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialize method. No arguments were needed
        """
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    def transform(self, X, **transform_params):
        # Splitting sample column according to standardized FORMAT column
        X[['GENOTYPE', 'ALLELIC DEPTH', 'DEPTH', 'GENOTYPE QUALITY', 'PHRED-SCALED LIKELIHOODS']] = X.iloc[:, 9].str.split(":", expand=True)

        # Retyping into integers to enable numeric operations
        X["GENOTYPE QUALITY"] = X["GENOTYPE QUALITY"].astype('int')
        X['POS'] = X['POS'].astype('int')

        # Finding allele frequencies and their respective copy number. Rewriting copy number into usable format
        temp = X['INFO'].str.split(";", expand=True)
        X['ALLELIC FREQUENCY'] = temp[1].str[3:]
        X['MINOR ALLELE COPY NUMBER'] = np.where(X['GENOTYPE'] == '0/1', 0, 1)
        X['MAJOR ALLELE COPY NUMBER'] = np.where(X['GENOTYPE'] == '0/1', 2, 1)
        return X
    
class FilterQualityTransformer(TransformerMixin):
    """
    A tranformer for filtering only quality reads.
    """
    
    def __init__(self, *args, percentage, **kwargs):
        """
        Initialize method.

        :param percentage: The limit that the read quality must meet. Expressed as percentages
        """
        self.percentage = percentage
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    def transform(self, X, **transform_params):
        """
        Function that uses GENOTYPE QUALITY to determine the read quality. 
        This function was more successful than the division of data into quantiles. 
        With a large amount of high-quality data, high-quality reads were lost because they did not fit into the quantile.
        """
        X = X[X["GENOTYPE QUALITY"] >= self.percentage]
        return X

class CopyCallsMergeTransformer(TransformerMixin):
    """
    A tranformer for merging copy number results with VCF file.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialize method.

        :param *args: dataset created from copy caller results
        """
        self.data = args[0]
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    """
    Custom function for finding suitable copy number segments.
    This function takes reads from the VCF file and searches CNV segments 
    to check if it finds a suitable segment which it can assign this read to.
    
    """
    def find_range(self, pos):
        for i, row in self.data.iterrows():
            if row[4] <= pos <= row[5]:
                return row[7]
        return None

    
    def transform(self, X, **transform_params):
        # Retyping formats into integers to enable numerical operations
        self.data['CNV level'] = self.data['CNV level'].apply(np.round) 
        self.data['CNV level'] = self.data['CNV level'].astype('int')
        self.data['CNV Region Start'] = self.data['CNV Region Start'].astype('int')
        self.data['CNV Region End'] = self.data['CNV Region End'].astype('int')

        # Filtering out segments with lost data
        self.data = self.data[self.data['CNV level'] != 0]

        # Finding suitable segments for reads
        X['COPY NUMBER'] = X['POS'].apply(self.find_range)

        # Filling segments without CNVs with default value
        X['COPY NUMBER'] = X['COPY NUMBER'].fillna(2.0)

        # Fitlering out unwanted chromosomes
        X = X.loc[X['CHROM'].isin(['chr1', 'chr2', 'chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10', 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chr22', 'chrX' , 'chrY', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10','11', '12', '13', '14', '15','16','17','18', '19', '20','21','22','X','Y','x','y', 'ch1', 'ch2', 'ch3','ch4','ch5','ch6','ch7','ch8','ch9', 'ch10', 'ch11','ch12','ch13','ch14','ch15','ch16','ch17','ch18','ch19','ch20','ch21','ch22','ch22', 'chX' , 'chY'])]
        return X
    
class PyCloneTransformer(TransformerMixin):
    """
    A tranformer for preparing input files for tools PyClone and PyClone-VI
    """
    
    def __init__(self, samples, *args, **kwargs):
        """
        Initialize method.

        :param samples: the number of samples that the transformer randomly selects for the PyClone input file
        """
        self.samples = samples
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    def transform(self, X, **transform_params):
        """
        Function to prepare dataframe suitable for PyClone and PyClone-VI.
        Subsequently this dataframes are exported into .tsv files.
        """

        pyclone_df = pd.DataFrame()
        # Create custom mutation_id according to PyClone docu
        pyclone_df['mutation_id'] = X['CHROM'] + ':' + X['POS'].astype(str)

        # Extract allelic counts and divide them into singular columns
        pyclone_df[['var_counts', 'ref_counts', 'none']] = X['ALLELIC DEPTH'].str.split(',',expand=True)

        # Retype into integer to allow numeric operations
        pyclone_df['var_counts'] = pyclone_df['var_counts'].astype(int)
        pyclone_df['ref_counts'] = pyclone_df['ref_counts'].astype(int)

        # Copy needed values
        pyclone_df['minor_cn'] = X['MINOR ALLELE COPY NUMBER'].values
        pyclone_df['major_cn'] = X['MAJOR ALLELE COPY NUMBER'].values
        pyclone_df['normal_cn'] = X['COPY NUMBER'].values
        pyclone_df['normal_cn'] = pyclone_df['normal_cn'].astype(int)

        # Count Variant frequency
        pyclone_df['variant_freq'] = pyclone_df['var_counts'] / (pyclone_df['ref_counts'] + pyclone_df['var_counts']) 

        # Reformat Genotype 
        pyclone_df['genotype'] = X['GENOTYPE']
        pyclone_df['genotype'] = np.where(X['GENOTYPE'] == '0/1', 'AB', 'BB')
        pyclone_df = pyclone_df.drop('none', axis=1)

        # Take 300 random samples
        sampled_data = pyclone_df.sample(n=self.samples)

        # Export into .tsv
        sampled_data.to_csv('./inputFiles/PyClone/PyCloneInput.tsv', sep="\t", index=False)

        # Add needed columns for PyClone-VI and rename existing ones
        pyclone_df['sample_id'] = 'R1'
        pyclone_df['alt_counts'] = pyclone_df['var_counts']
        pyclone_df = pyclone_df.drop('var_counts', axis=1)

        # Export into .tsv
        pyclone_df.to_csv('./inputFiles/PyCloneVI/PyCloneVIInput.tsv', sep="\t", index=False)
        return X
    
class SciCloneTransformer(TransformerMixin):
    """
    A tranformer for preparing input file for tool SciClone.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialize method.

        :param *args: dataset created from copy caller results
        """
        self.data = args[0]
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    def transform(self, X, **transform_params):
        # Creating empty dataframes
        sciclone_df = pd.DataFrame()
        sciclone_df2 = pd.DataFrame()
        temp = pd.DataFrame()

        # Copy values into dataframe
        sciclone_df['chr'] = X['CHROM'].values
        sciclone_df['pos'] = X['POS'].values

        # Create temporary dataframe for storing splitted data to prevent wrong column order
        temp[['var_counts', 'ref_counts', 'none']] = X['ALLELIC DEPTH'].str.split(',',expand=True)
        temp['ref_counts'] = temp['ref_counts'].astype(int)
        temp['var_counts'] = temp['var_counts'].astype(int)

        # Copy splitted data from temporary dataframe
        sciclone_df['ref_reads'] = temp['ref_counts'].values
        sciclone_df['var_reads'] = temp['var_counts'].values

        # Count variant frequency
        temp['variant_freq'] = temp['var_counts'] / (temp['ref_counts'] + temp['var_counts']) 
        sciclone_df['vaf'] = temp['variant_freq'].values * 100

        # Retype columns to integer
        sciclone_df['ref_reads'] = sciclone_df['ref_reads'].astype(int)
        sciclone_df['var_reads'] = sciclone_df['var_reads'].astype(int)

        # Exporting into .tsv file
        sciclone_df.to_csv('./inputFiles/SciClone/ScicloneVafFile.tsv', sep="\t", index=False, header=False)

        # Creating new dataframe with correct order of needed columns and correct names
        sciclone_df2['chr'] = self.data['chr'].values
        sciclone_df2['start'] = self.data['CNV Region Start'].values
        sciclone_df2['stop'] = self.data['CNV Region End'].values
        sciclone_df2['segment_mean'] = self.data['CNV level'].values
        sciclone_df2['segment_mean'] = sciclone_df2['segment_mean'].fillna(2.0)
        sciclone_df2['segment_mean'] = sciclone_df2['segment_mean'].astype('int')

        # Exporting into .tsv file
        sciclone_df2.to_csv('./inputFiles/SciClone/ScicloneCopyFile.tsv', sep="\t", header=False, index=False)
        return X
    
class CoverageFileTransformer(TransformerMixin):
    """
    A tranformer for creating coverage file. This file is later used in tool TitanCNA.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialize method.

        :param *args: dataset created from copy caller results
        """
        self.data = args[0]
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    def transform(self, X, **transform_params):
        # Create empty dataframes
        temp = pd.DataFrame()
        temp2 = pd.DataFrame()
        coverage_df = pd.DataFrame()

        # Copy needed values into temporary dataframe so we later switch this dataframe by one row
        temp['chr'] = self.data['chr'].values
        temp['start'] = self.data['CNV Region Start'].values
        temp['stop'] = self.data['CNV Region End'].values
        temp['segment_mean'] = self.data['CNV level'].values
        temp['segment_mean'] = temp['segment_mean'].fillna(2.0)
        temp['segment_mean'] = temp['segment_mean'].apply(np.round)
        temp['segment_mean'] = temp['segment_mean'].astype('int')

        # Filter segments without copy number information
        temp = temp[temp['segment_mean'] != 0]  
        """
        Shift segment ends by one row. This creates new segments where start represents stop of the segment and vice versa.
        By this quick approach we supply missing segments.
        """
        temp['prev_end'] = temp['stop'].shift()
        new_df = pd.DataFrame({'chr': temp['chr'],'end': temp['start'],'start': temp['prev_end'], })

        # Add original segments
        temp2['chr'] = temp['chr'].values
        temp2['start'] = temp['start'].values
        temp2['end'] = temp['stop'].values
        temp2['logR'] =  temp['segment_mean'].values

        # Merge segments into one dataframe and sort them
        merge_result = pd.concat([temp2, new_df])
        merge_result = merge_result.sort_values(by=['chr'])

        # Copy values into new dataframe in correct order
        coverage_df['chr'] = merge_result['chr'].values
        coverage_df['start'] = merge_result['start'].values
        coverage_df['end'] = merge_result['end'].values
        coverage_df['logR'] = merge_result['logR'].values
        coverage_df['logR'] = coverage_df['logR'].fillna(2.0)
        coverage_df['logR'] = coverage_df['logR'].astype('int')

        # Filter little segments and unwanted chromosomes
        coverage_df = coverage_df[coverage_df['end'] - coverage_df['start'] > 2]
        coverage_df = coverage_df.loc[coverage_df['chr'].isin(['chr1', 'chr2', 'chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10', 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chr22', 'chrX' , 'chrY', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10','11', '12', '13', '14', '15','16','17','18', '19', '20','21','22','X','Y','x','y', 'ch1', 'ch2', 'ch3','ch4','ch5','ch6','ch7','ch8','ch9', 'ch10', 'ch11','ch12','ch13','ch14','ch15','ch16','ch17','ch18','ch19','ch20','ch21','ch22','ch22', 'chX' , 'chY'])]
        
        # Export file into .cn format
        coverage_df.to_csv('./inputFiles/TitanCNA/Titancna.cn', sep="\t", index=False)
        return X

class TitanCNATransformer(TransformerMixin):
    """
    A tranformer for preapring input files for TitanCNA.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialize method. No parameters were needed
        """
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    def transform(self, X, **transform_params):
        # Create ampty dataframes
        titnacna_df = pd.DataFrame()
        temp = pd.DataFrame()
        
        # Extract allelic depth
        temp[['var_counts', 'ref_counts', 'none']] = X['ALLELIC DEPTH'].str.split(',',expand=True)

        # Retype values into integer to enable numeric operations
        temp['ref_counts'] = temp['ref_counts'].astype(int)
        temp['var_counts'] = temp['var_counts'].astype(int)

        # Copy all needed values
        titnacna_df['chr'] = X['CHROM'].values
        titnacna_df['posn'] = X['POS'].values
        titnacna_df['ref'] = X['REF'].values
        titnacna_df['refCount'] = temp['ref_counts'].values
        titnacna_df['Nref'] = X['ALT'].values
        titnacna_df['NrefCount'] = temp['var_counts'].values

        # Filter unwanted chromosomes
        titnacna_df = titnacna_df.loc[titnacna_df['chr'].isin(['chr1', 'chr2', 'chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10', 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chr22', 'chrX' , 'chrY', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10','11', '12', '13', '14', '15','16','17','18', '19', '20','21','22','X','Y','x','y', 'ch1', 'ch2', 'ch3','ch4','ch5','ch6','ch7','ch8','ch9', 'ch10', 'ch11','ch12','ch13','ch14','ch15','ch16','ch17','ch18','ch19','ch20','ch21','ch22','ch22', 'chX' , 'chY'])]
        
        # Export into .het file
        titnacna_df.to_csv('./inputFiles/TitanCNA/Titancna.het', sep="\t", index=False)
        return X

class FastCloneTransformer(TransformerMixin):
    """
    A tranformer to create input data for FastClone tool.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialize method.
        """
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    def transform(self, X, **transform_params):
        fastclone_df = pd.DataFrame()

        # Create custom mutation_id according to FastClone docu
        fastclone_df['mutation_id'] = X['CHROM'] + ':' + X['POS'].astype(str)
        fastclone_df[['ref_counts', 'var_counts', 'none']] = X['ALLELIC DEPTH'].str.split(',',expand=True)

        # Retype columns into integer to enable numeric operations
        fastclone_df['var_counts'] = fastclone_df['var_counts'].astype(int)
        fastclone_df['ref_counts'] = fastclone_df['ref_counts'].astype(int)
        
        # Copy needed values
        fastclone_df['minor_cn'] = X['MINOR ALLELE COPY NUMBER'].values
        fastclone_df['major_cn'] = X['MAJOR ALLELE COPY NUMBER'].values
        fastclone_df['normal_cn'] = X['COPY NUMBER'].values
        fastclone_df['normal_cn'] = fastclone_df['normal_cn'].astype(int)

        # Count variant frequency
        fastclone_df['variant_freq'] = fastclone_df['var_counts'] / (fastclone_df['ref_counts'] + fastclone_df['var_counts']) 

        # Copy and reformat genotype
        fastclone_df['genotype'] = X['GENOTYPE']
        fastclone_df['genotype'] = np.where(X['GENOTYPE'] == '0/1', 'AB', 'BB')
        fastclone_df = fastclone_df.drop('none', axis=1)

        # Export into .tsv file
        fastclone_df.to_csv('./inputFiles/FastClone/FastCloneInput.tsv', sep="\t", index=False)
        return X

    
class TitanCNA(TransformerMixin):
    """
    A transfomer which runs TitanCNA script.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialize method. No arguments were needed
        """
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    def transform(self, X, **transform_params):
        """
        This function runs TitanCNA R script and takes as input file prepared data files.
        Het file is created by TitanCNATransforme anc Cn file is created by CoverageFileTransformer.
        """
        os.system('Rscript ./scripts/titanCNA.R --id TitanCNA --hetFile ./inputFiles/TitanCNA/Titancna.het --cnFile ./inputFiles/TitanCNA/Titancna.cn --chrs "c(1:22, \\"X\\")" --estimatePloidy TRUE --outDir ./results/TitanCNA/')
        return X
    
class SciClone(TransformerMixin):
    """
    A transfomer which runs SciClone script.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialize method. No arguments were needed
        """
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    def transform(self, X, **transform_params):
        """
        This function runs SciClone R script and does not take any extra arguments. 
        R script loads all input files by itself.
        All input files were prepared by SciCloneTransfomer.
        """
        os.system('Rscript ./scripts/sciClone.R')
        return X
    
class FastClone(TransformerMixin):
    """
    A function that runs FastClone tool.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialize method. No arguments were needed.
        """
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    def transform(self, X, **transform_params):
        """
        This function runs FastClone and takes created data files as input.
        All input files were prepared by FastCloneTransformer.
        """

        os.system('fastclone load-pyclone prop ./inputFiles/FastClone/FastCloneInput.tsv None solve ./results/FastClone/')
        return X
    

 
class PyCloneVI(TransformerMixin):
    """
    A function that runs PyClone-VI tool.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialize method. No arguments were needed
        """
        
    def fit(self, X, y=None):
        """
        Fits transformer over data.
        """
        return self
    
    def transform(self, X, **transform_params):
        """
        This function runs PyClone-VI tool and takes as input prepared data files by PyCloneTransformer.
        Subsequently, this function prepares histplot which serves as a great display of cellular prevalence.
        """

        os.system('pyclone-vi fit -i ./inputFiles/PyCloneVI/PyCloneVIInput.tsv -o ./results/PyCloneVI/PyCloneVI.h5 -c 40 -d beta-binomial -r 10')
        
        # Reading results and creating histplot. Subsequently exporting histplot.
        data = pd.read_csv('./results/PyCloneVI/PyCloneVI.tsv', sep="\t")
        sns.histplot(data=data, x='cellular_prevalence', hue='cluster_id', multiple='stack', kde=False)
        plt.savefig('./results/PyCloneVI/PyCloneVI_Plot.png')
        return X
    