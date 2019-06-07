import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt
import math
import re
import csv
import argparse
import datetime
#from matplotlib import pylab as p

#test

#THIS Program was adjusted to analyse any amount of STD's en RAT samples (number of species not important)
#args are added
#new ways of filtering included : the script checks the loci of the jenamono's :
#if the highest number of reads for a loci is allocated to the correct mono...
#if not is A: filters it out or B: if the allocation is very 'pure' it renames the locus name to
#corresponding jenamono name
#it also checks if there are not to many reads allocted to the wrong jenamono causing higher background in
#zero species samples...

# names of mono's : jenamono1, jenamono2 etc...
# names of loci's : jenamono1_1, jenamono1_2, etc...
# This script uses two mono pools devined as pool 1 and 2 combined with two lists for renaming standaard to STD_
# The standards used for pool 1 or two are renamed from standaard to STD. only the STD are taken for analysis.
# Ratio samples are named ratio...
# species dict is the dictionary that is used to rename the jenomono's to original names.

pd.set_option('display.max_columns', 1000)
pd.options.mode.chained_assignment = None # getting rid of som error messages
#plt.style.use('bmh')

def parse_args():
    """Pass command line arguments"""
    parser = argparse.ArgumentParser(description='morph SGBS stats into species percentages (after filtering)')
    #input files
    parser.add_argument('-i', '--input',
                        help='tab delimited csv input file from msGBS_STATS.py')
    parser.add_argument('-op', '--output_prefix',
                        help='general output location (folder) and first part of name that endswith _'
                             'example: /Users/NielsWagemaker/SGBS_NIELS/pilot_3/mapping_month_refxx_mapXX_Exx_')
    parser.add_argument('-os', '--output_suffix',
                        help='final part of final output endswith .csv',default='filter1_filter2_minx_poolx.csv')
    parser.add_argument('-f1', '--filter_1',
                        help='filtering parameter_1, minumum reads threshold',default=8)
    parser.add_argument('-f2', '--filter_2',
                        help='filtering parameter_2, sensitivity threshold',default=15)
    parser.add_argument('-f3', '--filter_3',
                        help='filtering parameter_2, minimum read count per sample',default=1000)
    parser.add_argument('-p', '--pool',
                        help='pool_to_analyse',default=0)
    parser.add_argument('-e', '--extra',
                        help='write tussen.sum file, ', default='1')

    args = parser.parse_args()
    return args

def calculus(args):
    print "input file (created by msGBS_STATS.py) : ", args.input
    print
    print "output_prefix : ", args.output_prefix
    print "output_suffix : ", args.output_suffix
    print
    print "filter_1 : ", args.filter_1
    print "filter_2 : ", args.filter_2
    print "pool : ", args.pool
    print
    print "time started : ", datetime.datetime.now()

    '''reading input file'''

    df = pd.read_csv(args.input, sep='\t| |;', engine='python')

    ''' Sample types : monosamples : should start with: jenamono<number>
                     : standards used for standardization should begin with STD
                     : ratio samples; unknown samples that are known begin with ratio
                     : experimental samples; have yet to implemenmt'''

    df_jenamono = df.filter(regex='^Locus|^jenamono', axis=1)

    if args.extra:
        #as it is now it produces a list of loci and reads per mono
        csv_name = (args.output_prefix + 'e1_jenamono_test_%s_%s.csv' % (args.filter_1,args.filter_2))
        df_jenamono.to_csv(csv_name, sep='\t', decimal=',')
        print "written : %s" % csv_name
    else:
        print "extra files not written :" % csv_name
        print "Change --extra.args to 'yes' if this or other extra files are needed"

    df_jenamono.set_index('Locus', inplace=True)

    df_sum_sample_jenamono = df_jenamono.sum(axis=0)

    df_jenamono_div = df_jenamono.div(df_sum_sample_jenamono.squeeze())


    df_jenamono_div['max']=df_jenamono_div.idxmax(axis=1)
    getridof=[]

    df.set_index('Locus', inplace=True)


    '''FILTERING or renaming of loci based on mapping of monoculture samples: if too many reads 
     of a mono sample map to other mono reference.
     This is called emperical filtering.
     Basic parameters:
     
     Filter 1 : 2  # effectively this means that a minimum of 2 reads are needed to accept the loci
     Filter 2 : 15 # effectively this meand that Total reads (all mono's)> (reads mono +(reads mono / 15))
                     This reduces background signal.
                     If this filter fails it checks if the mono that has the highest relative number of reads meets
                     the loci filters requirements, if so the loci is renamed'''

    number =1000000 #renamed loci start at laci number 1000000
    tempcolumns = ["original_locus_name", "new locus name"]
    df_renamed_loci = pd.DataFrame(columns=tempcolumns)

    for index, row in df_jenamono_div.iterrows():
        locus = index
        max = row[-1]
        maxi = row[-1]+'_'
        number += 1
        SUM =row.filter(regex='jenamono').sum(axis=0)
        high = row.get(max)

        if maxi not in locus:
            reads_other_mono_sum = SUM - high
            #Here we check if the mono with the highest number of reads corresponds to the locus name
            #if not the locus is deleted or in some cases renamed
            #all tro be deleted loci are collected in list 'getridof'
            # #the df dataframe is exported to csv before and after filtering
            # if high > (10*reads_other_mono_sum):
            #     #This variable (10) was not varied yet... > filter 0
            #     new_locus = (maxi+str(number))
            #     df_renamed_loci.loc[number] = ["%s"%locus,"%s"%new_locus]
            #     df = df.rename(index={locus:new_locus})
            #     df.ix[new_locus, 'Species']=max
            #     locus = new_locus
            #
            #     mono = locus.split('_')
            #     mono = mono[0]
            #     summie = df_sum_sample_jenamono[mono]
            #     test = (int(args.filter_1) / float(summie))
            #
            #     if high < test:
            #         getridof.append(locus)
            #     if SUM > (high + (high / int(args.filter_2))):
            #         getridof.append(locus)
            #Onderstaand commando was een else:
            getridof.append(locus)
        else:
            mono = locus.split('_')
            mono = mono[0]
            summie = df_sum_sample_jenamono[mono]
            test =(int(args.filter_1)/float(summie))
            if high < test:
                getridof.append(locus)
            if SUM > (high + (high / int(args.filter_2))):
                getridof.append(locus)
    if args.extra:
        csv_name = (args.output_prefix+"e2_renamed_loci_filters"+args.output_suffix)
        csv_name2 = (args.output_prefix + "e3_before_filter_%s_%s.csv"%(args.filter_1, args.filter_2))
        df_renamed_loci.to_csv(csv_name, sep='\t', decimal=',')
        df.to_csv(csv_name2, sep='\t', decimal=',')
        print "written : %s" % csv_name
        print "written : %s" % csv_name2

    #Here i want to save the total sum per ratio sample before filtering
    df_ratio = df.filter(regex='^ratio', axis=1)
    df_ratio.loc['Before_filter'] = df_ratio[df_ratio.columns].sum(axis=0)
    df_ratio.reset_index(level=0, inplace=True)
    df_ratio = df_ratio[df_ratio.Locus.str.startswith('B')]
    df_ratio = df_ratio.transpose()
    df_ratio.reset_index(level=0, inplace=True)
    df_ratio.columns = df_ratio.iloc[0]
    df_ratio = df_ratio.drop(df_ratio.index[0])
    df_ratio[['ratio','ratio_number','exp']] =df_ratio['Locus'].str.split('_', expand=True)
    df_ratio = df_ratio.drop(['ratio','Locus','exp'], axis=1)

    #df_ratio is later combined with data after filtering and saved

    #write removed loci to file
    resultFyle = open((args.output_prefix + "e5_deleted_loci_p%s_%s_%s.csv" % (args.pool, args.filter_1, args.filter_2)),'wb')
    wr = csv.writer(resultFyle, dialect='excel')
    for item in getridof:
        item = ''.join(item)
        wr.writerow(item)

    #Here i remove the failed loci
    df = df.drop(getridof,axis=0)
    print "number of rows after: ",len(df.index)

    df_ratio2 = df.filter(regex='^ratio', axis=1)
    df_ratio2.loc['After_filter'] = df_ratio2[df_ratio2.columns].sum(axis=0)
    df_ratio2.reset_index(level=0, inplace=True)
    df_ratio2 = df_ratio2[df_ratio2.Locus.str.startswith('A')]
    df_ratio2 = df_ratio2.transpose()
    df_ratio2.reset_index(level=0, inplace=True)
    df_ratio2.columns = df_ratio2.iloc[0]
    df_ratio2 = df_ratio2.drop(df_ratio2.index[0])
    df_ratio2[['ratio','ratio_number2','exp']] =df_ratio2['Locus'].str.split('_', expand=True)
    df_ratio2 = df_ratio2.drop(['ratio','Locus'], axis=1)

    df_ratio_merged = pd.concat([df_ratio,df_ratio2], axis=1)

    if args.extra:
        csv_name = (args.output_prefix + 'e6_ratio_before_and_after_loci_filter_p%s_%s_%s.csv'%(args.pool, args.filter_1, args.filter_2))
        df_ratio_merged.to_csv(csv_name,sep='\t' , decimal =',')
        print "written : %s" % csv_name

    df_jenamono2 = df.filter(regex='^Locus|^jenamono',axis=1)
    df_jenamono2.sort_values(by='jenamono1',ascending=False)


    if args.extra:
        csv_name = (args.output_prefix + 'e7_jenamono_afterlocifilter_%s_%s.csv'%(args.filter_1, args.filter_2))
        df_jenamono2.to_csv(csv_name,sep='\t' , decimal =',')
        print "written : %s" % csv_name

    SAMPLE_NUM = len(df.columns)-2
    Pool = int(args.pool)
    Pool1_del = ['jenamono2', 'jenamono4', 'jenamono6', 'jenamono7', 'jenamono8']
    Pool2_del = ['jenamono3', 'jenamono5', 'jenamono9', 'jenamono10', 'jenamono11']

    #STD1 en 15 en 16 een eruit gewipt
    #renaming standaard to STD so that they are used for analysis

    Poolstd1_name = ['standaard_1_jena2016','standaard_2_jena2016','standaard_3_jena2016','standaard_4_jena2016',
                     'standaard_5_jena2016','standaard_6_jena2016','standaard_7_jena2016','standaard_8_jena2016',
                     'standaard_9_jena2016','standaard_10_jena2016']
    '''Standaard_1 is not replaced see below'''
    Poolstd1_replaCE = ['Standaard_1_jena2016','STD_2_jena2016','STD_3_jena2016','STD_4_jena2016','STD_5_jena2016',
                        'STD_6_jena2016','STD_7_jena2016','STD_8_jena2016','STD_9_jena2016','STD_10_jena2016']

    Poolstd2_name = ['standaard_11_jena2016','standaard_12_jena2016','standaard_13_jena2016','standaard_14_jena2016',
                     'standaard_15_jena2016','standaard_16_jena2016','standaard_17_jena2016','standaard_18_jena2016',
                     'standaard_19_jena2016','standaard_20_jena2016']
    Poolstd2_replaCE = ['STD_11_jena2016','STD_12_jena2016','STD_13_jena2016','STD_14_jena2016','Standaard_15_jena2016',
                        'STD_16_jena2016','STD_17_jena2016','Standaard_18_jena2016','STD_19_jena2016','STD_20_jena2016']

    if Pool ==1:
        Poollist = Pool1_del
        pooldict = {i:j for i,j in zip(Poolstd1_name,Poolstd1_replaCE)}
        df.rename(columns=pooldict, inplace=True)

    elif Pool ==2:
        Poollist = Pool2_del
        pooldict = {i:j for i,j in zip(Poolstd2_name,Poolstd2_replaCE)}
        df.rename(columns=pooldict, inplace=True)
    else:
        print "please use args.pool"
        print "default = 0"

    print

    #Here i filter for Loci of species in the Pool
    for p in Poollist:
        df = df[df.Species != p]

    df_NOS= df[df.columns[pd.Series(df.columns).str.startswith('STD')]]
    df_NOR= df[df.columns[2:]]

    #HET AANTAL SOORTEN
    SOORTLIJST = df.drop_duplicates(['Species'])
    print "REF       The",len(SOORTLIJST), "species are :"
    soortenlijst = list(SOORTLIJST['Species'])
    print soortenlijst

    Species_dict = {'jenamono1':'Plantago lanceolate mono 1',
                    'jenamono2':'Ranunculus acris mono 2',
                    'jenamono3':'Knautia arvensis mono 3',
                    'jenamono4':'Geranium pratense mono 4',
                    'jenamono5':'Centaurea jacea mono 5',
                    'jenamono6':'Dactylis glomerata mono 6',
                    'jenamono7':'Anthoxanthum odoratum mono 7',
                    'jenamono8':'Holcus lanatus mono 8',
                    'jenamono9':'Festuca rubra mono 9',
                    'jenamono10':'Avenula pubescens mono 10',
                    'jenamono11':'Poa pratensis mono 11',
                    'jenamono12':'Leucanthemum vulgare mono 12',
                    'jenamono13':'Phleum pratense mono 13'}

    new_soortenlijst =[]
    for item in soortenlijst:
        realname = Species_dict.get('%s'%item)
        new_soortenlijst.append(realname)

    print new_soortenlijst
    print
    NOS=len(df_NOS.columns)
    samples=len(df_NOR.columns)
    print 'The number of standards is               : ',NOS
    print 'The number of all samples (including standards) is         : ',samples
    print

    #here i replace the mono names to species names:

    df = df.replace(to_replace=soortenlijst,value=new_soortenlijst,regex=False )

    if args.extra:
        csv_name = (args.output_prefix + 'e8_tussen_pool%s.csv_%s_%s.csv'%(args.pool, args.filter_1, args.filter_2))
        df.to_csv(csv_name,sep='\t' , decimal =',')
        print "written : %s" % csv_name


    #Here the totals per species are culculated

    df = df.groupby('Species').sum()

    df = df.transpose()

    df['Total'] =df.sum(axis=1)

    if args.extra:
        csv_name = (args.output_prefix + 'e9_tussen_sum_p%s_%s_%s.csv'%(args.pool, args.filter_1, args.filter_2))
        df.to_csv(csv_name,sep='\t' , decimal =',')
        print "written : %s" % csv_name

    #Filter minimum number of reads per sample: need to be higher later?

    df = df[df['Total'] > int(args.filter_3)]

    #Here the reads are corrected for the total number of reads per sample
    for soort in new_soortenlijst:
        df['%s'%(soort+'_rel')]= df['%s'%soort] / df['Total']

    if args.extra:
        csv_name = (args.output_prefix + 'e10_tussen_sum_REL_p%s_%s_%s.csv'%(args.pool, args.filter_1, args.filter_2))
        df_e10 = df.transpose()
        df_e10 = df_e10.filter(regex=("^ratio_"))
        df_e10 = df_e10.filter(regex=("_jena2016$"))
        df_e10 = df_e10.rename(columns=lambda x: re.sub('ratio_','',x))
        df_e10 = df_e10.rename(columns=lambda x: re.sub('_jena2016', '', x))
        df_e10.columns = df_e10.columns.astype(int)
        df_e10 = df_e10.sort_index(axis=1)
        df_e10 = df_e10.transpose()
        df_e10['sample'] = df_e10.index
        if Pool == 1:
            df_e10 = df_e10[df_e10["sample"] <57]
        elif Pool == 2:
            df_e10 = df_e10[df_e10["sample"] >56]
        elif Pool == 0:
            print "no pool information given"

        df_e10.to_csv(csv_name,sep='\t' , decimal =',')
        print "written : %s" % csv_name

    STD_data = df.filter(regex=("^STD_"),axis=0)
    STD_data = STD_data.filter(regex=("rel$"),axis=1)

    standaardwaarden =list(STD_data.mean(axis=0))
    STD_dict = dict(zip(new_soortenlijst,standaardwaarden))

    #Here the standard corrected ratio values are calculated
    for key, value in STD_dict.iteritems():
        df['%s' % (key + '_stdcor')] = df['%s' %(key + '_rel')] / value

    STD_TMP = df.filter(regex=("stdcor$"))
    df['stdcor_total'] =STD_TMP.sum(axis=1)

    #Here the fractions (PCT) are calculated
    for soort in new_soortenlijst:
        df['%s'%(soort+'_PCT')]= df['%s_stdcor'%soort] / df['stdcor_total']
    STD_data2 = df.filter(regex=("PCT$"))

    # Here i filter for ratio samples only:

    STD_data2 = STD_data2.transpose()
    STD_data3 = STD_data2.filter(regex=("^ratio_"))
    STD_data4 = STD_data3.filter(regex=("_jena2016$"))
    STD_data4 = STD_data4.rename(columns=lambda x: re.sub('ratio_','',x))
    STD_data4 = STD_data4.rename(columns=lambda x: re.sub('_jena2016', '', x))
    STD_data4.columns = STD_data4.columns.astype(int)
    STD_data4 = STD_data4.sort_index(axis=1)

    STD_data4 = STD_data4.transpose()

    STD_data4['sample'] = STD_data4.index

    if Pool == 1:
        STD_data4 = STD_data4[STD_data4["sample"] <57]
    elif Pool == 2:
        STD_data4 = STD_data4[STD_data4["sample"] >56]
    elif Pool == 0:
        print "no pool information given"
    csv_name = (args.output_prefix + 'FINAL' + args.output_suffix)

    STD_data4.to_csv(csv_name,sep='\t' , decimal =',')
    print "written %s%s" % (args.output_prefix, args.output_suffix)
    #TODO: insert filter and pool options automaticly
    return args

##TODO: maybe filter out or make more files directly ? example: STD file / Ratio file / mono file / experiment file

def main():
    """main function loop"""
    #1 get command line arguments
    args = parse_args()

    args = calculus(args)
    print "time finished : ", datetime.datetime.now()

if __name__ == '__main__':
    main()


