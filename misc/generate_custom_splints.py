import sys
import argparse
import mappy

parser = argparse.ArgumentParser(description='Makes R2C2 splints from provided PCR primers and 16 predefined backbones. The splints will be able to circularize an amplicon/library amplified with the provided PCR primers',add_help=True)

parser.add_argument('--Primer1', '-p1', type=str, action='store',
                      help='First primer (should be at least 20nt)')
parser.add_argument('--Primer2', '-p2', type=str, action='store',
                      help='Second primer (should be at least 20nt)')

args = parser.parse_args()

backbones={}
backbones['UMI_Splint_1']='TGAGGCTGATGAGTTCCATANNNNNTATATNNNNNATCACTACTTAGTTTTTTGATAGCTTCAAGCCAGAGTTGTCTTTTTCTCTTTGCTGGCAGTAAAAGTATTGTGTACCTTTTGCTGGGTCAGGTTGTTCTTTAGGAGGAGTAAAAGGATCAAATGCACTAANNNNNTATATNNNNNGCGATCGAAAATATCCCTTT'
backbones['UMI_Splint_2']='TGCCGGTTGGGTATCAATAANNNNNTATATNNNNNATTGCCTTTATTCTATCTACTTAGTTTTGGCGATGTAGTCTACCTATCCTGATGCTGAATAAAGGCCTTCATGATATAGGCGTCAGAGAGTTGTAAGTTCTTGTTACGGGTGAAAAATTTTCGATGGCAGNNNNNTATATNNNNNCGTGATCCTAGAACCTAATT'
backbones['UMI_Splint_3']='GTCGTGATCAAACATTGGGTNNNNNTATATNNNNNTAAAAGTTTTCTGTGTCCATTACGTTTTTTGGAGACGGTCTCAACTATTCTTAATCTCGGCGAACTTGGGTTAATTTCCTAGAGATTTCGGGGACACTTTTCAGAGGAGGTAAGGTTTCGAACTGAATATNNNNNTATATNNNNNAGCTGTAATACGACTCACTT'
backbones['UMI_Splint_4']='TAGTGTAAGGTAGCATCCGTNNNNNTATATNNNNNTACGATCATATTATGTTCGTCTGTCTTTAACGGGCTATACTTCTTCTTTAGGAAAGGTCTGAACCTTCAATCTGGGTTCTACTTAAGGTTATGGAATCTTAAGGGTAGTAGTGGGGTCTAATACACAGGTNNNNNTATATNNNNNCGCCTTTACTCAAATAGTAG'
backbones['UMI_Splint_5']='GTGGCTCATGTCAATGATAGNNNNNTATATNNNNNTTCTGATTATATTTTGCTTTTGTTACTCCTTAGTCATTACTGGTCGAAGGATTCGGCCTAAAACAGTAACTACTTAAATATCCTCATGTAGTGGAGTTATAGTGTTGACGAGAAGTGGTGCTGATGGTTCNNNNNTATATNNNNNCTCTTCTCCGAAAATGTGAA'
backbones['UMI_Splint_6']='GAGTTTAGCACATGACTGGTNNNNNTATATNNNNNACGTCTCTGAACTTTTACTCTGCTTATTTATCTAGTTATTTAGCATGCGTAGATGGAGCTGATTACTTTAGATGTTTTAGGCCGGTTACAACTCTTAGTCGATTAAACTGGAGGATTGGAATAGTGGCTANNNNNTATATNNNNNCTAATATAAGCGTTCCCTAG'
backbones['UMI_Splint_7']='CCGAGAGGTGTATGCTTATANNNNNTATATNNNNNTGGACTTTTACCGTCTTTAATGTATCTGGCATATGTAGATTTAAGTATATGCTCCTGGTTCACCATGTTAGAGAGTCAATGAAGATCTATATAGAATCGTCTGGTATGTGTAGTGTAGTTCCTGCATCTGNNNNNTATATNNNNNCCTCTGACATACAAGTTGAT'
backbones['UMI_Splint_8']='CTGACATTTCGGTGGAGAATNNNNNTATATNNNNNTTCTCAACTGTAAAAGTCTCCCTGGATATATTTGTGTTTATGCTGATATTGGCATCCATGTTTGACGGAGGATTATCAGGTAGGTAAATTACTTCATTTGGAGATGAGGTGGTTGTACATTAACTTCCCTNNNNNTATATNNNNNAGCCTTCAACTGAAAGTTCT'
backbones['UMI_Splint_9']='ATCTGAAACTGTTATCGGGGNNNNNTATATNNNNNAAGATGATACCCTCTTTTCACTTTGTGTATGTAGAAATGGGTTCTGACTTTAAGTCTTTCTGTACCATGGTTGACACCGTCTATGTCTATCAGGGTTCAAAATGGTGTATTAGTAGAACGTGTATTAGGTNNNNNTATATNNNNNGTTACAGAATTCACTCCGAT'
backbones['UMI_Splint_10']='GGCGTGGCAAATTTATTCAGNNNNNTATATNNNNNATTAAACAATTTACAATCTGTGTCGGCAGATTGCTTTGGTGTTACCCCTCGATTTGTTCAGTTTATCAGAGCAGTTTTTGGGTGTATCGAAAATCGGGATTTTGAATACCTTGTTAGACATTGTAAGGCTNNNNNTATATNNNNNACCAGTATAGTCCCATTGTA'
backbones['UMI_Splint_11']='CTGATACGTGTTGAGCGATANNNNNTATATNNNNNATTCTTATGTTCGCAATCTTCATGGACTTTCTATCTGTCTAGTATATGGTCAGACTTAGACGATTGAGACGTGCTAAGTGTGATTTGTCTTGACAAGGCCTAATGGACATGAGGTTCTTGTATATTATGANNNNNTATATNNNNNCACCAATTAGACACTGTTTG'
backbones['UMI_Splint_12']='CACAGTATGTGTGTGCAAGTNNNNNTATATNNNNNCATACTGTGCCATTTTATCATCGGCAATGGTAAACTCTTTAAGTGGTTGTACGGTATTTTCTTATCTTAGCTAAATTACCTCGGTTAACTGTATATCGATCTGTAAGATAAGGAGAGCTGGGTGTTTGTGNNNNNTATATNNNNNGACTTGCTTCATATGCCAAA'
backbones['UMI_Splint_13']='AGGGTCCGATAATTGATCTGNNNNNTATATNNNNNATGCACGTATTGCCATTGATGGATTTTATCTGTATTCCGTACGTTTCCTTAGAACTTATGTTAGACTGGTCTGTATTGAAAGCTAAGGTTCGAATGGTGATAGTGGTTAACAATCTTTACGTGATCTACGNNNNNTATATNNNNNGTCATCAAGCATCACAGTTT'
backbones['UMI_Splint_14']='TTCTAGACTGGGAATTGGACNNNNNTATATNNNNNACTTCACTTGGACTTACTTTGTAAGAGGGGTATAATCTATCTATGTGTTTTCTTGCTTACCGACTATACTAGTGTGACATCTAGTATTGTTGTCAACGAGATTGATCAATTTCGAGCTGGGTGGGATTAANNNNNTATATNNNNNCCTCAGTTAACAATGCATTG'
backbones['UMI_Splint_15']='GGTCTTAGAATGGAACCTGTNNNNNTATATNNNNNCATGGATATAGTATTTATCATATTGAATCAGTTTTCTGTTACCAACGGTTCTTTCCTTTCGCGGGACCGGAATGAGGATGTTGGGATGTTGAGATCTGATGTATCATGTAGACTCTAGCTCTTTAAAATTNNNNNTATATNNNNNTACGAACTATACTGATCGTC'
backbones['UMI_Splint_16']='TATTCCTAGGTGGGAGTACANNNNNTATATNNNNNTCTATATAGTTCTTTGATTATAATATCGATTGCGGTTTACCGCATCAGGTCTCGCTATGATCTAGTTATACGGACAAACTCGTGAGTGACGATTTGTTGTGAATTCTTTTTGACATCGGGGTTTGAGAAANNNNNTATATNNNNNATTCTCTGATACCGGTAACA'

Oligo1={}
Oligo2={}
Oligo1['UMI_Splint_1']='TGAGGCTGATGAGTTCCATA NNNNNTATATNNNNN ATCACTACTTAGTTTTTTGATAGCTTCAAGCCAGAGTTGTCTTTTTCTCTTTGCTGGCAGTAAAAG'
Oligo2['UMI_Splint_1']='AAAGGGATATTTTCGATCGC NNNNNATATANNNNN TTAGTGCATTTGATCCTTTTACTCCTCCTAAAGAACAACCTGACCCAGCAAAAGGTACACAATA CTTTTACTGCCAGCAAAGAG'
Oligo1['UMI_Splint_2']='TGCCGGTTGGGTATCAATAA NNNNNTATATNNNNN ATTGCCTTTATTCTATCTACTTAGTTTTGGCGATGTAGTCTACCTATCCTGATGCTGAATAAAGGC'
Oligo2['UMI_Splint_2']='AATTAGGTTCTAGGATCACG NNNNNATATANNNNN CTGCCATCGAAAATTTTTCACCCGTAACAAGAACTTACAACTCTCTGACGCCTATATCATGAAG GCCTTTATTCAGCATCAGGA'
Oligo1['UMI_Splint_3']='GTCGTGATCAAACATTGGGT NNNNNTATATNNNNN TAAAAGTTTTCTGTGTCCATTACGTTTTTTGGAGACGGTCTCAACTATTCTTAATCTCGGCGAACT'
Oligo2['UMI_Splint_3']='AAGTGAGTCGTATTACAGCT NNNNNATATANNNNN ATATTCAGTTCGAAACCTTACCTCCTCTGAAAAGTGTCCCCGAAATCTCTAGGAAATTAACCCA AGTTCGCCGAGATTAAGAAT'
Oligo1['UMI_Splint_4']='TAGTGTAAGGTAGCATCCGT NNNNNTATATNNNNN TACGATCATATTATGTTCGTCTGTCTTTAACGGGCTATACTTCTTCTTTAGGAAAGGTCTGAACCT'
Oligo2['UMI_Splint_4']='CTACTATTTGAGTAAAGGCG NNNNNATATANNNNN ACCTGTGTATTAGACCCCACTACTACCCTTAAGATTCCATAACCTTAAGTAGAACCCAGATTGA AGGTTCAGACCTTTCCTAAA'
Oligo1['UMI_Splint_5']='GTGGCTCATGTCAATGATAG NNNNNTATATNNNNN TTCTGATTATATTTTGCTTTTGTTACTCCTTAGTCATTACTGGTCGAAGGATTCGGCCTAAAACAG'
Oligo2['UMI_Splint_5']='TTCACATTTTCGGAGAAGAG NNNNNATATANNNNN GAACCATCAGCACCACTTCTCGTCAACACTATAACTCCACTACATGAGGATATTTAAGTAGTTA CTGTTTTAGGCCGAATCCTT'
Oligo1['UMI_Splint_6']='GAGTTTAGCACATGACTGGT NNNNNTATATNNNNN ACGTCTCTGAACTTTTACTCTGCTTATTTATCTAGTTATTTAGCATGCGTAGATGGAGCTGATTAC'
Oligo2['UMI_Splint_6']='CTAGGGAACGCTTATATTAG NNNNNATATANNNNN TAGCCACTATTCCAATCCTCCAGTTTAATCGACTAAGAGTTGTAACCGGCCTAAAACATCTAAA GTAATCAGCTCCATCTACGC'
Oligo1['UMI_Splint_7']='CCGAGAGGTGTATGCTTATA NNNNNTATATNNNNN TGGACTTTTACCGTCTTTAATGTATCTGGCATATGTAGATTTAAGTATATGCTCCTGGTTCACCAT'
Oligo2['UMI_Splint_7']='ATCAACTTGTATGTCAGAGG NNNNNATATANNNNN CAGATGCAGGAACTACACTACACATACCAGACGATTCTATATAGATCTTCATTGACTCTCTAAC ATGGTGAACCAGGAGCATAT'
Oligo1['UMI_Splint_8']='CTGACATTTCGGTGGAGAAT NNNNNTATATNNNNN TTCTCAACTGTAAAAGTCTCCCTGGATATATTTGTGTTTATGCTGATATTGGCATCCATGTTTGAC'
Oligo2['UMI_Splint_8']='AGAACTTTCAGTTGAAGGCT NNNNNATATANNNNN AGGGAAGTTAATGTACAACCACCTCATCTCCAAATGAAGTAATTTACCTACCTGATAATCCTCC GTCAAACATGGATGCCAATA'
Oligo1['UMI_Splint_9']='ATCTGAAACTGTTATCGGGG NNNNNTATATNNNNN AAGATGATACCCTCTTTTCACTTTGTGTATGTAGAAATGGGTTCTGACTTTAAGTCTTTCTGTACC'
Oligo2['UMI_Splint_9']='ATCGGAGTGAATTCTGTAAC NNNNNATATANNNNN ACCTAATACACGTTCTACTAATACACCATTTTGAACCCTGATAGACATAGACGGTGTCAACCAT GGTACAGAAAGACTTAAAGT'
Oligo1['UMI_Splint_10']='GGCGTGGCAAATTTATTCAG NNNNNTATATNNNNN ATTAAACAATTTACAATCTGTGTCGGCAGATTGCTTTGGTGTTACCCCTCGATTTGTTCAGTTTAT'
Oligo2['UMI_Splint_10']='TACAATGGGACTATACTGGT NNNNNATATANNNNN AGCCTTACAATGTCTAACAAGGTATTCAAAATCCCGATTTTCGATACACCCAAAAACTGCTCTG ATAAACTGAACAAATCGAGG'
Oligo1['UMI_Splint_11']='CTGATACGTGTTGAGCGATA NNNNNTATATNNNNN ATTCTTATGTTCGCAATCTTCATGGACTTTCTATCTGTCTAGTATATGGTCAGACTTAGACGATTG'
Oligo2['UMI_Splint_11']='CAAACAGTGTCTAATTGGTG NNNNNATATANNNNN TCATAATATACAAGAACCTCATGTCCATTAGGCCTTGTCAAGACAAATCACACTTAGCACGTCT CAATCGTCTAAGTCTGACCA'
Oligo1['UMI_Splint_12']='CACAGTATGTGTGTGCAAGT NNNNNTATATNNNNN CATACTGTGCCATTTTATCATCGGCAATGGTAAACTCTTTAAGTGGTTGTACGGTATTTTCTTATC'
Oligo2['UMI_Splint_12']='TTTGGCATATGAAGCAAGTC NNNNNATATANNNNN CACAAACACCCAGCTCTCCTTATCTTACAGATCGATATACAGTTAACCGAGGTAATTTAGCTAA GATAAGAAAATACCGTACAA'
Oligo1['UMI_Splint_13']='AGGGTCCGATAATTGATCTG NNNNNTATATNNNNN ATGCACGTATTGCCATTGATGGATTTTATCTGTATTCCGTACGTTTCCTTAGAACTTATGTTAGAC'
Oligo2['UMI_Splint_13']='AAACTGTGATGCTTGATGAC NNNNNATATANNNNN CGTAGATCACGTAAAGATTGTTAACCACTATCACCATTCGAACCTTAGCTTTCAATACAGACCA GTCTAACATAAGTTCTAAGG'
Oligo1['UMI_Splint_14']='TTCTAGACTGGGAATTGGAC NNNNNTATATNNNNN ACTTCACTTGGACTTACTTTGTAAGAGGGGTATAATCTATCTATGTGTTTTCTTGCTTACCGACTA'
Oligo2['UMI_Splint_14']='CAATGCATTGTTAACTGAGG NNNNNATATANNNNN TTAATCCCACCCAGCTCGAAATTGATCAATCTCGTTGACAACAATACTAGATGTCACACTAGTA TAGTCGGTAAGCAAGAAAAC'
Oligo1['UMI_Splint_15']='GGTCTTAGAATGGAACCTGT NNNNNTATATNNNNN CATGGATATAGTATTTATCATATTGAATCAGTTTTCTGTTACCAACGGTTCTTTCCTTTCGCGGGA'
Oligo2['UMI_Splint_15']='GACGATCAGTATAGTTCGTA NNNNNATATANNNNN AATTTTAAAGAGCTAGAGTCTACATGATACATCAGATCTCAACATCCCAACATCCTCATTCCGG TCCCGCGAAAGGAAAGAACC'
Oligo1['UMI_Splint_16']='TATTCCTAGGTGGGAGTACA NNNNNTATATNNNNN TCTATATAGTTCTTTGATTATAATATCGATTGCGGTTTACCGCATCAGGTCTCGCTATGATCTAGT'
Oligo2['UMI_Splint_16']='TGTTACCGGTATCAGAGAAT NNNNNATATANNNNN TTTCTCAAACCCCGATGTCAAAAAGAATTCACAACAAATCGTCACTCACGAGTTTGTCCGTATA ACTAGATCATAGCGAGACCT'

primer1=args.Primer1
primer2=args.Primer2

print('\nCreating an R2C2 splint compatible with an amplicon amplified with the following provided primers\n')
print('Primer 1\t'+primer1)
print('Primer 2\t'+primer2)


for i in range(1,17,1):
    UMI='UMI_Splint_'+str(i)
    backbone=backbones[UMI]
    Oligo1seq=Oligo1[UMI]
    Oligo2seq=Oligo2[UMI]
    FullSplint=mappy.revcomp(primer1)+' '+backbone+' '+primer2
    FullOligo1=mappy.revcomp(primer1)+' '+Oligo1seq
    FullOligo2=mappy.revcomp(primer2)+' '+Oligo2seq
    print('\n\n### Splint'+' '+str(i)+' ###\n')

    print('>Splint '+str(i)+' Backbone')
    print(backbone)

    print('>Splint '+str(i)+' Full Sequence')
    print(FullSplint+'\n')

    print('Splint '+str(i)+' Oligo1\t'+FullOligo1)
    print('Splint '+str(i)+' Oligo2\t'+FullOligo2)
