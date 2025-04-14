#### Setup ####
mkdir -p ~/project/sh2b3/data
mkdir -p ~/project/sh2b3/ref
mkdir -p ~/project/sh2b3/results

s3fs dagr-sh2b3-results ~/project/sh2b3/data -o passwd_file=~/.passwd-s3fs \
  -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007
  
s3fs mouse-ref ~/project/sh2b3/ref -o passwd_file=~/.passwd-s3fs \
  -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007
  
#### Count exons #### 
for bam in ~/project/sh2b3/data/3_bam_filter/*bam;
do
    name=`echo "$(basename $bam)" | sed 's/_Aligned.sortedByCoord.filter.bam//'`
    echo $name
    
    featureCounts -T 10 -g exon_id -t exon -p \
        -a ~/project/sh2b3/ref/release106/STARref/Mus_musculus.GRCm39.106.gtf \
        -o ~/project/sh2b3/results/"$name"_feature_counts_exons.tsv $bam

done

#### Combine samples ####
python

import pandas as pd
import glob
import os

## List all counts files
all_count_file = glob.glob('/home/ec2-user//project/sh2b3/results/*_exons.tsv')

## Combine all files
### Blank df to hold results
all_count = pd.DataFrame(index=[])

for f in all_count_file:
    print(f)
    ### Read in tsv
    temp = pd.read_csv(f, delimiter="\t", skiprows=1, index_col='Geneid')
    ### Remove gene info columns
    temp = temp.drop(['Chr','Start','End','Strand','Length'], axis=1)
    ### Combine
    all_count = all_count.join(temp, how='outer')

## Clean sample names
all_count.columns = all_count.columns.str.replace('/home/ec2-user/project/sh2b3/data/3_bam_filter/', '', regex=True)
all_count.columns = all_count.columns.str.replace('_Aligned.sortedByCoord.filter.bam', '', regex=True)

## Save
all_count.to_csv("/home/ec2-user/project/sh2b3/results/combined_feature_counts_exon.tsv", index=True, encoding='utf-8-sig', sep="\t")

quit()

#### Save ####
aws s3 sync ~/project/sh2b3/results s3://dagr-sh2b3-results/6_exon_level/
