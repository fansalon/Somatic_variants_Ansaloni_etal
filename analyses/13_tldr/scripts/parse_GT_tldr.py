import sys
import argparse
import pandas as pd
import collections
import numpy as np

##########################################################################################################################################################################################
##########################################################################################################################################################################################
def getArgs():
  parser = argparse.ArgumentParser()

  # create argument list
  parser.add_argument('--input',type=str, help='input file, resulting from tldr [required]', required=True)
  parser.add_argument('--refs',type=str, help='fai index file from TE consensus fasta sequences [required]', required=True)

  # create arguments
  arg = parser.parse_args()

  return arg




def parse_GT(dataf,col):

  dictio = {}
  dc_keys = ["WT_2_CTX","WT_2_HIP","WT_2_KID",
                    "WT_3_CTX","WT_3_HIP","WT_3_KID",
                    "WT_4_CTX","WT_4_HIP","WT_4_KID",
                    "WT_5_CTX","WT_5_HIP","WT_5_KID",
                    "WT_6_CTX","WT_6_HIP","WT_6_KID"]

  # open file
  with open(dataf) as inpudf:
    header = inpudf.readline()
    for line in inpudf:
      line = line.strip().split("\t")

      ####### get column of interest

      # 1 split by ,
      sup_reads = line[col].split(",")
      dictio_tmp = {}

      # sometimes sup_reads can be NA
      if sup_reads[0] == 'NA':
        dictio_tmp = dictio_tmp
      else:
        # 2 for each sample, split by | getting sample name and number of reads --> write to dictio
        for i in range(0,len(sup_reads)):
          sample       = sup_reads[i].split("|")[0]
          sample_reads = sup_reads[i].split("|")[1]
          # add to list
          dictio_tmp[sample] = sample_reads


      # 3 iterate across each known keys / samples and add sup reads for each sample for the MEI in analysis
      keyerror_c = 0
      for k in dc_keys:

        # get supporting reads from dictio tmp of the sample k in analysis
        # if key/sample not found, it means that the supporting reads are 0
        try:
          a = int(dictio_tmp[k])
        except KeyError:
          a = 0
          keyerror_c += 1

        # add number read value to original dictio
        if k in dictio.keys():
          dictio[k].append(a)
        else:
          dictio[k] = [a]

  return dictio


def getMEIType(dataf):

  # tissues and mice
  tissues = [ "CTX", "HIP", "KID"]
  mice = [ "2", "3", "4", "5", "6" ]

  # Step 1
  # count, for each mouse, the number of tissues where the MEI is found
  new_df = pd.DataFrame(dataf["ID"])
  for i in mice:
      samples_interest = []
      # get samples of interest
      for tx in tissues:
        samples_interest.append("WT_" + i + "_" + tx + "_ALT")

      # REF reads for the mouse in analysis
      d1 = dataf[samples_interest]
      # count in how many tissues, for each mouse, the MEI is supported by at least 1 read
      new_df["n_tx_mouse_"+str(i)] = (d1>=1).astype(int).sum(axis=1)

  # Step 2
  # if the MEI is found in 3 tissues in at least one mouse is classified as common, otherwise as tissue specific
  new_df.index = new_df["ID"]
  new_df = new_df.drop(["ID"], axis=1)
  final_type = []

  # read file row by row (i.e., MEI by MI)
  for i in range(0,len(new_df.iloc[:,0])):
    # get row
    row = new_df.iloc[i,:]

    # count in how many mice the MEI is found in 3/3 tissues
    common = (row>=3).astype(int).sum()

    # if the MEI is found in all the 3 tissues in 0 mice --> tissue-specific
    if common == 0:
      type = "Tissue-specific"
      # if tissue specific, count in how many mice the MEI is found
      tspec = (row>=1).astype(int).sum()
      # if found in 1 moue only --> tissue-specific private
      if tspec == 1:
        type2 = "private"
      # if found in all the mice --> tissue-specific recurrent
      elif tspec == 5:
        type2 = "recurrent"
      else:
        type2 = "unclassified"

    # if the MEI is found in all the 3 tissues, in at least one mouse --> common
    elif common >= 1:
      type = "Common"
      # if the MEI is found in all the 3 tissues, in just one mouse --> common private
      if common == 1:
        type2 = "private"
      # if the MEI is found in all the 3 tissues, in all the 5 mice --> common recurrent
      elif common == 5:
        type2 = "recurrent"
      else:
        type2 = "unclassified"

    # assign classification
    final_type.append(type+"_"+type2)

  new_df["Type"] = final_type
  new_df["ID"] = new_df.index
  new_df = new_df.reset_index(drop=True)
  return(new_df[["ID","Type"]])


# define in MEI is full length or not - it is considered full lenght if no more than 5% bp of the TE consensus length are missing from both 5' and 3'
def MEIfullLength(df,refs_l):

  # read file
  refs_l = pd.read_csv(refs_l,sep='\t',header=None)
  refs_l.columns = ["TE","len","o1","o2","o3"]

  # format TE col
  refs_l[["TE", "Subfamily"]] = refs_l["TE"].str.split(":", n=1, expand=True)
  # drop cols
  refs_l = refs_l.drop(["TE","o1","o2","o3"], axis=1)
  # add to file
  df = pd.merge(df,refs_l,on="Subfamily")

  ### define if MEI is full length or not
  # 1 calc differences between TE len and end coord
  df["missing_end"] = abs(df["EndTE"] - df["len"])

  # 2 calc "tolerance" in length difference, say 5%
  df["tolerance"] = 5 * df["len"] / 100

  # 3 full length if both StartTE and missing_end are <= tolerance
  df["full_length"] = (df["StartTE"] <= df["tolerance"]) & (df["missing_end"] <= df["tolerance"])

  ### clean-up
  df = df.drop(["missing_end","len","tolerance"], axis=1)

  #print(df[["Subfamily","StartTE","EndTE","len","missing_end","tolerance","full_length"]].to_string())

  return df



def main():

  # get args
  args = getArgs()
  df = args.input
  pd_df = pd.read_csv(df,sep="\t",header=0)

  # drop some columns
  pd_df = pd_df.drop(["UUID","Consensus"], axis=1)
  # generate ID
  pd_df["ID"] = pd_df["Chrom"] + ":" + pd_df["Start"].apply(str) + "-" + pd_df["End"].apply(str)


  # parse the columns containing supporting and empppty reads info
  # the columns are formatted as follow: "WT_2_CTX|17,WT_2_HIP|16,WT_3_CTX|19"

  # 1 - supporting reads
  #print("Parsing number of supporting reads")
  sup_reads_dictio = parse_GT(df,17)

  # 2 - empty reads
  #print("Parsing number of empty reads")
  neg_reads_dictio = parse_GT(df,18)
  # convert to dataf
  sup_reads_df = pd.DataFrame.from_dict(sup_reads_dictio)
  sup_reads_df.index = pd_df["ID"]
  neg_reads_df = pd.DataFrame.from_dict(neg_reads_dictio)
  neg_reads_df.index = pd_df["ID"]

  # 3 - calc maf
  maf_df = round( sup_reads_df/ (neg_reads_df+sup_reads_df), 3)

  # 4 - is MEI somatic?
  somatic_df = (maf_df<0.05).astype(int).sum(axis=1)

  # 5 - calc mean maf - do it after having counted the somatic
  maf_df["mean"] = round(maf_df.mean(skipna=True,axis=1),3)
  # add is.somatic to mean maf
  maf_df["NumSomaticSamples"] = somatic_df

  # 6 - rename columns
  sup_reads_df.columns = sup_reads_df.columns+"_ALT"
  neg_reads_df.columns = neg_reads_df.columns+"_REF"
  maf_df.columns = maf_df.columns+"_MAF"

  # 7 - get ID column and reset index, avoiding writing it to new column
  maf_df["ID"]       = maf_df.index
  sup_reads_df["ID"] = sup_reads_df.index
  neg_reads_df["ID"] = neg_reads_df.index
  sup_reads_df       = sup_reads_df.reset_index(drop=True)
  neg_reads_df       = neg_reads_df.reset_index(drop=True)
  maf_df             = maf_df.reset_index(drop=True)

  # 8 - classify MEI as tissue specific, common, private, recurrent
  mei_type = getMEIType(sup_reads_df)

  # 9 - define in MEI is full-length or not
  pd_df = MEIfullLength(pd_df,args.refs)


  ### Prepare for output

  # get file 1, MEI information
  # 1 - convert tsd to bool
  pd_df["TSD"] = pd_df["TSD"].notna()
  # 2 - get file 1 --> info on MEI + mean MAF + is.somatic + type
  f1 = pd.merge(pd_df, maf_df, on = "ID")                                              # merge 1
  ff = pd.merge(f1,mei_type, on = "ID")                                                # merge, add also type (common/recurrent, etc)
  ff = ff.rename(columns={"NumSomaticSamples_MAF": "NumSomaticSamples"})               # rename col (_MAF was added to all names in some previous step)
  ff["is.somatic"] = ff["NumSomaticSamples"]==ff["NumSamples"]                         # define if somatic (number of samples with maf<0.05 equals to tot numb samples)
  ff["TE"] = ff["Family"] + ":" + ff["Subfamily"]                                      # TE ID
  ff2 = ff[ ["Chrom","Start","End","ID","Strand","TE","StartTE","EndTE","LengthIns","Inversion","NonRef","TSD","Filter","NumSamples","NumSomaticSamples","is.somatic","mean_MAF","Type","full_length"] ]

  ff2.to_csv(df+".bed",sep='\t',index=False)

  #print("Number of somatic events (MEI has MAF<0.05 in all the samples):")
  #print("\t%s" % str( (ff2["is.somatic"]==True).astype(int).sum()))
  #print("\t%s" % collections.Counter(ff2[ff2["is.somatic"]==True]["Type"]))


  # get file 2 - supporting reads ALT, REF and maf
  # 1 - merge ALT and REF
  ff3_1 = pd.merge(sup_reads_df, neg_reads_df, on = "ID")
  # add maf
  maf_df = maf_df.drop(["NumSomaticSamples_MAF","mean_MAF"], axis=1)
  ff3_2 = pd.merge(ff3_1,maf_df, on = "ID")
  # add type
  ff3 = pd.merge(ff3_2,mei_type, on = "ID")
  # move ID col to first place
  id_col = ff3.pop("ID")
  ff3.insert(0,"ID",id_col)
  ff3.to_csv(df.replace(".txt","")+"_reads_summary.txt",sep='\t',index=False)















if __name__ == '__main__':
  main()
