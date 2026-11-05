import sys
import argparse
import statistics
import pandas as pd
from pysam import VariantFile
import scipy.stats as st
import statsmodels.stats.multitest as multi
import datetime
import numpy as np

##########################################################################################################################################################################################
##########################################################################################################################################################################################
def getArgs():
  parser = argparse.ArgumentParser()

  # create argument list
  parser.add_argument('--sample',type=str, help='sample file (e.g., something like WT_2_CTX [required]', required=True)
  parser.add_argument('--SR',type=str, help='SR input summary file [required]', required=True)
  parser.add_argument('--nanorepeat_dir',type=str, help='directory where the nanorepeat output are stored [required]', required=True)
  parser.add_argument('--outdir',type=str, help='directory where to write results [required]', required=True)
  parser.add_argument('--combined_vcf',type=str, help='path to combined vcf [required]', required=True)
  parser.add_argument('--zscore',type=str, help='path to combined vcf [required]', required=False,default=2.5)

  # create arguments
  arg = parser.parse_args()

  return arg


def recallGT(row_list,zscore_threshold):
  # if only one genotype is found, try to recall it to see if some somatic reads can be identified

  # convert list to pd
  read_name = []
  rep_size = []
  allele_id = []
  phasing = []
  for i in range(0,len(row_list)):
    # skip header element
    if i!=0:
      line_split = row_list[i].split(";")
      read_name.append(line_split[0])
      rep_size.append(line_split[1])
      allele_id.append(line_split[2])
      phasing.append(line_split[3])

  # iter across all the rep_size looking for outlier:
  # every round calculate mean and sd without the analysed value
  # then calc zscore
  rep_size = [float(i) for i in rep_size]

  if len(rep_size) >= 1:
    # calc IQR
    q1 = np.quantile(rep_size,q=0.25)
    q3 = np.quantile(rep_size,q=0.75)
    IQR = q3 - q1
    range_min = q1 - (1.5*IQR)
    range_max = q3 + (1.5*IQR)
    for j in range(0,len(rep_size)):
      if rep_size[j] < range_min or rep_size[j] > range_max:
        allele_id[j] = "-1"

  # there must be at least 3 elements, in order to calculate mean and sd on at least 2 since 1 element will be removed in each iteration
 # if len(rep_size) >= 3:
#    for j in range(0,len(rep_size)):
#      # select all the items except the analysed item
#      new_l = rep_size[0:j] + rep_size[j+1:]
#
#      # calc mean
#      media = np.mean(new_l)
#      # calc sd
#      sd = np.std(new_l)
#      if sd == 0:
#        sd = 0.00000000001
#      # calc zscore
#      zscore = abs((rep_size[j] - media) / sd)
#
#      # if zscore is > than threshold, convert genotype to "-1"
#      if zscore >= zscore_threshold:
#        allele_id[j] = "-1"

  # format as original list
  out_list = []
  # add "header"
  out_list.append('Read_Name;Read_Repeat_Size;Read_Allele_ID;PhasingConfidence')
  # add elements
  for i in range(0,len(read_name)):
    new_line = read_name[i]+";"+str(rep_size[i])+";"+allele_id[i]+";"+phasing[i]
    out_list.append(new_line)

  return out_list


# count how many genotypes are found
def countGT(ll):
  gt_list = []

  # read the list element by element
  for i in range(0, len(ll)):
    # discard the first element which contains a "header" and not data
    if i!=0:

      # get SR "genotype"
      gt_list.append(ll[i].split(";")[2])

  # make it unique
  gt_list_un = list(set(gt_list))

  return gt_list_un



def getSRlengthGT(ll):
  ln_germline = []
  ln_somatic = []

  # read the list element by element
  for i in range(0, len(ll)):
    # discard the first element which contains a "header" and not data
    if i!=0:

      # get SR length
      sr_l = ll[i].split(";")[1]

      # get SR "genotype"
      sr_gt = ll[i].split(";")[2]

      # append SR length to two separate lists (germline if GT is 1/2/... and somatic if GT is -1)
      if sr_gt == "1" or sr_gt == "2" or sr_gt == "3":
        ln_germline.append(float(sr_l))
      elif sr_gt == "-1":
        ln_somatic.append(float(sr_l))
      else:
        print("ERROR: detected genotype %s" % sr_gt)
        print(ll)
        sys.exit()

  return ln_germline, ln_somatic



def getSomaticLength(ref_len,lis):

  del_l = []
  ins_l = []
  # for each item in list, divide in del and ins
  for i in range(0,len(lis)):
    if lis[i] > ref_len:
      ins_l.append(lis[i])
    else:
      del_l.append(lis[i])

  # calc mean
  if len(del_l)>=2:
    mean_del = round( statistics.mean(del_l) - ref_len )
  elif len(del_l)==0:
    mean_del = "NA"
  else:
    mean_del = round( del_l[0] - ref_len )

  # do the same with SD
  if len(del_l)>=2:
    sd_del = round(statistics.stdev(del_l),2)
  elif len(del_l)==0:
    sd_del = "NA"
  else:
    sd_del = 0

  return mean_del, sd_del, len(del_l)


def parseVCF(pathvcf):
  chrom_ls = []
  start_ls = []
  len_ls = []
  svtype_ls = []

  # open vcf
  vcf = VariantFile(pathvcf)
  # read it line by line
  for rec in vcf:
    chrom_ls.append(rec.contig) # chrom
    start_ls.append(rec.pos) # start coord
    len_ls.append(round(rec.info["STDEV_LEN"],2)) # SVLEN stdev
    svtype_ls.append(rec.info["SVTYPE"]) # SVLEN stdev

  # write to list and the to pd df
  df_list = {
                "chr" : chrom_ls,
                "start" : start_ls,
                "SVLEN_sd" : len_ls,
                "SVTYPE" : svtype_ls
                }
  outvcf_df = pd.DataFrame(df_list)
  outvcf_df["id"] = outvcf_df["chr"] + "-" + outvcf_df["start"].apply(str) + "-" + outvcf_df["SVTYPE"]

  # select only columns of interest
  outvcf_df2 = outvcf_df.loc[:,["id","SVLEN_sd"]]

  return outvcf_df2

##########################################################################################################################################################################################
##########################################################################################################################################################################################



def main():

  # parse input args
  args = getArgs()

  # get args
  sample = args.sample
  sample2 = sample.split("_")
  sample2 = sample2[0] + "_" + sample2[1]

  # get input file
  df = args.nanorepeat_dir + "/" + sample + ".NanoRepeat_output.tsv"

  # get SR annotation file
  inp = pd.read_csv(args.SR,sep="\t",header=None,names=["chr","start","end","SR_consensus","id_sv"])
  inp["ID"] = inp["chr"] + ":" + inp["start"].apply(str) + "-" + inp["end"].apply(str)

  print("Analysing sample: %s" % sample)
  print("Analysing file: %s" % df)


  # start reading the file with nanorepeat results
  # read it line by line
  outrow = []
  validated = 0
  n_validated = 0
  tot = 0
  no_genot = 0
  single_genot = 0
  selected = 0
  multiple_genot = 0
  no_somatic = 0
  with open (df) as inpf:
    for line in inpf:

      tot += 1
      #print(tot)

      # get coord
      id = line.split("\t")[0] + ":" + line.split("\t")[1] + "-" + line.split("\t")[2]

      # get last column
      line_s = line.split("\t")[-1]

      # separate by |
      line_s2 = line_s.split("|")

      # re-call somatic reads - I am not entirely happy with the way nanorepeat calls the somatic reads
      # so I am here repeating the analysis. Basically, I iter across all the repeat length, calculating
      # every time mean and sd on all the values except the one in analysis. Then I calculate the zscore for the value in analysis
      # If abs(zscore) >2.5 --> the read is somatically supporting a variant
      line_s2 = recallGT(line_s2,args.zscore)

      # count genotypes (-1 is somatic, but then other genotypes are possibile)
      count_gt = countGT(line_s2)

      # continue with analysis only if 2 genotypes are found
      continue_analysis = False
      if len(count_gt) == 0: # no genotype --> discard
        no_genot += 1
      elif len(count_gt) > 2: # more than 2 genotypes --> discard
        multiple_genot += 1
      elif len(count_gt) == 1: # 1 genotype --> not validated
        single_genot += 1
      elif len(count_gt) == 2: # 2 genotypes
        if "-1" not in count_gt: # but no somatic --> meaning multiple genotypes!
          no_somatic += 1
        else: # 2 genotypes, and some somatic reads --> proceed!
          selected += 1
          continue_analysis = True


      if continue_analysis:
        #print(line_s2)

        # get SR length in two lists (one for somatic and one for germline)
        sr_len_germ, sr_len_som = getSRlengthGT(line_s2)

        # get SR mean length
        if len(sr_len_germ) == 1:
          sr_m_len = sr_len_germ[0]
        elif len(sr_len_germ) >=2:
          sr_m_len = round(statistics.mean(sr_len_germ))
        else:
          print("No genotype")
          sys.exit()

        # get: somatic SR DEL length, its sd and the number of reads supporting the somatic allele
        # sr_m_len_som_del is NA if no reads supporting the DEL are found, otherwise is the DEL len
        sr_m_len_som_del, sr_m_sd_som_del, del_reads = getSomaticLength(sr_m_len,sr_len_som)
        # get MAF of the somatic DEL
        tot_reads = len(line_s2)-1 # do not count header item
        maf = round(del_reads / tot_reads,3)


        # define if DEL is validated or not
        if sr_m_len_som_del == "NA":
          type = "not_validated"
          n_validated += 1
        elif sr_m_len_som_del != "NA":
          type = "DEL_validated"
          validated += 1
        else:
          print("What is this?")
          sys.exit()

        # write to list: id, INS len, INS reads, DEL len, DEL reads, DEL maf
        outrow.append([id, str(sr_m_len_som_del), str(sr_m_sd_som_del), del_reads, maf, type])


  # convert res to dataframe
  res = pd.DataFrame(outrow, columns = ['ID','DEL_len','DEL_len_sd','DEL_reads','DEL_MAF','type' ])


  # merge res with original input file
  merged = pd.merge(inp,res,on="ID",how="right")
  merged.drop("ID",inplace=True,axis=1)
  # modify output in a DEL-centric manner
  merged["chr_variant"] = merged["id_sv"].str.split("-",n=7).str[0]
  merged["start_variant"] = merged["id_sv"].str.split("-",n=7).str[1]
  merged["end_variant"] = merged["id_sv"].str.split("-",n=7).str[2]
  merged["SVTYPE"] = merged["id_sv"].str.split("-",n=7).str[3]
  merged["SVLEN"] = merged["id_sv"].str.split("-",n=7).str[4]
  merged["SV_class"] = merged["id_sv"].str.split("-",n=7).str[5]
  merged["Sample"] = sample
  merged["SR"] = merged["chr"] + "-" + merged["start"].apply(str) + "-" + merged["end"].apply(str) + "-" + merged["SR_consensus"]
  merged2 = merged.loc[:,["chr_variant","start_variant","end_variant","SVTYPE","SVLEN","SV_class","Sample","SR","DEL_len","DEL_len_sd","DEL_reads","DEL_MAF","type"]]
  # select only validated
  validated_df = merged2[merged2["type"]=="DEL_validated"]
  # add id for next step
  validated_df = validated_df.copy() # make a copy to avoid stupid pandas warning
  validated_df["id"] = validated_df["chr_variant"] + "-" + validated_df["start_variant"].apply(str) + "-" + validated_df["SVTYPE"]


  # now, add the sd of the SV length --> extract this from the combined vcf
  # this pd contains 2columns: col1: id "chr-start-SVTYPE"; col2: SVLEN sd
  vcf_data = parseVCF(args.combined_vcf)
  # merge
  validated_df2 = pd.merge(validated_df,vcf_data,on="id",how="left")
  validated_df2.drop("id",inplace=True,axis=1)


  # reorder and write to out
  validated_final = validated_df2.loc[:,["chr_variant","start_variant","end_variant","SVTYPE","SVLEN","SVLEN_sd","SV_class","Sample","SR","DEL_len","DEL_len_sd","DEL_reads","DEL_MAF","type"]]
  validated_final.to_csv(args.outdir+"/validation_res_"+sample+".txt",index=False,sep='\t')



  # create summary file
  summary_df = pd.DataFrame({
  "no_GT": ["0","1","2","2",">2"],
  "Count": [no_genot, single_genot, validated, n_validated, no_somatic+multiple_genot],
  "Supported_SVTYPE": ["NA","NA","DEL","INS","NA"],
  "Validated": [False,False,True,False,False]   })
  summary_df["Sample"] = sample
  # write to out
  summary_df.to_csv(args.outdir+"/summary_"+sample+".txt",index=False,sep='\t')

  validation_rate = validated / (single_genot+validated+n_validated) * 100
  print("\nValidation rate: %s/%s (%s%%)" % (validated,single_genot+validated+n_validated,round(validation_rate)))


##########################################################################################################################################################################################
##########################################################################################################################################################################################


if __name__ == '__main__':
  main()
