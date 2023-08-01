import subprocess

class GWAS:
    """GWAS class."""

    def vcf_to_bed(self, vcf_file, id_file, file_out, maf, geno):
        """Converts the vcf to bed files."""
        if id_file == None:
            process = subprocess.Popen(["plink", "--vcf", vcf_file, "--make-bed", "--out", file_out,
                                        "--allow-extra-chr", "--set-missing-var-ids", "@:#", "--maf", maf,
                                        "--geno", geno, "--double-id"])

        else:
            process = subprocess.Popen(["plink", "--vcf", vcf_file, "--make-bed", "--out", file_out,
                                        "--allow-extra-chr", "--set-missing-var-ids", "@:#", "--maf", maf,
                                        "--geno", geno, "--double-id", "--keep", id_file])

        process.wait()
        # read the log file
        plink_log = open(file_out + '.log').read()
        return plink_log

    def filter_out_missing(self, bed):
        import numpy as np
        import logging
        from pysnptools.util import log_in_place
        sid_batch_size = 1000  # number of SNPs to read at a time, decrease in case we run oput of memory
        # list of all NaN columns
        all_nan = []
        with log_in_place("read snp #", logging.INFO) as updater:
            for sid_start in range(0, bed.sid_count, sid_batch_size):
                updater(f"sid {sid_start:,} of {bed.sid_count:,}")
                # read a batch of SNPs
                snp_data_batch = bed[:, sid_start:sid_start + sid_batch_size].read()
                # find the columns that are all NaN (missing)
                nan_in_batch = np.isnan(snp_data_batch.val).all(axis=0)
                all_nan.append(nan_in_batch)

        # concatenate the list of all NaN columns into a single boolean array
        all_nan = np.concatenate(all_nan, axis=0)
        # print the number of all NaN columns
        logging.info(f"number of all missing columns is {np.sum(all_nan):,}. They are: ")
        logging.info(bed.sid[all_nan])
        # remove these SNPs from consideration
        bed_fixed = bed[:, ~all_nan]
        return bed_fixed

    def start_gwas(self, bed_file, pheno_file):
        from fastlmm.association import single_snp
        from pysnptools.snpreader import Bed, Pheno
        import pysnptools.util as pstutil
        bed = Bed(str(bed_file), count_A1=False)
        pheno = Pheno(str(pheno_file))
        print(f"Before, bed is {bed.iid_count:,} x {bed.sid_count:,}, pheno is {pheno.iid_count:,}")
        # replace original bed with one that has the same iids as the pheno
        bed, pheno = pstutil.intersect_apply([bed, pheno])

        print(f"After intersection, bed is {bed.iid_count} x {bed.sid_count:,}, pheno is {pheno.iid_count:,}")
        bed_fixed = self.filter_out_missing(bed)
        # format numbers with commas and no decimals
        print(f"After removing missing, bed is {bed.iid_count:,} x {bed.sid_count:,}, pheno is {pheno.iid_count:,}")
        # run single_snp with the fixed file
        single_snp(bed_fixed, pheno, output_file_name="single_snp.csv")

    def plot_gwas(self, gwas_file, limit):
        """Manhatten and qq-plot."""
        import matplotlib.pyplot as plt
        import geneview as gv
        import pandas as pd
        dataset2 = pd.read_csv(gwas_file, delimiter='\t')  # Take your df from wherever
        dataset = dataset2[['SNP', 'Chr', 'ChrPos', 'PValue']]
        dataset = dataset.head(limit)
        dataset = dataset.sort_values(by=['Chr', 'ChrPos'])

        #print (dataset)
        ax = gv.manhattanplot(data=dataset, chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP")
        plt.savefig("manhatten.png", dpi=100)
        ax = gv.qqplot(data=dataset["PValue"])
        plt.savefig("qq.png", dpi=100)
        #plt.show()






