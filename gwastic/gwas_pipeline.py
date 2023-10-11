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

    def validate_gwas_input_files(self, bed_file, pheno_file):
        """Validate the input file.
        Phenotypic file must be ID_space_ID_value.
        fam file ids must match phenotype ids."""

        fam_data = open(bed_file.replace('.bed', '.fam'), 'r').readlines()
        fam_ids = []
        for line in fam_data:
            line = line.strip()
            line = line.split(' ')
            fam_ids.append(line[0])

        pheno_data = open(pheno_file, 'r').readlines()
        pheno_ids = []
        columns = None
        for line2 in pheno_data:
            line2 = line2.strip()
            line2 = line2.split(' ')
            columns = len(line2)
            if len(line2) > 1:
                pheno_ids.append(line2[0])
                #return (True, "Input files validated.")
            else:
                return (False, "Phenotpic file is not space delimited.")

        # Check whether all pheno ids are in the fam file
        check_ids = all(x in pheno_ids for x in fam_ids)

        if check_ids and columns == 3:
            return (True, "Input files validated.")
        elif columns != 3:
            return (False, "Invalid phenotypic file. File should contain 3 coloumns (FID IID Value)")
        else:
            return (False, "Phenotpic ID's does not match with .fam file IDs.")

    def start_gwas(self, bed_file, pheno_file ,chrom_mapping, algorithm, add_log):
        from fastlmm.association import single_snp, single_snp_linreg
        from pysnptools.snpreader import Bed, Pheno
        import pysnptools.util as pstutil
        import fastlmm.util.util as flutil
        import time
        # First we have to validate the input files
        check_input_data = self.validate_gwas_input_files(bed_file, pheno_file)

        if check_input_data[0]:
            bed = Bed(str(bed_file), count_A1=False, chrom_map=chrom_mapping)
            pheno = Pheno(str(pheno_file))
            s1 = 'Raw BED: Nr. of IDs: ' + str(bed.iid_count) + ' Nr.of SNPs: ' + str(bed.sid_count) + ' Nr. of Phenoytpic IDs: ' + str(pheno.iid_count)
            add_log(s1, warn=True)
            # replace original bed with one that has the same iids as the pheno
            bed, pheno = pstutil.intersect_apply([bed, pheno])
            s2 = "After intersection:" + 'bed ids: ' + str(bed.iid_count) + ' SNPs: ' + str(bed.sid_count) + ' Pheno IDs: ' + str(pheno.iid_count)
            add_log(s2, warn=True)
            bed_fixed = self.filter_out_missing(bed)

            ### save data for NN
            import numpy as np
            np.save('snp', bed_fixed.read().val)
            np.save('pheno', pheno.read().val)
            ###

            # format numbers with commas and no decimals
            s3 = "After removing missing:" + 'bed ids: ' + str(bed.iid_count) + ' SNPs: ' + str(bed.sid_count) + ' Pheno IDs: ' + str(pheno.iid_count)
            add_log(s3, warn=True)
            # run single_snp with the fixed file
            add_log('Starting GWAS Analysis, this might take a while...')
            #t1 = time.process_time()
            if algorithm == 'FaST-LMM':
                df = single_snp(bed_fixed, pheno, output_file_name="single_snp.csv")
            elif algorithm == 'Linear regression':
                df = single_snp_linreg(test_snps=bed_fixed, pheno=pheno, output_file_name="single_snp.csv")
            exchanged_dict = {v: k for k, v in chrom_mapping.items()}
            df['Chr'] = df['Chr'].replace(exchanged_dict)
            return df
            #single_snp_linreg(test_snps=bed_fixed, pheno=pheno, output_file_name="single_snp.csv")
            #t2 = time.process_time()
            #print(t2-t1)
        else:
            add_log(check_input_data[1], error=True)

    def plot_gwas(self, df, limit):
        """Manhatten and qq-plot."""
        import matplotlib.pyplot as plt
        import geneview as gv

        dataset2 = df  # Take your df from wherever
        dataset = dataset2[['SNP', 'Chr', 'ChrPos', 'PValue']]
        dataset = dataset.head(limit)
        dataset = dataset.sort_values(by=['Chr', 'ChrPos'])

        #print (dataset)
        ax = gv.manhattanplot(data=dataset, chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP")
        plt.savefig("manhatten.png", dpi=100)

        ax = gv.qqplot(data=dataset["PValue"])
        ax.set(ylim=(0, 10), xlim=(0, 10))
        plt.savefig("qq.png", dpi=100)
        #plt.show()






