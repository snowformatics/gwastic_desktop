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

