import dearpygui.dearpygui as dpg
dpg.create_context()
#import dearpygui_ext.themes as dpg_ext
from dearpygui_ext import logger
import os
import webbrowser
from gwastic_desktop.gwas_pipeline import GWAS
from gwastic_desktop.gp_pipeline import GenomicPrediction
from gwastic_desktop.helpers import HELPERS
from gwastic_desktop.plot_pipeline import Plot
from pysnptools.snpreader import Bed, Pheno
import pysnptools.util as pstutil

#poetry build
#poetry publish
# pipline nr 214


def main():
    app = GWASApp()
    app.run()


class GWASApp:
    def __init__(self):
        self.gwas = GWAS()
        self.helper = HELPERS()
        self.genomic_predict_class = GenomicPrediction()
        self.plot_class = Plot()

        with dpg.font_registry():
            script_dir = os.path.dirname(__file__)  # <-- absolute dir the script is in
            rel_path = "test.ttf"
            abs_file_path = os.path.join(script_dir, rel_path)
            self.font = dpg.add_font(abs_file_path, 20*2, tag="ttf-font")

        with dpg.theme() as self.our_theme:
            with dpg.theme_component(dpg.mvAll):
                dpg.add_theme_color(dpg.mvThemeCol_Button, (79, 143, 192, 255))

        self.vcf_file = ''
        self.pheno_file = ''
        self.vcf_app_data = None
        self.variants_app_data = None
        self.results_directory = None
        self.bed_app_data = None
        self.pheno_app_data = None
        self.default_path = 'Z://gwas_test_data/gwastic_paper/small_set'
        self.gwas_result_name = "gwas_results.csv"
        self.gwas_result_name_top = "gwas_results_top10000.csv"
        self.genomic_predict_name = "genomic_prediction_results.csv"
        self.manhatten_plot_name = "manhatten_plot.png"
        self.qq_plot_name = "qq_plot.png"
        self.gp_plot_name = "Bland_Altman_plot.png"
        self.gp_plot_name_scatter = "GP_scatter_plot.png"
        self.pheno_stats_name = 'pheno_statistics.pdf'
        self.geno_stats_name = 'geno_statistics.pdf'

        self.log_win = dpg.add_window(label="Log", pos=(0, 635), width=1000, height=500, horizontal_scrollbar=True)
        self.logz = logger.mvLogger(self.log_win)
        table_window = None
        # Set up GUI components and callbacks
        self.setup_gui()

    def setup_gui(self):
        dpg.create_viewport(title='GWAStic Desktop Software', width=2000, height=1200, resizable=True)

    # Menu bar
    #     with dpg.viewport_menu_bar():
    #         with dpg.menu(label="Help"):
    #             dpg.add_button(label="Documentation", callback=lambda: webbrowser.open("https://snowformatics.gitbook.io/product-docs/"))

        # File dialogs
        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_vcf, file_count=3, tag="file_dialog_vcf",
                             width=700, height=400, default_path=self.default_path):

            dpg.add_file_extension("Source files (*.vcf *.gz){.vcf,.gz}", color=(255, 255, 0, 255))

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_variants, file_count=3,
                             tag="file_dialog_variants", width=700, height=400, default_path=self.default_path):
            dpg.add_file_extension("Source files (*.txt *.csv){.txt,.csv}", color=(255, 255, 0, 255))
            dpg.add_file_extension(".*")

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_pheno, file_count=3,
                             tag="file_dialog_pheno", width=700, height=400, default_path=self.default_path):
            dpg.add_file_extension("Source files (*.txt *.csv){.txt,.csv}", color=(255, 255, 0, 255))
            dpg.add_file_extension(".*")

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_bed, file_count=3,
                             tag="file_dialog_bed", width=700, height=400, default_path=self.default_path):
            dpg.add_file_extension(".bed", color=(255, 150, 150, 255))
            dpg.add_file_extension(".*")

        dpg.add_file_dialog(
            directory_selector=True, show=False, callback=self.callback_save_results, tag="select_directory",
            cancel_callback=self.cancel_callback_directory, width=700, height=400, default_path=self.default_path)

        with dpg.window(label="GWAStic", width=1000, height=600, no_close=True, horizontal_scrollbar=True):
            with dpg.tab_bar(label='tabbar'):
                with dpg.tab(label='GWAS Analysis'):
                    dpg.add_text("\nStart GWAS Analysis", indent=50)
                    dpg.add_spacer(height=20)
                    geno = dpg.add_button(label="Choose a BED file", callback=lambda: dpg.show_item("file_dialog_bed"), indent=50, tag= 'tooltip_bed')
                    pheno = dpg.add_button(label="Choose a phenotype file", callback=lambda: dpg.show_item("file_dialog_pheno"), indent=50, tag= 'tooltip_pheno')
                    dpg.add_spacer(height=20)
                    self.gwas_combo = dpg.add_combo(label="Select algorithm", items=["FaST-LMM", "Linear regression", "Random Forest (AI)", "XGBoost (AI)"], indent=50, width=200, default_value="FaST-LMM", tag= 'tooltip_algorithm')
                    dpg.add_spacer(height=20)
                    gwas_btn = dpg.add_button(label="Run GWAS", callback=self.run_gwas, user_data=[geno, pheno], indent=50)
                    dpg.bind_item_theme(gwas_btn, self.our_theme)

                with dpg.tab(label='Genomic Prediction Analysis'):
                    dpg.add_text("\nStart Genomic Prediction", indent=50)
                    dpg.add_spacer(height=20)
                    geno = dpg.add_button(label="Choose a BED file", callback=lambda: dpg.show_item("file_dialog_bed"),
                                          indent=50)
                    pheno = dpg.add_button(label="Choose a phenotype file",
                                           callback=lambda: dpg.show_item("file_dialog_pheno"), indent=50)
                    dpg.add_spacer(height=20)
                    self.gwas_gp = dpg.add_combo(label="Select Algorithm", items=["XGBoost (AI)", "Random Forest (AI)", 'GP_LMM'],#, 'Test'],
                                  indent=50, width=200, default_value="XGBoost (AI)")
                    dpg.add_spacer(height=20)
                    gwas_btn = dpg.add_button(label="Run Genomic Prediction", callback=self.run_genomic_prediction, user_data=[geno, pheno],
                                              indent=50)
                    dpg.bind_item_theme(gwas_btn, self.our_theme)

                with dpg.tab(label='Convert VCF'):
                    dpg.add_text("\nConvert a VCF file into BED file format and apply MAF\nor missing genotype filter.", indent=50)
                    dpg.add_spacer(height=20)
                    dpg.add_text("Select files:", indent=50)
                    vcf = dpg.add_button(label="Choose a VCF file", callback=lambda: dpg.show_item("file_dialog_vcf"),
                                         indent=50, tag='tooltip_vcf')

                    variant_ids = dpg.add_button(label="Choose a variants file (optional)", tag= 'tooltip_variant',
                                                 callback=lambda: dpg.show_item("file_dialog_variants"), indent=50)
                    dpg.add_spacer(height=20)
                    dpg.add_text("Apply filters:", indent=50)
                    maf_input = dpg.add_input_float(label="Minor allele frequency (MAF)", width=150, default_value=0.05,step=0.005, indent=50, tag= 'tooltip_maf')
                    geno_input = dpg.add_input_float(label="Missing genotype rate", width=150, default_value=0.1, step=0.005, indent=50, tag= 'tooltip_missing')
                    dpg.add_spacer(height=20)
                    convert_btn = dpg.add_button(label="Convert VCF", callback=self.convert_vcf, user_data=[maf_input, geno_input, vcf, variant_ids], indent=50)
                    dpg.bind_item_theme(convert_btn, self.our_theme)

                with dpg.tab(label='Settings'):
                    dpg.add_spacer(height=10)
                    dpg.add_text("Machine Learning Settings", indent=50, color=(72,138,199))
                    dpg.add_spacer(height=10)
                    #self.std_set = dpg.add_checkbox(label="Apply Standardization", indent=50, tag= 'tooltip_stand')
                    #dpg.add_spacer(height=10)
                    self.train_size_set = dpg.add_input_int(label="Training size", width=150, default_value=70,step=10, indent=50,
                                      min_value=0, max_value=100, min_clamped=True, max_clamped=True, tag= 'tooltip_training')
                    dpg.add_spacer(height=10)
                    self.model_nr = dpg.add_input_int(label="Nr. of models", width=150, default_value=1, step=1,
                                                            indent=50,
                                                            min_value=1, max_value=50, min_clamped=True,
                                                            max_clamped=True, tag='tooltip_model')
                    dpg.add_spacer(height=10)
                    self.estim_set = dpg.add_input_int(label="Number of trees", width=150, default_value=200, step=10, indent=50,
                                      min_value=1, min_clamped=True, tag= 'tooltip_trees')
                    dpg.add_spacer(height=10)
                    self.max_dep_set = dpg.add_input_int(label="Max depth", width=150, default_value=3, step=10, indent=50,
                                      min_value=0, max_value=100, min_clamped=True, max_clamped=True, tag= 'tooltip_depth')
                    dpg.add_spacer(height=10)

            # Tooltips
            with dpg.tooltip("tooltip_vcf"):
                dpg.add_text("Click to upload a Variant Call Format (VCF) file (.gz or .vcf).", color=[79,128,226])
            with dpg.tooltip("tooltip_variant"):
                dpg.add_text("Choose a variant file.\nVariant IDs must match with IDs in the VCF file.\nBy using a variant file, you can create a subset of your VCF file.\n Recommended if only a subset with phenotypic data is available.", color=[79,128,226])
            with dpg.tooltip("tooltip_maf"):
                dpg.add_text("Set the Minor Allele Frequency (MAF) threshold for filtering variants.\nVariants with a MAF below this value will be excluded.", color=[79, 128, 226])
            with dpg.tooltip("tooltip_missing"):
                dpg.add_text("Enter the maximum allowable rate of missing genotypes per variant. \nVariants with missing data above this rate will be excluded.", color=[79, 128, 226])
            with dpg.tooltip("tooltip_bed"):
                dpg.add_text("Click to select a Binary BED file containing genotype data.\nMust be accompanied by .bim and .fam files.\n", color=[79, 128, 226])
            with dpg.tooltip("tooltip_pheno"):
                dpg.add_text("Click to select a file with phenotype data that will be used in the GWAS analysis.\nVariant IDs must match with IDs in the .fam file.\nMust be space seperated.\nExample.\nID1 ID1 0.25\nID2 ID2 0.89\nImportant:ID's must not contain spaces", color=[79, 128, 226])
            with dpg.tooltip("tooltip_algorithm"):
                dpg.add_text("Select the algorithm to be used for the analysis.", color=[79, 128, 226])
           # with dpg.tooltip("tooltip_pvalue"):
                #dpg.add_text("All output rows with p-values less than this threshold will be included.\nBy default, all rows are included.\nThis is used to exclude rows for large datasets.", color=[79, 128, 226])
            #with dpg.tooltip("tooltip_chrom"):
                #dpg.add_text("Perform single SNP GWAS via cross validation over the chromosomes.\nDefault to True.\nWarning: setting False can cause proximal contamination.", color=[79, 128, 226])
            #with dpg.tooltip("tooltip_stand"):
                #dpg.add_text("Check this to standardize features by removing the mean and scaling to unit variance,\noften required for machine learning algorithms.", color=[79, 128, 226])
            with dpg.tooltip("tooltip_training"):
                dpg.add_text("Set the percentage of the dataset to be used for training the model.\nThe rest will be used for testing.", color=[79, 128, 226])
            with dpg.tooltip("tooltip_trees"):
                dpg.add_text("Specify the number of trees to be used in the forest.\nMore trees can increase accuracy but also computation time.", color=[79, 128, 226])
            with dpg.tooltip("tooltip_model"):
                dpg.add_text("Specify the number of models to be used in the analysis.\nMore models can increase accuracy but also computation time.", color=[79, 128, 226])
            with dpg.tooltip("tooltip_depth"):
                dpg.add_text("Determine the maximum depth of the trees.\nDeeper trees can model more complex relationships.", color=[79, 128, 226])

            dpg.bind_font(self.font)
            dpg.set_global_font_scale(0.6)

    def callback_vcf(self, sender, app_data):
        """Get vcf file path selected from the user."""
        self.vcf_app_data = app_data
        vcf_path, current_path = self.get_selection_path(self.vcf_app_data)
        self.add_log('VCF Selected: ' + vcf_path)

    def callback_bed(self, sender, app_data):
        """Get vcf file path selected from the user."""
        self.bed_app_data = app_data
        try:
            bed_path, current_path = self.get_selection_path(self.bed_app_data)
            self.add_log('BED file Selected: ' + bed_path)
        except TypeError:
            self.add_log('Invalid BED file Selected', error=True)

    def callback_variants(self, sender, app_data):
        """Get variant file path selected from the user."""
        self.variants_app_data = app_data
        variants_path, current_path = self.get_selection_path(self.variants_app_data)
        self.add_log('Variant File Selected: ' + variants_path)

    def callback_pheno(self, sender, app_data):
        """Get phenotype file path selected from the user."""
        self.pheno_app_data = app_data
        try:
            pheno_path, current_path = self.get_selection_path(self.pheno_app_data)
            self.add_log('Pheno File Selected: ' + pheno_path)
        except TypeError:
            self.add_log('Wrong Pheno File Selected: ' + pheno_path, error=True)

    def get_selection_path(self, app_data):
        """Extract path from the app_data dictionary selections key."""

        current_path = app_data['current_path'] + '/'
        k = app_data['selections']
        try:
            for key, value in k.items():
                file_path = value
            return file_path, current_path
        except UnboundLocalError:
            pass

    def callback_save_results(self, sender, app_data):
        """Save the results inside the folder including tables as csv and plots as png. """
        self.results_directory = app_data
        results_path, current_path = self.get_selection_path(self.results_directory)
        save_dir = self.helper.save_results(os.getcwd(), current_path, self.gwas_result_name, self.gwas_result_name_top,
                                            self.manhatten_plot_name, self.qq_plot_name, self.algorithm,
                                            self.genomic_predict_name, self.gp_plot_name, self.gp_plot_name_scatter,
                                            self.add_log, self.settings_lst, self.pheno_stats_name, self.geno_stats_name )
        self.add_log('Results saved in: ' + save_dir)

    def cancel_callback_directory(self, sender, app_data):
        self.add_log('Process Canceled')

    def save_default_path(self, sender, data, user_data):
        """Overwrite the default path. Restart necessary."""
        default_path = str(dpg.get_value(user_data[0]))
        self.helper.save_settings(default_path)
        self.add_log('Settings saved. Please restart the software.', warn=True)

    def delete_files(self):
        """Delete the temporary files after analysis."""
        dpg.delete_item("manhatten_image")
        dpg.delete_item("manhatten_tag")
        dpg.delete_item("qq_image")
        dpg.delete_item("qq_tag")
        dpg.delete_item("table_gwas")
        dpg.delete_item("table_gp")
        dpg.delete_item("ba_tag")
        dpg.delete_item("ba_tag2")
        dpg.delete_item("ba_image")
        dpg.delete_item("ba_image2")

        # delete files
        file_names = [self.gwas_result_name, self.gwas_result_name_top, self.genomic_predict_name, self.gp_plot_name,
                      self.manhatten_plot_name, self.qq_plot_name, self.gp_plot_name_scatter,
                      self.manhatten_plot_name.replace('manhatten_plot', 'manhatten_plot_high'),
                     self.qq_plot_name.replace('qq_plot', 'qq_plot_high'),
                      self.gp_plot_name_scatter.replace('GP_scatter_plot', 'GP_scatter_plot_high'),
                      self.gp_plot_name.replace('Bland_Altman_plot', 'Bland_Altman_plot_high'),
                      self.genomic_predict_name.replace('.csv', '_valdation.csv'),
                      self.pheno_stats_name, self.geno_stats_name
                      ]
        for f in file_names:
            if os.path.exists(f):
                os.remove(f)

    def add_log(self, message, warn=False, error=False):
        """Adds a log message."""
        if warn:
            self.logz.log_warning(message)
        elif error:
            self.logz.log_error(message)
        else:
            self.logz.log_info(message)

    def convert_vcf(self, sender, data, user_data):
        """Converts vcf to bed file using PLINK."""
        maf = str(dpg.get_value(user_data[0]))
        geno = str(dpg.get_value(user_data[1]))
        variant_ids = self.variants_app_data
        vcf_path, current_path = self.get_selection_path(self.vcf_app_data)

        if variant_ids is None:
            variants_path = None
        else:
            variants_path, current_path = self.get_selection_path(self.variants_app_data)
        self.add_log('Start converting VCF to BED...')
        plink_log = self.gwas.vcf_to_bed(vcf_path, variants_path,
                                         vcf_path.split('.')[0] +'_maf' + str(round(float(maf),2)) +'_geno' + str(round(float(geno),2)),
                                         maf, geno)
        self.add_log(plink_log)

    def run_gwas(self, sender, data, user_data):
        """Starts the GWAS pipeline."""

        self.delete_files()

        # Get all settings
        train_size_set = (100-dpg.get_value(self.train_size_set))/100
        estimators = dpg.get_value(self.estim_set)
        model_nr = dpg.get_value(self.model_nr)
        #pvalue_set = dpg.get_value(self.pvalue_set)
        #std_set = dpg.get_value(self.std_set)
        max_dep_set = dpg.get_value(self.max_dep_set)
        self.algorithm = dpg.get_value(self.gwas_combo)

        try:
            self.add_log('Reading files...')
            bed_path, current_path1 = self.get_selection_path(self.bed_app_data)
            pheno_path, current_path2 = self.get_selection_path(self.pheno_app_data)
            self.add_log('Validating files...')
            check_input_data = self.gwas.validate_gwas_input_files(bed_path, pheno_path)
            # Replace chromosome names, they need to be numbers
            chrom_mapping = self.helper.replace_with_integers(bed_path.replace('.bed', '.bim'))
            self.settings_lst = [self.algorithm, bed_path, pheno_path, train_size_set, estimators, model_nr, max_dep_set]
            if check_input_data[0]:
                bed = Bed(str(bed_path), count_A1=False, chrom_map=chrom_mapping)
                pheno = Pheno(str(pheno_path))
                # replace original bed with one that has the same iids as the pheno
                bed, pheno = pstutil.intersect_apply([bed, pheno])
                bed_fixed = self.gwas.filter_out_missing(bed)
                # format numbers with commas and no decimals
                s3 = "Dataset after intersection:" + ' SNPs: ' + str(bed.sid_count) + ' Pheno IDs: ' + str(
                    pheno.iid_count)
                self.add_log(s3, warn=True)
                # run single_snp with the fixed file
                self.add_log('Starting Analysis, this might take a while...')
                if self.algorithm == 'FaST-LMM' or self.algorithm == 'Linear regression':
                    gwas_df, df_plot =self.gwas.run_gwas_lmm(bed_fixed, pheno, chrom_mapping, self.add_log
                                                             , self.gwas_result_name, self.algorithm, bed_path)
                elif  self.algorithm == 'Random Forest (AI)':
                    gwas_df, df_plot = self.gwas.run_gwas_rf(bed_fixed, pheno, bed_path, train_size_set,
                                                            estimators, self.gwas_result_name, chrom_mapping,
                                                            self.add_log, model_nr)
                elif self.algorithm == 'XGBoost (AI)':
                    gwas_df, df_plot = self.gwas.run_gwas_xg(bed_fixed, pheno, bed_path, train_size_set, estimators,
                                                             self.gwas_result_name, chrom_mapping, self.add_log,
                                                             model_nr, max_dep_set)
            else:
                self.add_log(check_input_data[1], error=True)

            if gwas_df is not None:
                self.add_log('GWAS Analysis done.')
                self.add_log('GWAS Results Plotting...')
                self.plot_class.plot_pheno_statistics(pheno_path, self.pheno_stats_name)
                self.plot_class.plot_geno_statistics(bed_fixed, pheno, self.geno_stats_name)
                self.gwas.plot_gwas(df_plot, 10000, self.algorithm, self.manhatten_plot_name, self.qq_plot_name, chrom_mapping)

                self.add_log('Done...')
                self.show_results_window(gwas_df, self.algorithm, genomic_predict=False)

            else:
                self.add_log('Error, GWAS Analysis could not be started.', error=True)

        except TypeError:
            self.add_log('Please select a phenotype and genotype file. ', error=True)

    def run_genomic_prediction(self, sender, data, user_data):
        """Starts the GP pipeline."""

        self.delete_files()
        self.add_log('Reading Bed file...')

        # Get settings
        self.algorithm = dpg.get_value(self.gwas_gp)
        test_size = (100 - dpg.get_value(self.train_size_set)) / 100
        estimators = dpg.get_value(self.estim_set)
        max_dep_set = dpg.get_value(self.max_dep_set)
        model_nr = dpg.get_value(self.model_nr)
        try:
            self.add_log('Reading files...')
            bed_path, current_path1 = self.get_selection_path(self.bed_app_data)
            pheno_path, current_path2 = self.get_selection_path(self.pheno_app_data)
            self.add_log('Validating files...')
            check_input_data = self.gwas.validate_gwas_input_files(bed_path, pheno_path)
            # Replace chromosome names, they need to be numbers
            chrom_mapping = self.helper.replace_with_integers(bed_path.replace('.bed', '.bim'))
            self.settings_lst = [self.algorithm, bed_path, pheno_path, test_size, estimators, model_nr, max_dep_set]
            if check_input_data[0]:
                bed = Bed(str(bed_path), count_A1=False, chrom_map=chrom_mapping)
                pheno = Pheno(str(pheno_path))
                # replace original bed with one that has the same iids as the pheno
                bed, pheno = pstutil.intersect_apply([bed, pheno])
                bed_fixed = self.gwas.filter_out_missing(bed)
                # format numbers with commas and no decimals
                s3 = "Dataset after intersection:" + ' SNPs: ' + str(bed.sid_count) + ' Pheno IDs: ' + str(
                    pheno.iid_count)
                self.add_log(s3, warn=True)
                self.add_log('Starting Analysis, this might take a while...')
                if self.algorithm == 'GP_LMM':
                    gp_df = self.genomic_predict_class.run_lmm_gp(bed_fixed, pheno, self.genomic_predict_name, model_nr, self.add_log,
                                                 bed_path, chrom_mapping)
                elif self.algorithm == 'Random Forest (AI)':
                    gp_df = self.genomic_predict_class.run_gp_rf(bed_fixed, pheno, bed_path, test_size, estimators,
                                                self.genomic_predict_name, chrom_mapping, self.add_log,model_nr)
                elif self.algorithm == 'XGBoost (AI)':
                    gp_df, df_val = self.genomic_predict_class.run_gp_xg(bed_fixed, pheno, bed_path, test_size, estimators,
                                                self.genomic_predict_name, chrom_mapping, self.add_log,model_nr, max_dep_set)
                else:
                    self.genomic_predict_class.model_validation(bed_fixed, pheno, bed_path, test_size, estimators,
                                                                self.genomic_predict_name,  chrom_mapping, self.add_log,
                                                                 model_nr, max_dep_set, validation_size=0.1)
            else:
                self.add_log(check_input_data[1], error=True)

            if gp_df is not None:
                self.add_log('Genomic Prediction done.')
                self.add_log('Genomic Prediction Plotting...')
                self.genomic_predict_class.plot_gp(gp_df, self.gp_plot_name, self.algorithm)
                self.genomic_predict_class.plot_gp_scatter(gp_df, self.gp_plot_name_scatter, self.algorithm)
                try:
                    if df_val is not None:
                        self.genomic_predict_class.plot_gp_scatter(df_val, 'validation_' + self.gp_plot_name_scatter, self.algorithm)
                except UnboundLocalError:
                    pass
                self.add_log('Done...')
                self.show_results_window(gp_df, self.algorithm, genomic_predict=True)
            else:
                self.add_log('Error, GWAS Analysis could not be started.', error=True)

        except TypeError:
            self.add_log('Please select a phenotype and genotype file. ', error=True)

    def show_results_window(self, df, algorithm, genomic_predict):
        """Create a new window to show the plots and tables."""
        with dpg.window(label="Results", width=975, height=600, horizontal_scrollbar=True, pos=(1000, 35)):
            dpg.add_button(label="Export Results", pos =(400, 40), callback=lambda: dpg.show_item("select_directory"))
            dpg.add_spacer(height=60)

            if genomic_predict:
                width, height, channels, data = dpg.load_image(self.gp_plot_name)
                with dpg.texture_registry(show=False):
                    dpg.add_static_texture(width=width, height=height, default_value=data, tag="ba_tag")

                width, height, channels, data = dpg.load_image(self.gp_plot_name_scatter)
                with dpg.texture_registry(show=False):
                    dpg.add_static_texture(width=width, height=height, default_value=data, tag="ba_tag2")

                with dpg.tab_bar(label='tabbar'):
                    with dpg.tab(label="Genomic Prediction Results"):
                        df = df[['ID1', 'BED_ID2_x', 'Mean_Predicted_Value', 'Pheno_Value', 'Difference']]
                        df.columns = ['FID', 'IID', 'Predicted_Value', 'Pheno_Value', 'Difference']
                        # df = df.sort_values(by=['PValue'], ascending=True)
                        with dpg.table(label='DatasetTable2', row_background=True, borders_innerH=True,
                                       borders_outerH=True, borders_innerV=True, borders_outerV=True, tag='table_gp',
                                       sortable=True, resizable=True):
                            for i in range(df.shape[1]):
                                dpg.add_table_column(label=df.columns[i], parent='table_gp')
                            for i in range(len(df)):
                                with dpg.table_row():
                                    for j in range(df.shape[1]):
                                        #dpg.add_text(f"{df.iloc[i, j]}")
                                        value = df.iloc[i, j]
                                        # Format the value as a string with desired precision, e.g., 2 decimal places
                                        formatted_value = f"{value:.2f}" if isinstance(value, float) else str(value)
                                        dpg.add_text(formatted_value)
                    with dpg.tab(label="Correlation Plot (Predicted vs. Phenotype)"):
                        dpg.add_image(texture_tag="ba_tag2", tag="ba_image2", width=750, height=450)
                    with dpg.tab(label="Bland-Altman Plot (Model Accuracy)"):
                        dpg.add_image(texture_tag="ba_tag", tag="ba_image", width=750, height=450)

            else:
                width, height, channels, data = dpg.load_image(self.manhatten_plot_name)
                with dpg.texture_registry(show=False):
                    dpg.add_static_texture(width=width, height=height, default_value=data, tag="manhatten_tag")
                with dpg.tab_bar(label='tabbar'):
                    with dpg.tab(label="Manhatten Plot"):
                        if algorithm == "FaST-LMM" or algorithm == "Linear regression":
                            dpg.add_image(texture_tag="manhatten_tag", tag="manhatten_image", width=950, height=400)
                        else:
                            dpg.add_image(texture_tag="manhatten_tag", tag="manhatten_image", width=900, height=300)
                    if algorithm == "FaST-LMM" or algorithm == "Linear regression":
                        width2, height2, channels2, data2 = dpg.load_image(self.qq_plot_name)
                        with dpg.texture_registry(show=False):
                            dpg.add_static_texture(width=width2, height=height2, default_value=data2, tag="qq_tag")
                        df = df.sort_values(by=['PValue'], ascending=True)
                        with dpg.tab(label="QQ-Plot"):
                            dpg.add_image(texture_tag="qq_tag", tag="qq_image", height=450, width=450)
                    else:

                        df = df.sort_values(by=['SNP effect'], ascending=False)

                    with dpg.tab(label="GWAS Results (Top 500)"):
                        if algorithm == "FaST-LMM" or algorithm == "Linear regression":
                            df = df[['SNP', 'Chr', 'ChrPos', 'PValue']]
                        else:
                            df.columns = df.columns.str.replace('SNP effect_sd', 'SNP effect SD')

                        max_rows= len(df)
                        if max_rows > 501:
                            max_rows = 500
                        with dpg.table(label='DatasetTable',row_background=True, borders_innerH=True,
                                       borders_outerH=True, borders_innerV=True, borders_outerV=True, tag='table_gwas',
                                       sortable=True):
                            for i in range(df.shape[1]):
                                dpg.add_table_column(label=df.columns[i], parent='table_gwas')
                            for i in range(max_rows):
                                with dpg.table_row():
                                    for j in range(df.shape[1]):
                                        value = df.iloc[i, j]
                                        # Format the value as a string with desired precision, e.g., 2 decimal places
                                        #formatted_value = f"{value:.6f}" if isinstance(value, float) else str(value)
                                        dpg.add_text(value)

    def run(self):
        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.start_dearpygui()
        dpg.destroy_context()


if __name__ == "__main__":
    main()

