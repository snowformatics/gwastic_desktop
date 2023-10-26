import dearpygui.dearpygui as dpg
import dearpygui_ext.themes as dpg_ext
from dearpygui_ext import logger
from gwastic_desktop.gwas_pipeline import GWAS
from gwastic_desktop.helpers import HELPERS
import os

def main():
    app = GWASApp()
    app.run()


class GWASApp:
    def __init__(self):
        self.gwas = GWAS()
        self.helper = HELPERS()
        dpg.create_context()

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
        self.algorithm = "FaST-LMM"
        self.default_path = "C:/gwas_test_data/test/vcf2gwas/"
        self.gwas_result_name = "gwas_results.csv"
        self.gwas_result_name_top = "gwas_results_top10000.csv"
        self.genomic_predict_name = "genomic_prediction_results.csv"
        self.manhatten_plot_name = "manhatten_plot.png"
        self.qq_plot_name = "qq_plot.png"
        self.test_size = 0.2
        self.estimators = 100

        self.log_win = dpg.add_window(label="Log", pos=(0, 635), width=1000, height=500)
        self.logz = logger.mvLogger(self.log_win)
        table_window = None
        # Set up GUI components and callbacks
        self.setup_gui()

    def save_callback(self):
        print("Save Clicked")

    def setup_gui(self):
        dpg.create_viewport(title='GWAStic Desktop Software', width=2000, height=1200)

        def print_me(sender):
            print(f"Menu Item: {sender}")


    # Menu bar
        with dpg.viewport_menu_bar():
            # with dpg.menu(label="Settings"):
            #     with dpg.menu(label="Methods"):
            #         dpg.add_menu_item(label="Decision Tree Options", callback=print_me, check=True)
            #         dpg.add_input_float(label="Test Size", default_value=0.2, min_value=0.0, max_value=1.0,  width=150)
            #         dpg.add_input_int(label="Number of trees", default_value=100, width=150)
            #         dpg.add_input_int(label="Depth of the tree", width=150)
            #         dpg.add_spacer(height=20)
            #         dpg.add_menu_item(label="Other", callback=print_me, check=True)

            with dpg.menu(label="Help"):
                dpg.add_menu_item(label="Documentation", callback=print_me, check=True)
                dpg.add_menu_item(label="Tutorials", callback=print_me, check=True)

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
            cancel_callback=self.cancel_callback_directory, width=700, height=400)

        with dpg.window(label="GWAStic", width=1000, height=600, no_close=True):
            with dpg.tab_bar(label='tabbar'):
                with dpg.tab(label='GWAS Analysis'):
                    dpg.add_text("\nStart GWAS Analysis", indent=50)
                    dpg.add_spacer(height=20)
                    geno = dpg.add_button(label="Choose BED", callback=lambda: dpg.show_item("file_dialog_bed"), indent=50, tag= 'tooltip_bed')
                    pheno = dpg.add_button(label="Choose Phenotype", callback=lambda: dpg.show_item("file_dialog_pheno"), indent=50, tag= 'tooltip_pheno')
                    dpg.add_spacer(height=20)
                    dpg.add_combo(label="Algorithm", items=["FaST-LMM", "Linear regression", "Random Forest (AI)", "XGBoost (AI)"], indent=50, width=200, default_value="FaST-LMM", callback=self.get_algorithm)
                    dpg.add_spacer(height=20)
                    gwas_btn = dpg.add_button(label="Run GWAS", callback=self.run_gwas, user_data=[geno, pheno], indent=50)
                    dpg.bind_item_theme(gwas_btn, self.our_theme)

                with dpg.tab(label='Genomic Prediction Analysis'):
                    dpg.add_text("\nStart Genomic Prediction", indent=50)
                    dpg.add_spacer(height=20)
                    geno = dpg.add_button(label="Choose BED", callback=lambda: dpg.show_item("file_dialog_bed"),
                                          indent=50)
                    pheno = dpg.add_button(label="Choose Phenotype",
                                           callback=lambda: dpg.show_item("file_dialog_pheno"), indent=50)
                    dpg.add_spacer(height=20)
                    self.gp_combo = dpg.add_combo(label="Algorithm",
                                  items=["Random Forest (AI)", "XGBoost (AI)"],
                                  indent=50, width=200, default_value="Random Forest (AI)", callback=self.get_algorithm)
                    dpg.add_spacer(height=20)
                    gwas_btn = dpg.add_button(label="Run Genomic Prediction", callback=self.run_genomic_prediction, user_data=[geno, pheno],
                                              indent=50)
                    dpg.bind_item_theme(gwas_btn, self.our_theme)

                with dpg.tab(label='Convert VCF'):
                    dpg.add_text("\nConvert a VCF file into BED file format and apply MAF\nor missing genotype filter.", indent=50)
                    dpg.add_spacer(height=20)
                    dpg.add_text("Select files:", indent=50)
                    vcf = dpg.add_button(label="Choose VCF", callback=lambda: dpg.show_item("file_dialog_vcf"),
                                         indent=50, tag='tooltip_vcf')

                    variant_ids = dpg.add_button(label="Choose Variants (optional)", tag= 'tooltip_variant',
                                                 callback=lambda: dpg.show_item("file_dialog_variants"), indent=50)
                    dpg.add_spacer(height=20)
                    dpg.add_text("Apply filters:", indent=50)
                    maf_input = dpg.add_input_float(label="MAF", width=150, default_value=0.05,step=0.005, indent=50, tag= 'tooltip_maf')
                    geno_input = dpg.add_input_float(label="Missing genotype rate", width=150, default_value=0.1, step=0.005, indent=50, tag= 'tooltip_missing')
                    dpg.add_spacer(height=20)
                    convert_btn = dpg.add_button(label="Convert", callback=self.convert_vcf, user_data=[maf_input, geno_input, vcf, variant_ids], indent=50)
                    dpg.bind_item_theme(convert_btn, self.our_theme)

            # Tooltips
            with dpg.tooltip("tooltip_vcf"):
                dpg.add_text("Choose a VCF file (.gz or .vcf).", color=[79,128,226])
            with dpg.tooltip("tooltip_variant"):
                dpg.add_text("Choose a variant file.\nVariant IDs must match with IDs in the VCF file.\n By using a variant file, you can create a subset of your VCF file.\n Highly recommended if only a subset with phenotypic data is available).", color=[79,128,226])
            with dpg.tooltip("tooltip_maf"):
                dpg.add_text("Filters out all variants with minor allele frequency below the provided threshold.", color=[79, 128, 226])
            with dpg.tooltip("tooltip_missing"):
                dpg.add_text("Filters out all variants with missing call rates exceeding the provided value to be removed.", color=[79, 128, 226])
            with dpg.tooltip("tooltip_bed"):
                dpg.add_text("Choose a BED file.\nPrimary representation of genotype calls at biallelic variants. Must be accompanied by .bim and .fam files.\n", color=[79, 128, 226])
            with dpg.tooltip("tooltip_pheno"):
                dpg.add_text("Choose a phenotype file.\nVariant IDs must match with IDs in the .fam file.\nMust be space seperated.\nExample.\nID1 ID1 0.25\nID2 ID2 0.89\nImportant:ID's must not contain spaces", color=[79, 128, 226])

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
        bed_path, current_path = self.get_selection_path(self.bed_app_data)
        self.add_log('BED file Selected: ' + bed_path)

    def callback_variants(self, sender, app_data):
        """Get variant file path selected from the user."""
        self.variants_app_data = app_data
        variants_path, current_path = self.get_selection_path(self.variants_app_data)
        self.add_log('Variant File Selected: ' + variants_path)

    def callback_pheno(self, sender, app_data):
        """Get phenotype file path selected from the user."""
        self.pheno_app_data = app_data
        pheno_path, current_path = self.get_selection_path(self.pheno_app_data)
        self.add_log('Pheno File Selected: ' + pheno_path)

    def get_selection_path(self, app_data):
        """Extract path from the app_data dictionary selections key."""

        current_path = app_data['current_path'] + '/'
        k = app_data['selections']
        for key, value in k.items():
            file_path = value
        return file_path, current_path

    def callback_save_results(self, sender, app_data):
        """"""
        self.results_directory = app_data
        results_path, current_path = self.get_selection_path(self.results_directory)
        self.helper.save_results(os.getcwd(), current_path, self.gwas_result_name, self.gwas_result_name_top,
                                 self.manhatten_plot_name, self.qq_plot_name, self.algorithm, self.genomic_predict_name)
        self.add_log('Results saved in: ' + current_path)

    def cancel_callback_directory(self, sender, app_data):
        self.add_log('Process Canceled')

    def get_algorithm(self, sender, data):

        self.algorithm = data

    def delete_files(self, genomic_predict):

        dpg.delete_item("manhatten_image")
        dpg.delete_item("manhatten_tag")
        dpg.delete_item("qq_image")
        dpg.delete_item("qq_tag")
        dpg.delete_item("table_gwas")
        dpg.delete_item("table_gp")

        # delete files
        file_names = [self.gwas_result_name, self.gwas_result_name_top, self.genomic_predict_name,
                      self.manhatten_plot_name, self.qq_plot_name]
        for f in file_names:
            if os.path.exists(f):
                print (f)
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
        plink_log = self.gwas.vcf_to_bed(vcf_path, variants_path, vcf_path.split('.')[0], maf, geno)

        self.add_log(plink_log)
        #self.logz.log_debug('Converting done!')

    def run_gwas(self, sender, data, user_data):
        self.delete_files(genomic_predict = False)

        self.add_log('Reading Bed file...')
        try:
            bed_path, current_path1 = self.get_selection_path(self.bed_app_data)
            self.add_log('Reading Phenotypic file...')
            pheno_path, current_path2 = self.get_selection_path(self.pheno_app_data)
            # Replace chromosome names, they need to be numbers
            chrom_mapping = self.helper.replace_with_integers(bed_path.replace('.bed', '.bim'))
            gwas_df = self.gwas.start_gwas(bed_path, pheno_path, chrom_mapping, self.algorithm, self.add_log,
                                           self.test_size, self.estimators, self.gwas_result_name, False, None)
            if gwas_df is not None:
                self.add_log('GWAS Analysis done.')
                self.add_log('GWAS Results Plotting...')
                self.gwas.plot_gwas(gwas_df, 10000, self.algorithm, self.manhatten_plot_name, self.qq_plot_name)
                self.add_log('Done...')
                self.show_results_window(gwas_df, self.algorithm, genomic_predict=False)
            else:
                self.add_log('Error, GWAS Analysis could not be started.', error=True)

        except TypeError:
            self.add_log('Please select a phenotype and genotype file. ', error=True)

    def run_genomic_prediction(self, sender, data, user_data):
        self.delete_files(genomic_predict = True)
        self.add_log('Reading Bed file...')
        self.algorithm = dpg.get_value(self.gp_combo)
        #print (self.algorithm)

        try:
            bed_path, current_path1 = self.get_selection_path(self.bed_app_data)
            self.add_log('Reading Phenotypic file...')
            pheno_path, current_path2 = self.get_selection_path(self.pheno_app_data)
            # Replace chromosome names, they need to be numbers
            chrom_mapping = self.helper.replace_with_integers(bed_path.replace('.bed', '.bim'))
            gp_df = self.gwas.start_gwas(bed_path, pheno_path, chrom_mapping, self.algorithm, self.add_log,
                                         self.test_size, self.estimators, self.gwas_result_name, True,
                                         self.genomic_predict_name)

            if gp_df is not None:
                self.add_log('Genomic Prediction done.')
                self.add_log('Genomic Prediction Plotting...')
                #self.gwas.plot_gwas(gp_df, 10000, self.algorithm, self.manhatten_plot_name, self.qq_plot_name)
                self.add_log('Done...')
                self.show_results_window(gp_df, self.algorithm, genomic_predict=True)
            else:
                self.add_log('Error, GWAS Analysis could not be started.', error=True)

        except TypeError:
            self.add_log('Please select a phenotype and genotype file. ', error=True)

    def show_results_window(self, df, algorithm, genomic_predict):

        with dpg.window(label="Results", width=975, height=600, horizontal_scrollbar=True, pos=(1000, 35)):
            dpg.add_button(label="Download Results", pos =(400, 40), callback=lambda: dpg.show_item("select_directory"))
            dpg.add_spacer(height=60)

            if genomic_predict:
                with dpg.tab_bar(label='tabbar'):
                    with dpg.tab(label="Genomic Prediction "):
                        df = df[['ID1', 'BED_ID2', 'Predicted_Value', 'Pheno_Value']]
                        df.columns = ['FID', 'IID', 'Predicted_Value', 'Pheno_Value']
                        # df = df.sort_values(by=['PValue'], ascending=True)
                        with dpg.table(label='DatasetTable2', row_background=True, borders_innerH=True,
                                       borders_outerH=True, borders_innerV=True, borders_outerV=True, tag='table_gp'):
                            for i in range(df.shape[1]):
                                dpg.add_table_column(label=df.columns[i], parent='table_gp')
                            for i in range(len(df)):
                                with dpg.table_row():
                                    for j in range(df.shape[1]):
                                        dpg.add_text(f"{df.iloc[i, j]}")

            else:

                width, height, channels, data = dpg.load_image(self.manhatten_plot_name)
                with dpg.texture_registry(show=False):
                    dpg.add_static_texture(width=width, height=height, default_value=data, tag="manhatten_tag")

                with dpg.tab_bar(label='tabbar'):
                    with dpg.tab(label="Manhatten Plot"):
                        dpg.add_image(texture_tag="manhatten_tag", tag="manhatten_image", width=950, height=400)
                    if algorithm == "FaST-LMM" or algorithm == "Linear regression":
                        width2, height2, channels2, data2 = dpg.load_image(self.qq_plot_name)
                        with dpg.texture_registry(show=False):
                            dpg.add_static_texture(width=width2, height=height2, default_value=data2, tag="qq_tag")
                        df = df.sort_values(by=['PValue'], ascending=True)
                        with dpg.tab(label="QQ-Plot"):
                            dpg.add_image(texture_tag="qq_tag", tag="qq_image", height=450, width=450)
                    else:
                        df = df.sort_values(by=['PValue'], ascending=False)

                    with dpg.tab(label="GWAS Results (Top 500)"):
                        df = df[['SNP', 'Chr', 'ChrPos', 'PValue']]
                        #df = df.sort_values(by=['PValue'], ascending=True)
                        with dpg.table(label='DatasetTable',row_background=True, borders_innerH=True,
                                       borders_outerH=True, borders_innerV=True, borders_outerV=True, tag='table_gwas'):
                            for i in range(df.shape[1]):
                                dpg.add_table_column(label=df.columns[i], parent='table_gwas')
                            for i in range(500):
                                with dpg.table_row():
                                    for j in range(df.shape[1]):
                                        dpg.add_text(f"{df.iloc[i, j]}")


    def run(self):
        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.start_dearpygui()
        dpg.destroy_context()


if __name__ == "__main__":
    main()

