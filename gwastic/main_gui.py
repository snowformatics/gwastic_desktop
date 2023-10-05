import dearpygui.dearpygui as dpg
from dearpygui_ext import logger
from gwas_pipeline import GWAS
from helpers import HELPERS
#import pandas as pd
#import matplotlib.pyplot as plt
#import geneview as gv
#from math import sin

# sindatax = []
# sindatay = []
# for i in range(0, 100):
#     sindatax.append(i / 100)
#     sindatay.append(0.5 + 0.5 * sin(50 * i / 100))

# dataset2 = pd.read_csv("single_snp.csv", delimiter='\t')  # Take your df from wherever
# dataset = dataset2[['SNP', 'Chr', 'ChrPos', 'PValue']]
# dataset = dataset.head(1000)
# dataset = dataset.sort_values(by=['Chr', 'ChrPos'])
# dataset = dataset[dataset["Chr"] == 1.0]
# dataset["Chr"] = dataset["Chr"].replace(1.0, 'Chr1')
# print (dataset)
# sindatax = dataset['ChrPos'].tolist()
# sindatay = dataset['PValue'].tolist()


class GWASApp:
    def __init__(self):
        self.gwas = GWAS()
        self.helper = HELPERS()
        dpg.create_context()

        with dpg.font_registry():
            self.font = dpg.add_font("test.ttf", 20*2, tag="ttf-font")

        with dpg.theme() as self.our_theme:
            with dpg.theme_component(dpg.mvAll):
                dpg.add_theme_color(dpg.mvThemeCol_Button, (79, 143, 192, 255))

        self.vcf_file = ''
        self.pheno_file = ''

        self.vcf_app_data = None
        self.variants_app_data = None
        self.bed_app_data = None
        self.pheno_app_data = None
        self.default_path = "C:/gwas_test_data/test/"

        self.log_win = dpg.add_window(label="Log", pos=(0, 600), width=1000, height=500)
        self.logz = logger.mvLogger(self.log_win)

        # Set up GUI components and callbacks
        self.setup_gui()

    def save_callback(self):
        print("Save Clicked")

    def setup_gui(self):
        dpg.create_viewport(title='Custom Title', width=2000, height=1200)

        def print_me(sender):
            print(f"Menu Item: {sender}")
            #self.default_path = "C:/gwas_test_data/test/"


    # Menu bar
        with dpg.viewport_menu_bar():
            with dpg.menu(label="Settings"):
                dpg.add_menu_item(label="Setting 1", callback=print_me, check=True)
                dpg.add_menu_item(label="Setting 2", callback=print_me)
            dpg.add_menu_item(label="Help", callback=print_me)

        # File dialogs
        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_vcf, file_count=3, tag="file_dialog_vcf",
                             width=700, height=400, default_path=self.default_path):

            dpg.add_file_extension(".vcf", color=(255, 150, 150, 255))
            dpg.add_file_extension(".gz", color=(255, 255, 0, 255))
            #dpg.add_file_extension(".*")

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_variants, file_count=3,
                             tag="file_dialog_variants", width=700, height=400, default_path=self.default_path):
            dpg.add_file_extension(".txt", color=(255, 150, 150, 255))
            dpg.add_file_extension(".*")
            dpg.add_file_extension(".csv", color=(255, 255, 0, 255))

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_pheno, file_count=3,
                             tag="file_dialog_pheno", width=700, height=400, default_path=self.default_path):
            dpg.add_file_extension(".txt", color=(255, 150, 150, 255))
            dpg.add_file_extension(".*")
            dpg.add_file_extension(".csv", color=(255, 255, 0, 255))

        with dpg.file_dialog(directory_selector=False, show=False, callback=self.callback_bed, file_count=3,
                             tag="file_dialog_bed", width=700, height=400, default_path=self.default_path):
            dpg.add_file_extension(".bed", color=(255, 150, 150, 255))
            dpg.add_file_extension(".*")

        with dpg.window(label="GWAStic", width=1000, height=600):
            with dpg.tab_bar(label='tabbar'):
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
                    maf_input = dpg.add_input_float(label="MAF", width=150, default_value=0.05, indent=50, tag= 'tooltip_maf')
                    geno_input = dpg.add_input_float(label="Missing genotype rate", width=150, default_value=0.1, indent=50, tag= 'tooltip_missing')
                    dpg.add_spacer(height=20)
                    convert_btn = dpg.add_button(label="Convert", callback=self.convert_vcf, user_data=[maf_input, geno_input, vcf, variant_ids], indent=50)
                    dpg.bind_item_theme(convert_btn, self.our_theme)

                with dpg.tab(label='GWAS Analysis'):
                    dpg.add_text("\nStart Fast-LMM GWAS", indent=50)
                    dpg.add_spacer(height=20)
                    geno = dpg.add_button(label="Choose BED", callback=lambda: dpg.show_item("file_dialog_bed"), indent=50, tag= 'tooltip_bed')
                    pheno = dpg.add_button(label="Choose Phenotype", callback=lambda: dpg.show_item("file_dialog_pheno"), indent=50, tag= 'tooltip_pheno')
                    dpg.add_spacer(height=20)
                    #dpg.add_checkbox(label="Replace Chromosome Labels", callback=self.callback_checkbox, indent=50)
                    #dpg.add_spacer(height=20)
                    gwas_btn = dpg.add_button(label="Run GWAS", callback=self.run_gwas, user_data=[geno, pheno], indent=50)
                    dpg.bind_item_theme(gwas_btn, self.our_theme)

                with dpg.tab(label='Genomic Prediction Analysis'):
                    dpg.add_button(label="Comming soon", callback=self.retrieve_callback, user_data=[maf_input, geno_input, vcf, pheno])

            # Tooltips
            with dpg.tooltip("tooltip_vcf"):
                dpg.add_text("Choose a VCF file (.gz or .vcf).", color=[79,128,226])
            with dpg.tooltip("tooltip_variant"):
                dpg.add_text("Choose a variant file.\nVariant IDs must match with IDs in the VCF file.\n By using a variant file, you can create a subset of your VCF file.\n Highly recommended if only a subset with phenotypic data is available).", color=[79,128,226])
            with dpg.tooltip("tooltip_maf"):
                dpg.add_text("Filters out all variants with minor allele frequency below the provided threshold.", color=[79, 128, 226])
            with dpg.tooltip("tooltip_missing"):
                dpg.add_text("Filters out all variants with missing call rates exceeding the provided value to be removed.", color=[79, 128, 226])
            with dpg.tooltip("tooltip_missing"):
                dpg.add_text("Filters out all variants with missing call rates exceeding the provided value to be removed.", color=[79, 128, 226])
            with dpg.tooltip("tooltip_bed"):
                dpg.add_text("Choose a BED file.\nPrimary representation of genotype calls at biallelic variants. Must be accompanied by .bim and .fam files.\n", color=[79, 128, 226])
            with dpg.tooltip("tooltip_pheno"):
                dpg.add_text("Choose a phenotype file.\nVariant IDs must match with IDs in the .fam file.\nMust be space seperated.\nExample.\nID1 0.25\nID2 0.89\nImportant:ID's must not contain spaces", color=[79, 128, 226])



            dpg.bind_font(self.font)
            dpg.set_global_font_scale(0.6)

    def show_plot(self, df):
        width, height, channels, data = dpg.load_image("manhatten.png")
        width2, height2, channels2, data2 = dpg.load_image("qq.png")
        with dpg.texture_registry(show=False):
            dpg.add_static_texture(width=width, height=height, default_value=data, tag="manhatten_tag")
            dpg.add_static_texture(width=width2, height=height2, default_value=data2, tag="qq_tag")

        with dpg.window(label="Results", width=950, height=600, pos=[1000, 1]):
            with dpg.tab_bar(label='tabbar'):
                with dpg.tab(label="Manhatten Plot"):
                    dpg.add_image("manhatten_tag")
                with dpg.tab(label="QQ-Plot"):
                    dpg.add_image("qq_tag")
                with dpg.tab(label="GWAS Results"):
                    dataset2 = df
                    dataset = dataset2[['SNP', 'Chr', 'ChrPos', 'PValue']]
                    with dpg.table(label='DatasetTable',row_background=True, borders_innerH=True,
                                   borders_outerH=True, borders_innerV=True, borders_outerV=True):
                        for i in range(dataset.shape[1]):  # Generates the correct amount of columns
                            dpg.add_table_column(label=dataset.columns[i])  # Adds the headers
                        for i in range(500):  # Shows the first n rows
                            with dpg.table_row():
                                for j in range(dataset.shape[1]):
                                    dpg.add_text(f"{dataset.iloc[i, j]}")  # Displays the value of
                                    # each row/column combination

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

    def add_log(self, message, warn=False):
        """Adds a log message."""
        if warn:
            self.logz.log_warning(message)
        else:
            self.logz.log_info(message)

    def callback_checkbox(self, sender):
        print(f"Menu Item: {sender}")
        value = dpg.get_value(sender)
        print(value)

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

        self.add_log('Starting converting VCF to BED...')
        plink_log = self.gwas.vcf_to_bed(vcf_path, variants_path, current_path + 'test', maf, geno)
        self.add_log(plink_log)
        self.logz.log_debug('Converting done!')

    def run_gwas(self, sender, data, user_data):
        self.add_log('Reading Bed file...')
        bed_path, current_path1 = self.get_selection_path(self.bed_app_data)
        self.add_log('Reading Phenotypic file...')
        pheno_path, current_path2 = self.get_selection_path(self.pheno_app_data)
        self.add_log('Replacing chromosome names...')
        chrom_mapping = self.helper.replace_with_integers(bed_path.replace('.bed', '.bim'))
        gwas_df = self.gwas.start_gwas(bed_path, pheno_path, chrom_mapping, self.add_log)
        if gwas_df:
            self.add_log('GWAS Analysis done.')
            self.add_log('GWAS Results Plotting...')
            self.gwas.plot_gwas(gwas_df, 10000)
            self.add_log('Done...')
            self.show_plot(gwas_df)
        else:
            self.add_log('GWAS Analysis done.')

    def retrieve_callback(self, sender, data, user_data):
        pass

    def run(self):
        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.start_dearpygui()
        dpg.destroy_context()


if __name__ == "__main__":
    app = GWASApp()
    app.run()
