import dearpygui.dearpygui as dpg
from dearpygui_ext import logger
from gwas_pipeline import GWAS
from helpers import HELPERS

gwas = GWAS()
helper = HELPERS()
dpg.create_context()

with dpg.font_registry():
    font = dpg.add_font("test.ttf", 20*2, tag="ttf-font")

with dpg.theme() as our_theme:
    with dpg.theme_component(dpg.mvAll):
        dpg.add_theme_color(dpg.mvThemeCol_Button, (79, 143, 192, 255))


vcf_file = ''
pheno_file = ''

vcf_app_data = None
variants_app_data = None
bed_app_data = None
pheno_app_data = None

log_win = dpg.add_window(label="Log", pos=(0, 600), width=1000, height=500)
logz = logger.mvLogger(log_win)


# Callbacks
def get_selection_path(app_data):
    """Extract path from the app_data dictionary selections key."""
    current_path = app_data['current_path'] + '/'
    k = app_data['selections']
    for key, value in k.items():
        file_path = value

    return file_path, current_path


def add_log(message):
    """Adds a log message."""
    logz.log_info(message)


def callback_vcf(sender, app_data):
    """Get vcf file path selected from the user."""
    global vcf_app_data  # Use the global variable to store the value
    vcf_app_data = app_data
    vcf_path, current_path = get_selection_path(vcf_app_data)
    add_log('VCF Selected: ' + vcf_path)


def callback_bed(sender, app_data):
    """Get vcf file path selected from the user."""
    global bed_app_data  # Use the global variable to store the value
    bed_app_data = app_data
    bed_path, current_path = get_selection_path(bed_app_data)
    add_log('BED file Selected: ' + bed_path)


def callback_variants(sender, app_data):
    """Get variant file path selected from the user."""
    global variants_app_data  # Use the global variable to store the value
    variants_app_data = app_data
    variants_path, current_path = get_selection_path(variants_app_data)
    add_log('Variant File Selected: ' + variants_path)


def callback_pheno(sender, app_data):
    """Get phenotype file path selected from the user."""
    global pheno_app_data  # Use the global variable to store the value
    pheno_app_data = app_data
    pheno_path, current_path = get_selection_path(pheno_app_data)
    add_log('Pheno File Selected: ' + pheno_path)



def retrieve_callback(sender, callback):
    pass
    #dpg.show_logger()
    #dpg.log_info(dpg.get_value("MAF##inputfloat"))


# Get functions
def convert_vcf(sender, data, user_data):
    """Converts vcf to bed file using PLINK."""
    maf = str(dpg.get_value(user_data[0]))
    geno = str(dpg.get_value(user_data[1]))
    variant_ids = variants_app_data
    vcf_path, current_path = get_selection_path(vcf_app_data)

    if variant_ids == None:
        variants_path = None
    else:
        variants_path, current_path = get_selection_path(variants_app_data)

    add_log('Staring converting VCF to BED...')
    plink_log = gwas.vcf_to_bed(vcf_path, variants_path, current_path + 'test', maf, geno)
    add_log(plink_log)
    logz.log_debug('Converting done!')


def run_gwas(sender, data, user_data):
    bed_path, current_path1 = get_selection_path(bed_app_data)
    pheno_path, current_path2 = get_selection_path(pheno_app_data)
    helper.replace_with_integers(bed_path.replace('.bed', '.bim'))
    gwas.start_gwas(bed_path, pheno_path)
    #print(bed_path, current_path1, pheno_path, current_path2)


def get_data_gp(sender, data, user_data):
    pass


# File diaologs
with dpg.file_dialog(directory_selector=False, show=False, callback=callback_vcf, file_count=3, tag="file_dialog_vcf",
                     width=700 ,height=400, default_path="E:/GWAS_data/test/"):
    dpg.add_file_extension(".vcf", color=(255, 150, 150, 255))
    dpg.add_file_extension(".gz", color=(255, 255, 0, 255))
    dpg.add_file_extension(".*")

with dpg.file_dialog(directory_selector=False, show=False, callback=callback_variants, file_count=3,
                     tag="file_dialog_variants", width=700 ,height=400, default_path="E:/GWAS_data/test/"):
    dpg.add_file_extension(".txt", color=(255, 150, 150, 255))
    dpg.add_file_extension(".*")
    dpg.add_file_extension(".csv", color=(255, 255, 0, 255))

with dpg.file_dialog(directory_selector=False, show=False, callback=callback_pheno, file_count=3,
                     tag="file_dialog_pheno", width=700 ,height=400, default_path="E:/GWAS_data/test/"):
    dpg.add_file_extension(".txt", color=(255, 150, 150, 255))
    dpg.add_file_extension(".*")
    dpg.add_file_extension(".csv", color=(255, 255, 0, 255))

with dpg.file_dialog(directory_selector=False, show=False, callback=callback_bed, file_count=3,
                     tag="file_dialog_bed", width=700 ,height=400,default_path="E:/GWAS_data/test/"):
    dpg.add_file_extension(".bed", color=(255, 150, 150, 255))
    dpg.add_file_extension(".*")


with dpg.window(label="GWAStic", width=1000, height=600):
    with dpg.tab_bar(label='tabbar'):
        with dpg.tab(label='Convert VCF'):
            dpg.add_text("\nConvert a VCF file into BED file format and apply MAF\nor missing genotype filter.", indent=50)
            dpg.add_spacer(height=20)
            dpg.add_text("Select files:", indent=50)
            vcf = dpg.add_button(label="Choose VCF", callback=lambda: dpg.show_item("file_dialog_vcf"), indent=50)
            variant_ids = dpg.add_button(label="Choose Variants",
                                         callback=lambda: dpg.show_item("file_dialog_variants"), indent=50)
            dpg.add_spacer(height=20)
            dpg.add_text("Apply filters:", indent=50)
            maf_input = dpg.add_input_float(label="MAF", width=150, default_value=0.05, indent=50)
            geno_input = dpg.add_input_float(label="Missing genotype rate", width=150, default_value=0.1, indent=50)
            dpg.add_spacer(height=20)
            convert_btn = dpg.add_button(label="Convert", callback=convert_vcf, user_data=[maf_input, geno_input, vcf, variant_ids], indent=50)
            dpg.bind_item_theme(convert_btn, our_theme)

        with dpg.tab(label='GWAS Analysis'):
            dpg.add_text("\nStart Fast-LMM GWAS.",
                         indent=50)
            dpg.add_spacer(height=20)
            geno = dpg.add_button(label="Choose BED", callback=lambda: dpg.show_item("file_dialog_bed"), indent=50)
            pheno = dpg.add_button(label="Choose Phenotype", callback=lambda: dpg.show_item("file_dialog_pheno"), indent=50)
            dpg.add_spacer(height=20)
            gwas_btn = dpg.add_button(label="Run GWAS", callback=run_gwas, user_data=[geno, pheno], indent=50)
            dpg.bind_item_theme(gwas_btn, our_theme)

        with dpg.tab(label='Genomic Prediction Analysis'):
            dpg.add_button(label="Run Genomic Prediction", callback=retrieve_callback, user_data=[maf_input, geno_input, vcf, pheno])

    dpg.bind_font(font)
    dpg.set_global_font_scale(0.6)

with dpg.window(label="Manhatten Plot", width=800, height=400, pos=[1000,1]):
    dpg.add_simple_plot(label="Simpleplot1", default_value=(0.3, 0.9, 0.5, 0.3), height=300)



dpg.create_viewport(title='Custom Title', width=2000, height=1200)
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()