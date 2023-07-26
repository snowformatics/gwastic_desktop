import dearpygui.dearpygui as dpg
from dearpygui_ext import logger
from run_gwas import GWAS

gwas = GWAS()

dpg.create_context()

with dpg.font_registry():
    font = dpg.add_font("test.ttf", 20*2, tag="ttf-font")

with dpg.theme() as our_theme:
    with dpg.theme_component(dpg.mvAll):
        dpg.add_theme_color(dpg.mvThemeCol_Button, (79, 143, 192, 255))


vcf_file = ''
pheno_file = ''

vcf_app_data = None

log_win = dpg.add_window(label="Log", pos=(0, 600), width=1000, height=300)
logz = logger.mvLogger(log_win)


# Callbacks
def get_selection_path(app_data):
    k = app_data['selections']
    for key, value in k.items():
        return value


def add_log(message):
    logz.log_info(message)


def callback_vcf(sender, app_data):
    global vcf_app_data  # Use the global variable to store the value
    vcf_app_data = app_data
    vcf_path = get_selection_path(vcf_app_data)
    add_log('VCF Selected: ' + vcf_path)


def callback_pheno(sender, app_data):
    add_log('Phenotype File Selected: ')
    #print("Sender: ", sender)
    #print("App Data: ", app_data)


def retrieve_callback(sender, callback):
    pass
    #dpg.show_logger()
    #dpg.log_info(dpg.get_value("MAF##inputfloat"))


# Get functions
def convert_vcf(sender, data, user_data):

    maf = str(dpg.get_value(user_data[0]))
    geno = str(dpg.get_value(user_data[1]))
    vcf_path = get_selection_path(vcf_app_data)
    add_log('Staring converting VCF to BED...')
    gwas.vcf_to_bed(vcf_path, None, 'test', maf, geno)
    #print (maf, geno, value)



def get_data_gwas(sender, data, user_data):
    pass
    #print(dpg.get_value(user_data[0]), dpg.get_value(user_data[1]))


def get_data_gp(sender, data, user_data):
    pass


# File diaologs
with dpg.file_dialog(directory_selector=False, show=False, callback=callback_vcf, file_count=3, tag="file_dialog_vcf",
                     width=700 ,height=400, default_path="E:/GWAS_data/test/"):
    dpg.add_file_extension(".vcf", color=(255, 150, 150, 255))
    dpg.add_file_extension(".gz", color=(255, 255, 0, 255))
    dpg.add_file_extension(".*")


with dpg.file_dialog(directory_selector=False, show=False, callback=callback_pheno, file_count=3, tag="file_dialog_pheno", width=700 ,height=400):
    dpg.add_file_extension(".*")
    dpg.add_file_extension(".txt", color=(255, 150, 150, 255))
    dpg.add_file_extension(".csv", color=(255, 255, 0, 255))

with dpg.file_dialog(directory_selector=False, show=False, callback=callback_vcf, file_count=3, tag="file_dialog_bed", width=700 ,height=400):
    dpg.add_file_extension(".*")
    dpg.add_file_extension(".bed", color=(255, 150, 150, 255))


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
            convert_btn = dpg.add_button(label="Convert", callback=convert_vcf, user_data=[maf_input, geno_input, vcf], indent=50)
            dpg.bind_item_theme(convert_btn, our_theme)

        with dpg.tab(label='GWAS Analysis'):
            geno = dpg.add_button(label="Choose BED", callback=lambda: dpg.show_item("file_dialog_bed"))
            pheno = dpg.add_button(label="Choose Phenotype", callback=lambda: dpg.show_item("file_dialog_pheno"))
            dpg.add_button(label="Run GWAS", callback=get_data_gwas, user_data=[geno, pheno])

        with dpg.tab(label='Genomic Prediction Analysis'):
            dpg.add_button(label="Run Genomic Prediction", callback=retrieve_callback, user_data=[maf_input, geno_input, vcf, pheno])

    dpg.bind_font(font)
    dpg.set_global_font_scale(0.6)

with dpg.window(label="Manhatten Plot", width=800, height=400, pos=[1000,1]):
    dpg.add_simple_plot(label="Simpleplot1", default_value=(0.3, 0.9, 0.5, 0.3), height=300)





dpg.create_viewport(title='Custom Title', width=2000, height=1000)
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()