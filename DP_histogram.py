from ._anvil_designer import DP_histogramTemplate
from anvil import *
import anvil.server
import anvil.tables as tables
import anvil.tables.query as q
from anvil.tables import app_tables


class DP_histogram(DP_histogramTemplate):
    def __init__(self, **properties):
        # Set Form properties and Data Bindings.
        self.init_components(**properties)
        self.card_nr_PA_options.visible = False
        self.strength_box.visible = False
        self.FA_pattern_box.visible = False
        self.card_enzyme_options.visible = False
        self.card_weight_monomers.visible = False
        self.card_nano2_options.visible = False
        self.download_link_plot.visible = False
        self.download_link_data.visible = False
    
    
    def submit_button_click(self, **event_args):
        # check if radio buttons were clicked
        if (self.radio_button_enzyme.selected == False and
            self.radio_button_nano2.selected == False):
                Notification("Choose cleavage type: 'enzymatic' or 'NaNO2'",
                             title="Input missing!",
                            timeout=5,
                            style="danger").show()
                return
        if (self.radio_button_weight_fraction.selected == False and
            self.radio_button_molar_fraction.selected == False):
                Notification("Choose quantification of products: 'molar fraction' or 'weight fraction'",
                             title="Input missing!",
                            timeout=5,
                            style="danger").show()
                return
        # save input parameters
        if self.DP_box.text == "":
            DP = 100
        else:
            try:
                DP = int(self.DP_box.text)
                if DP < 50:
                    Notification("The specified units per molecule are rather low, the generated output might not be meaningful (see home, section 'Polymer library')",
                             title="Potentially unsuitable input!",
                            timeout=5,
                            style="warning").show()
                elif DP > 10000:
                    Notification("The specified units per molecule are rather high, the algorithm might need longer than allowed on the server (see home, section 'Polymer library')",
                             title="Potentially unsuitable input!",
                            timeout=5,
                            style="warning").show()
            except:
                Notification("Only integers are accepted as input for units per molecule, e.g. 100",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
        if self.overall_FA_box.text == "":
            overall_FA = 0.32
        else:
            try:
                overall_FA = float(self.overall_FA_box.text)
                if overall_FA < 0 or overall_FA > 1:
                    Notification("Only numbers from 0-1 are accepted as input for average fraction of unit A, e.g. 0.32",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                    return
            except:
                Notification("Only numbers from 0-1 are accepted as input for average fraction of unit A, e.g. 0.32",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
        if self.molecules_box.text == "":
            molecules = 100
        else:
            try:
                molecules = int(self.molecules_box.text)
                if molecules < 50:
                    Notification("The specified number of molecules is rather low, the generated output might not be meaningful (see home, section 'Polymer library')",
                             title="Potentially unsuitable input!",
                            timeout=5,
                            style="warning").show()
                elif molecules > 10000:
                    Notification("The specified number of molecules is rather high, the algorithm might need longer than allowed on the server (see home, section 'Polymer library')",
                             title="Potentially unsuitable input!",
                            timeout=5,
                            style="warning").show()
            except:
                Notification("Only integers are accepted as input for number of molecules, e.g. 100",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
        if self.pattern_box.text == "" or self.card_nr_PA_options.visible == False:
            pattern = 3
        else:
            try:
                pattern = int(self.pattern_box.text)
                if pattern < 2:
                    Notification("A frequency of overrepresented A-units < 2 is not sensible",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                    return
                elif pattern > 100:
                    Notification("The frequency of overrepresented A-units is rather high, the generated output might not be meaningful (see home, section 'Polymer library')",
                             title="Potentially unsuitable input!",
                            timeout=5,
                            style="warning").show()
            except:
                Notification("Only integers are accepted as input for the frequency of overrepresented A-units, e.g. 3",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
        if self.FA_pattern_box.text == "" or self.card_nr_PA_options.visible == False:
            FA_pattern = overall_FA
        else:
            try:
                FA_pattern = float(self.FA_pattern_box.text)
                if FA_pattern < 0 or FA_pattern > 1:
                    Notification("Only numbers from 0-1 are accepted as input for average fraction of unit A within the units with overrepresentation of A, e.g. 0.32",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                    return
                elif FA_pattern < overall_FA:
                    Notification("The specified average fraction of unit A within the units with overrepresentation of A is smaller than the overall average fraction of unit A",
                             title="Potentially unsuitable input!",
                            timeout=5,
                            style="warning").show()
            except:
                Notification("Only numbers from 0-1 are accepted as input for average fraction of unit A within the units with overrepresentation of A, e.g. 0.32",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
        if self.strength_box.text == "" or self.card_nr_PA_options.visible == False:
            strength = 0
        else:
            try:
                strength = float(self.strength_box.text)
                if strength < 0 or strength > 1:
                    Notification("Only numbers from 0-1 are accepted as input for strength of overrepresentation of A-units, e.g. 0.5",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                    return
            except:
                Notification("Only numbers from 0-1 are accepted as input for strength of overrepresentation of A-units, e.g. 0.5",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
        if self.minus_specificity_box.text == "":
            minus_specificity = "_BA"
        else:
            minus_specificity = self.minus_specificity_box.text
            if len(minus_specificity) != 3:
                Notification("The input for minus subsite specificity must consist of three characters",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
            else:
                for char in minus_specificity:
                    if char not in "ABX_":
                        Notification("Valid characters to specify the minus subsite specificity are only A, B, X and _",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                        return
        if self.plus_specificity_box.text == "":
            plus_specificity = "XX_"
        else:
            plus_specificity = self.plus_specificity_box.text
            if len(plus_specificity) != 3:
                Notification("The input for plus subsite specificity must consist of three characters",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
            else:
                for char in plus_specificity:
                    if char not in "ABX_":
                        Notification("Valid characters to specify the plus subsite specificity are only A, B, X and _",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                        return
        if self.efficiency_box.text == "":
            efficiency = 1
        else:
            try:
                efficiency = float(self.efficiency_box.text)
                if efficiency < 0 or efficiency > 1:
                    Notification("Only numbers from 0-1 are accepted as input for efficiency, e.g. 0.5",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                    return
            except:
                Notification("Only numbers from 0-1 are accepted as input for efficiency, e.g. 0.5",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
        if self.cuts_box.text == "":
            cuts = 5
        else:
            try:
                cuts = int(self.cuts_box.text)
                if cuts < 1:
                    Notification("The specified number of cleavages per chitosan molecule < 1 is not sensible",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                    return
                elif cuts > 1000:
                    Notification("The specified number of cleavages per chitosan molecule is rather high, the generated output might not be meaningful",
                             title="Potentially unsuitable input!",
                            timeout=5,
                            style="warning").show()
            except:
                Notification("Only integers are accepted as input for number of cleavages per chitosan molecule, e.g. 5",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
        if self.mass_A_box.text == "":
            mass_A = 221.097
        else:
            try:
                mass_A = float(self.mass_A_box.text)
                if mass_A < 0:
                    Notification("The specified mass of unit A of < 0 is not sensible",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                    return
                elif mass_A > 10000:
                    Notification("The specified mass of unit A is rather high, please check for accuracy",
                             title="Potentially unsuitable input!",
                            timeout=5,
                            style="warning").show()
            except:
                Notification("Only numbers are accepted as input for mass of unit A, e.g. 221.097",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
        if self.mass_B_box.text == "":
            mass_B = 179.087
        else:
            try:
                mass_B = float(self.mass_B_box.text)
                if mass_B < 0:
                    Notification("The specified mass of unit B of < 0 is not sensible",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                    return
                elif mass_B > 10000:
                    Notification("The specified mass of unit B is rather high, please check for accuracy",
                             title="Potentially unsuitable input!",
                            timeout=5,
                            style="warning").show()
            except:
                Notification("Only numbers are accepted as input for mass of unit B, e.g. 179.087",
                             title="Non-accepted input!",
                            timeout=5,
                            style="danger").show()
                return
        # show pop up window
        Notification("The simulation is running, it might take up to 20 seconds until the output is shown.",
                             title="Input submitted!",
                            timeout=5,
                            style="success").show()
        # calculate data for plot and csv file
        if self.radio_button_molar_fraction.selected == True:
            y_label = "molar fraction"
            if self.radio_button_enzyme.selected == True:
                if self.radio_button_FAp.selected == True:
                    filename = "hist_enz_molar_DP%s_FA%s_mol%s_FAp%s_p%s_-%s_+%s_e%s" % (DP, overall_FA,
                                                                                            molecules, FA_pattern,
                                                                                            pattern, minus_specificity,
                                                                                            plus_specificity, efficiency)
                    data_DP_histogram = anvil.server.call('histogram_DP_FAp_molar_enzyme',
                                                    DP, overall_FA,
                                                    pattern, FA_pattern,
                                                    molecules,
                                                    minus_specificity,
                                                    plus_specificity,
                                                    efficiency)
                else:
                    filename = "hist_enz_molar_DP%s_FA%s_mol%s_s%s_p%s_-%s_+%s_e%s" % (DP, overall_FA,
                                                                                            molecules, strength,
                                                                                            pattern, minus_specificity,
                                                                                            plus_specificity, efficiency)
                    data_DP_histogram = anvil.server.call('histogram_DP_strength_molar_enzyme',
                                                    DP, overall_FA,
                                                    pattern, strength,
                                                    molecules,
                                                    minus_specificity,
                                                    plus_specificity,
                                                    efficiency)
            elif self.radio_button_nano2.selected == True:
                if self.radio_button_FAp.selected == True:
                    filename = "hist_nano2_molar_DP%s_FA%s_mol%s_FAp%s_p%s_c%s" % (DP, overall_FA,
                                                                                            molecules, FA_pattern,
                                                                                            pattern, cuts)
                    data_DP_histogram = anvil.server.call('histogram_DP_FAp_molar_nano2',
                                                    DP, overall_FA,
                                                    pattern, FA_pattern,
                                                    molecules, cuts)
                else:
                    filename = "hist_nano2_molar_DP%s_FA%s_mol%s_s%s_p%s_c%s" % (DP, overall_FA,
                                                                                            molecules, strength,
                                                                                            pattern, cuts)
                    data_DP_histogram = anvil.server.call('histogram_DP_strength_molar_nano2',
                                                    DP, overall_FA,
                                                    pattern, strength,
                                                    molecules, cuts)
        elif self.radio_button_weight_fraction.selected == True:
            y_label = "weight fraction"
            if self.radio_button_enzyme.selected == True:
                if self.radio_button_FAp.selected == True:
                    filename = "hist_enz_weight_DP%s_FA%s_mol%s_FAp%s_p%s_-%s_+%s_e%s" % (DP, overall_FA,
                                                                                            molecules, FA_pattern,
                                                                                            pattern, minus_specificity,
                                                                                            plus_specificity, efficiency)
                    data_DP_histogram = anvil.server.call('histogram_DP_FAp_weight_enzyme',
                                                    DP, overall_FA,
                                                    pattern, FA_pattern,
                                                    molecules,
                                                    minus_specificity,
                                                    plus_specificity,
                                                    efficiency,
                                                    mass_A, mass_B)
                else:
                    filename = "hist_enz_weight_DP%s_FA%s_mol%s_s%s_p%s_-%s_+%s_e%s" % (DP, overall_FA,
                                                                                            molecules, strength,
                                                                                            pattern, minus_specificity,
                                                                                            plus_specificity, efficiency)
                    data_DP_histogram = anvil.server.call('histogram_DP_strength_weight_enzyme',
                                                    DP, overall_FA,
                                                    pattern, strength,
                                                    molecules,
                                                    minus_specificity,
                                                    plus_specificity,
                                                    efficiency,
                                                    mass_A, mass_B)
            elif self.radio_button_nano2.selected == True:
                if self.radio_button_FAp.selected == True:
                    filename = "hist_nano2_weight_DP%s_FA%s_mol%s_FAp%s_p%s_c%s" % (DP, overall_FA,
                                                                                            molecules, FA_pattern,
                                                                                            pattern, cuts)
                    data_DP_histogram = anvil.server.call('histogram_DP_FAp_weight_nano2',
                                                    DP, overall_FA,
                                                    pattern, FA_pattern,
                                                    molecules, cuts,
                                                    mass_A, mass_B)
                else:
                    filename = "hist_nano2_weight_DP%s_FA%s_mol%s_s%s_p%s_c%s" % (DP, overall_FA,
                                                                                            molecules, strength,
                                                                                            pattern, cuts)
                    data_DP_histogram = anvil.server.call('histogram_DP_strength_weight_nano2',
                                                    DP, overall_FA,
                                                    pattern, strength,
                                                    molecules, cuts,
                                                    mass_A, mass_B)     
        # create plot and show options to download
        media_obj = anvil.server.call('make_histogram',
                                      data_DP_histogram,
                                      y_label, filename)
        self.image_1.source = media_obj
        self.download_link_plot.url = media_obj
        self.download_link_data.url = anvil.server.call('histogram_csv',
                                                        data_DP_histogram,
                                                        y_label, filename)
        self.download_link_plot.visible = True
        self.download_link_data.visible = True
        self.image_1.scroll_into_view()
    

    def check_box_nr_PA_change(self, **event_args):
        """This method is called when this checkbox is checked or unchecked"""
        if self.check_box_nr_PA.checked == True:
            self.card_nr_PA_options.visible = True
            #self.card_nr_PA_options.scroll_into_view()
        else:
            self.card_nr_PA_options.visible = False

    
    def radio_button_FAp_clicked(self, **event_args):
        """This method is called when this radio button is selected"""
        self.strength_box.visible = False
        self.FA_pattern_box.visible = True

    
    def radio_button_strength_clicked(self, **event_args):
        """This method is called when this radio button is selected"""
        self.strength_box.visible = True
        self.FA_pattern_box.visible = False


    def radio_button_enzyme_clicked(self, **event_args):
        """This method is called when this radio button is selected"""
        self.card_enzyme_options.visible = True
        self.card_nano2_options.visible = False
        #self.card_enzyme_options.scroll_into_view()

    def radio_button_nano2_clicked(self, **event_args):
        """This method is called when this radio button is selected"""
        self.card_enzyme_options.visible = False
        self.card_nano2_options.visible = True
        #self.card_nano2_options.scroll_into_view()

    def radio_button_molar_fraction_clicked(self, **event_args):
        """This method is called when this radio button is selected"""
        self.card_weight_monomers.visible = False

    def radio_button_weight_fraction_clicked(self, **event_args):
        """This method is called when this radio button is selected"""
        self.card_weight_monomers.visible = True
        #self.card_weight_monomers.scroll_into_view()

    def DP_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def overall_FA_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def molecules_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def pattern_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def FA_pattern_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def strength_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def minus_specificity_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def plus_specificity_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def efficiency_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def cuts_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def mass_A_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def mass_B_box_pressed_enter(self, **event_args):
        self.submit_button_click()
