from ._anvil_designer import block_sizesTemplate
from anvil import *
import anvil.server
import anvil.tables as tables
import anvil.tables.query as q
from anvil.tables import app_tables


class block_sizes(block_sizesTemplate):
    def __init__(self, **properties):
        # Set Form properties and Data Bindings.
        self.init_components(**properties)
        self.card_nr_PA_options.visible = False
        self.eff_hydrolysis_label.visible = False
        self.efficiency_box.visible = False
        self.card_lim_ms_options.visible = False
        self.download_link_plot.visible = False
        self.download_link_data.visible = False
        self.label_7.visible = False
        self.box_overall_FA.visible = False

    def submit_button_click(self, **event_args):
        # check if radio buttons were clicked
        if (self.radio_button_plot.selected == False and
            self.radio_button_single.selected == False):
                Notification("Choose whole range for average fraction of unit A or single specific value",
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
        if self.box_overall_FA.visible == False:
            overall_FA = None     
            if DP * molecules > 20000:
                Notification("The specified combination of units per molecule and number of molecules results in a high number of simulated monomers (units per molecule x number of molecules > 20,000), the algorithm might need longer than allowed on the server (see home, section 'Polymer library')",
                    title="Potentially unsuitable input!",
                    timeout=5,
                    style="warning").show()
        elif self.box_overall_FA.text == "":
            overall_FA = 0.32
        else:
            try:
                overall_FA = float(self.box_overall_FA.text)
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
        if self.efficiency_box.visible == False:
            efficiency = 1
        elif self.efficiency_box.text == "":
            efficiency = 0.5
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
        if self.card_lim_ms_options.visible == False:
            ion_eff = False
            DP_cutoff = False
            A_cutoff = False
        else:
            if self.check_box_ion_eff.checked == True:
                ion_eff = True
            else:
                ion_eff = False
            if self.DP_cutoff_box.text == "":
                DP_cutoff = 14
            else:
                try:
                    DP_cutoff = int(self.DP_cutoff_box.text)
                    if DP_cutoff < 4:
                        Notification("The specified DP limit up to which products are considered is rather low, the generated output might not be meaningful",
                                title="Potentially unsuitable input!",
                                timeout=5,
                                style="warning").show()
                except:
                    Notification("Only integers are accepted as input for the DP limit up to which products are considered, e.g. 14",
                                title="Non-accepted input!",
                                timeout=5,
                                style="danger").show()
                    return  
            if self.A_cutoff_box.text == "":
                A_cutoff = 8
            else:
                try:
                    A_cutoff = int(self.A_cutoff_box.text)
                    if A_cutoff < 3:
                        Notification("The specified limit for A-units up to which products are considered is rather low, the generated output might not be meaningful",
                                title="Potentially unsuitable input!",
                                timeout=5,
                                style="warning").show()
                except:
                    Notification("Only integers are accepted as input for the limit for A-units up to which products are considered, e.g. 8",
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
        if overall_FA == None:
            filename = "blocks_DP%s_FA0-1_mol%s_s%s_p%s_e%s_Ac%s_DPc%s_ie%s" % (
                DP, molecules, strength, pattern,
                efficiency, A_cutoff, DP_cutoff, ion_eff)
        else:
            filename = "blocks_DP%s_FA%s_mol%s_s%s_p%s_e%s_Ac%s_DPc%s_ie%s" % (
                DP, overall_FA, molecules, strength, pattern,
                efficiency, A_cutoff, DP_cutoff, ion_eff)
        data_blocks = anvil.server.call('block_sizes', DP, overall_FA, molecules, strength, pattern,
                                        efficiency, A_cutoff, DP_cutoff, ion_eff)
        # create plot and show options to download
        media_obj = anvil.server.call('make_blocks_plot',
                                      data_blocks,
                                      overall_FA,
                                      filename)
        self.image_1.source = media_obj
        self.download_link_plot.url = media_obj
        self.download_link_data.url = anvil.server.call('blocks_csv',
                                                        data_blocks,
                                                        filename)
        self.download_link_plot.visible = True
        self.download_link_data.visible = True
        self.image_1.scroll_into_view()


    def check_box_nr_PA_change(self, **event_args):
        """This method is called when this checkbox is checked or unchecked"""
        if self.check_box_nr_PA.checked == True:
            self.card_nr_PA_options.visible = True
        else:
            self.card_nr_PA_options.visible = False

    def radio_button_plot_clicked(self, **event_args):
        """This method is called when this radio button is selected"""
        self.label_7.visible = False
        self.box_overall_FA.visible = False

    def radio_button_single_clicked(self, **event_args):
        """This method is called when this radio button is selected"""
        self.label_7.visible = True
        self.box_overall_FA.visible = True

    def check_box_incomp_lys_change(self, **event_args):
        if self.check_box_incomp_lys.checked == True:
            self.eff_hydrolysis_label.visible = True
            self.efficiency_box.visible = True
        else:
            self.eff_hydrolysis_label.visible = False
            self.efficiency_box.visible = False

    def check_box_lim_ms_change(self, **event_args):
        if self.check_box_lim_ms.checked == True:
            self.card_lim_ms_options.visible = True
        else:
            self.card_lim_ms_options.visible = False

    def DP_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def molecules_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def box_overall_FA_pressed_enter(self, **event_args):
        self.submit_button_click()
    
    def pattern_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def strength_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def efficiency_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def DP_cutoff_box_pressed_enter(self, **event_args):
        self.submit_button_click()

    def A_cutoff_box_pressed_enter(self, **event_args):
        self.submit_button_click()
