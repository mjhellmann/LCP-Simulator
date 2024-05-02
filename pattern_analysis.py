from ._anvil_designer import pattern_analysisTemplate
from anvil import *
import anvil.server
import anvil.tables as tables
import anvil.tables.query as q
from anvil.tables import app_tables


class pattern_analysis(pattern_analysisTemplate):
    def __init__(self, **properties):
        # Set Form properties and Data Bindings.
        self.init_components(**properties)
        self.card_nr_PA_options.visible = False
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
        # show pop up window
        Notification("The simulation is running, it might take up to 20 seconds until the output is shown.",
                             title="Input submitted!",
                            timeout=5,
                            style="success").show()
        # calculate data for plot and csv file
        if overall_FA == None:
            filename = "nmr_DP%s_FA0-1_mol%s_s%s_p%s" % (DP, molecules, strength, pattern)
        else:
            filename = "nmr_DP%s_FA%s_mol%s_s%s_p%s" % (DP, overall_FA, molecules, strength, pattern)
        data_nmr = anvil.server.call('diads_triads', DP, overall_FA, molecules, strength, pattern)
        # create plot and show options to download
        media_obj = anvil.server.call('make_PA_plot',
                                      data_nmr,
                                      overall_FA,
                                      filename)
        self.image_1.source = media_obj
        self.download_link_plot.url = media_obj
        self.download_link_data.url = anvil.server.call('nmr_csv',
                                                        data_nmr,
                                                        filename)
        self.download_link_plot.visible = True
        self.download_link_data.visible = True
        self.image_1.scroll_into_view()


    def check_box_nr_PA_change(self, **event_args):
        """This method is called when this checkbox is checked or unchecked"""
        if self.check_box_nr_PA.checked == True:
            self.card_nr_PA_options.visible = True
            self.card_nr_PA_options.scroll_into_view()
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
