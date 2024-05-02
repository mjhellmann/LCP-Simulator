from ._anvil_designer import app_frameTemplate
from anvil import *
import anvil.server
import anvil.tables as tables
import anvil.tables.query as q
from anvil.tables import app_tables
from ..DP_histogram import DP_histogram
from ..product_profile import product_profile
from ..pattern_analysis import pattern_analysis
from ..home import home
from ..block_sizes import block_sizes

class app_frame(app_frameTemplate):
    def __init__(self, **properties):
        # Set Form properties and Data Bindings.
        self.init_components(**properties)
        self.card_home.background = "#68a39b"
        hm = home()
        self.content_panel.add_component(hm)
        
    
    def link_DP_histogram_click(self, **event_args):
        self.card_block_sizes.background = "#333333"
        self.card_DP_histogram.background = "grey"
        self.card_home.background = "#00796b"
        self.card_NMR.background = "#333333"
        self.card_product_profiles.background = "#333333"
        dphis = DP_histogram()
        self.content_panel.clear()
        self.content_panel.add_component(dphis)

    
    def link_product_profile_click(self, **event_args):
        self.card_block_sizes.background = "#333333"
        self.card_DP_histogram.background = "#333333"
        self.card_home.background = "#00796b"
        self.card_NMR.background = "#333333"
        self.card_product_profiles.background = "grey"
        proprof = product_profile()
        self.content_panel.clear()
        self.content_panel.add_component(proprof)


    def link_nmr_click(self, **event_args):
        self.card_block_sizes.background = "#333333"
        self.card_DP_histogram.background = "#333333"
        self.card_home.background = "#00796b"
        self.card_NMR.background = "grey"
        self.card_product_profiles.background = "#333333"
        nmr = pattern_analysis()
        self.content_panel.clear()
        self.content_panel.add_component(nmr)

    
    def link_block_sizes_click(self, **event_args):
        self.card_block_sizes.background = "grey"
        self.card_DP_histogram.background = "#333333"
        self.card_home.background = "#00796b"
        self.card_NMR.background = "#333333"
        self.card_product_profiles.background = "#333333"
        b_sizes = block_sizes()
        self.content_panel.clear()
        self.content_panel.add_component(b_sizes)

    def link_home_click(self, **event_args):
        self.card_block_sizes.background = "#333333"
        self.card_DP_histogram.background = "#333333"
        self.card_home.background = "#68a39b"
        self.card_NMR.background = "#333333"
        self.card_product_profiles.background = "#333333"
        hm = home()
        self.content_panel.clear()
        self.content_panel.add_component(hm)
