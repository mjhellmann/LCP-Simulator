from ._anvil_designer import homeTemplate
from anvil import *
import anvil.server
import anvil.tables as tables
import anvil.tables.query as q
from anvil.tables import app_tables


class home(homeTemplate):
    def __init__(self, **properties):
        # Set Form properties and Data Bindings.
        self.init_components(**properties)
        self.library_text.visible = False
        self.library_less_button.visible = False
        self.histogram_text.visible = False
        self.histogram_less_button.visible = False
        self.profile_text.visible = False
        self.profile_less_button.visible = False
        self.nmr_text.visible = False
        self.nmr_less_button.visible = False
        self.blocks_text.visible = False
        self.blocks_less_button.visible = False

    def library_more_button_click(self, **event_args):
        self.library_more_button.visible = False
        self.library_less_button.visible = True
        self.library_text.visible = True

    def library_less_button_click(self, **event_args):
        self.library_less_button.visible = False
        self.library_more_button.visible = True
        self.library_text.visible = False

    def histogram_more_button_click(self, **event_args):
        self.histogram_more_button.visible = False
        self.histogram_less_button.visible = True
        self.histogram_text.visible = True

    def histogram_less_button_click(self, **event_args):
        self.histogram_less_button.visible = False
        self.histogram_more_button.visible = True
        self.histogram_text.visible = False

    def profile_more_button_click(self, **event_args):
        self.profile_more_button.visible = False
        self.profile_less_button.visible = True
        self.profile_text.visible = True

    def profile_less_button_click(self, **event_args):
        self.profile_less_button.visible = False
        self.profile_more_button.visible = True
        self.profile_text.visible = False

    def nmr_more_button_click(self, **event_args):
        self.nmr_more_button.visible = False
        self.nmr_less_button.visible = True
        self.nmr_text.visible = True

    def nmr_less_button_click(self, **event_args):
        self.nmr_less_button.visible = False
        self.nmr_more_button.visible = True
        self.nmr_text.visible = False

    def blocks_more_button_click(self, **event_args):
        self.blocks_more_button.visible = False
        self.blocks_less_button.visible = True
        self.blocks_text.visible = True

    def blocks_less_button_click(self, **event_args):
        self.blocks_less_button.visible = False
        self.blocks_more_button.visible = True
        self.blocks_text.visible = False
